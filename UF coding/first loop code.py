#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:18:12 2023

@author: xiaomeihai
"""

import geopandas as gpd
import matplotlib.pyplot as plt 
#Used for dealing with our dataframe
import pandas as pd
#Using for finding lat and long coords X nautical miles out
import geopy 
import geopy.distance
from geopy.distance import geodesic
#Used for drawing circles
from numpy import cos,sin,arccos
import numpy as np
#Used for polygon shenanigans
import shapely
from shapely.geometry import Polygon
from shapely.geometry import shape
from shapely.ops import unary_union
#Used for visualizing
import geopandas as gpd
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
# Two packages I have added to interpolate the gaps between observed hurricane
from statsmodels.nonparametric.smoothers_lowess import lowess  # version 0.13.2
from scipy.interpolate import interp1d # version 1.8.1
import warnings
import math
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as sm

## create and load Function -- Hurricance Shapefile
def parametric_circle(t,xc,yc,R):
    x = xc + R*cos(t)
    y = yc + R*sin(t)
    return x,y

def inv_parametric_circle(x,xc,R):
    t = np.arccos(((x-xc)/R))
    return t

def hurricane_shapefile(df, storm_id, wind_speed = 34, quadrant_points = 30):
    '''
    Parameters
    ----------
    df : pandas dataframe
        Hurdat 2 data.
    storm_id : string
        storm ID from hurdat.
    wind_speed : int, possible values are 34, 50, and 64
        DESCRIPTION. The default is 34.
    quadrant_points : INT, # of points in each arc
        DESCRIPTION. The default is 30.

    Returns
    -------
    List of non-empty polygons for the given hurricane
    
    '''
    #Empty list to store our shapes
    storm_shapes = []
    #Bearing for our wind directions
    bearings = [0, 90, 180, 270]
    
    #Subset for our storm
    sub_df = df[df['id'] == storm_id]
    
    for i in range(len(sub_df)):
        #Get the center of our storm
        latitude = sub_df['latitude'].iloc[i]
        longitude = sub_df['longitude'].iloc[i]
        center = geopy.Point(latitude, longitude)
        
        wind_columns = [f'{wind_speed}kt_wind_radii_NE_quad', 
                        f'{wind_speed}kt_wind_radii_SE_quad',
                        f'{wind_speed}kt_wind_radii_SW_qud', 
                        f'{wind_speed}kt_wind_radii_NW_qud']
        
        #Empty list to store the quadrant shapefiles
        quadrant_shapes = []
        for j in range(4):
            wind_miles = sub_df[wind_columns[j]].iloc[i]
            
            if wind_miles == 0:
                continue
            
            #Get the two points that make up the ends of our arc
            p1 = geodesic(nautical = wind_miles).destination(point=geopy.Point(latitude, longitude), bearing= bearings[j])
            p2 = geodesic(nautical = wind_miles).destination(point=geopy.Point(latitude, longitude), bearing= (bearings[j] + 90))
            
            #Treat our radius different depending on the quadrant
            if j == 0:
                radius = (p1.latitude - center.latitude) 
            if j ==1:
                radius = (p2.latitude - center.latitude) * -1
            elif j==2:
                radius = (p1.latitude - center.latitude) 
            elif j==3:
                radius = (p2.latitude - center.latitude) * -1
            
            #Convert our start and end points into tuples
            start_point = tuple([p1.latitude, p1.longitude])
            end_point = tuple([p2.latitude, p2.longitude])
            #Draw our semi-circle
            start_t = inv_parametric_circle(start_point[0], latitude, radius)
            end_t = inv_parametric_circle(end_point[0], latitude, radius)
            arc_T = np.linspace(start_t, end_t, quadrant_points)
            #Get the points from our circle
            X,Y = parametric_circle(arc_T, latitude, longitude, radius)
            points  = list(tuple(zip(Y, X)))
            #Insert the center of the storm as a first point so we can make the polygon
            points.insert(0, (longitude, latitude))
            #Create our polygon
            sub_poly = shapely.geometry.shape({'type': "Polygon", "coordinates": [points]})
            #Append to our list
            quadrant_shapes.append(sub_poly)
            
        #Union together each quadrant to form a circle-ish
        storm_shapes.append(unary_union(quadrant_shapes))
    #Remove any records that had 0 for all wind radii
    storm_shapes = [i for i in storm_shapes if not i.is_empty]
    
    return(storm_shapes)

# Load the Hurricane Dataset
hurricane_df = pd.read_csv(r'https://raw.githubusercontent.com/JackOgozaly/Hurricane_Crop_Acreage/main/Data/historical_hurricane_date.csv')

hurricane_missing_val = hurricane_df.groupby("id")["64kt_wind_radii_SE_quad"].apply(lambda x: x.isnull().sum()).to_frame(name = "na_count")
hurricane_numb_obs = hurricane_df.groupby("id")["34kt_wind_radii_SE_quad"].apply(lambda x: len(x)).to_frame(name = "tot_count")
hurricane_no_pass = hurricane_df.groupby("id")["34kt_wind_radii_SE_quad"].apply(lambda x: sum(x==0)).to_frame(name = "zero_count")
filter_h = hurricane_missing_val.merge(hurricane_no_pass.merge(hurricane_numb_obs, how="left", on="id"), how="left", on="id").assign(
    prop_zero = lambda x: x.zero_count/x.tot_count
    )
# hurricane_ids = [i for i in hurricane_missing_val[hurricane_missing_val==0].index]
hurricane_ids = [i for i in filter_h[(filter_h.na_count==0)&(filter_h.prop_zero != 1)].index]

h_id = "AL012010"
def get_percentage(h_id):
    sub_hurricane_df = hurricane_df[hurricane_df["id"] == h_id].reset_index()
    sub_hurricane_df["ind"] = range(len(sub_hurricane_df))
    ## Interpolation
    df_new = pd.DataFrame()
    for i in range(len(sub_hurricane_df)-1):
        x = sub_hurricane_df.longitude[i:i+2]
        y = sub_hurricane_df.latitude[i:i+2]
        if x[i] != x[i+1]:
            y_hat = lowess(y, x, frac=1/6)
            xnew = np.linspace(x[i], x[i+1], 100)
            f_linear = interp1d(y_hat[:,0], y=y_hat[:,1], bounds_error=False, kind='linear', fill_value='extrapolate') 
            ynew_linear = f_linear(xnew)
        else:
            xnew = np.linspace(x[i], x[i+1], 100)
            ynew_linear = np.linspace(y[i], y[i+1], 100)
        df_new = pd.concat([df_new, pd.DataFrame({
            'ind': sub_hurricane_df.ind[i],
            'longitude': xnew,
            'latitude': ynew_linear})], axis=0)
    
    # Subset the ordinal directions
    ord = sub_hurricane_df[["ind","64kt_wind_radii_NE_quad", "64kt_wind_radii_SE_quad", "64kt_wind_radii_SW_qud", "64kt_wind_radii_NW_qud","50kt_wind_radii_NE_quad", "50kt_wind_radii_SE_quad", "50kt_wind_radii_SW_qud", "50kt_wind_radii_NW_qud","34kt_wind_radii_NE_quad", "34kt_wind_radii_SE_quad", "34kt_wind_radii_SW_qud", "34kt_wind_radii_NW_qud"]]
    # Add 1 to the "ind" so that we use the "next" ordinal directions for interpolation
    ord["ind"] -= 1
    df_new = df_new.merge(ord, on="ind", how="left")
    # Add the last observation to complete the dataset.
    df_new = pd.concat([df_new, sub_hurricane_df[-1:]], axis=0)[df_new.columns].reset_index() # reset index just in case
    # Compare the observed vs. Interpolated long/lat
    df_new["id"] = h_id
    hurricane_n = hurricane_shapefile(df_new, h_id, wind_speed=64)
    hurricane_tcn = hurricane_shapefile(df_new, h_id, wind_speed=50)
    hurricane_tsn = hurricane_shapefile(df_new, h_id, wind_speed=34)
    ## county shapefile and plot
    county_shapes = gpd.read_file('https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip')
    #Create our county shapefiles and plot them
    county_shapes_12 = county_shapes[county_shapes['STATEFP'] == '12'].reset_index(drop=True)
    county_shapes_12["full_area"] = county_shapes_12.geometry.area
    county_shapes_12 = county_shapes_12.drop(columns=["COUNTYNS", "AFFGEOID", "GEOID", "LSAD", "ALAND", "AWATER"])
    shape_n = unary_union(hurricane_n)
    p = gpd.GeoSeries(shape_n)
    shape_tcn = unary_union(hurricane_tcn)
    p1= gpd.GeoSeries(shape_tcn)
    shape_tsn = unary_union(hurricane_tsn)
    p2 = gpd.GeoSeries(shape_tsn)
    # changing series to dataframe
    def create_gdf(point_series):
        gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(point_series))
        gdf.set_crs("EPSG:4269", allow_override=True, inplace=True)
        return gdf

    pgdf = create_gdf(p)
    p1gdf = create_gdf(p1)
    p2gdf = create_gdf(p2)

    wind_speeds = ['64kt Winds', '50kt Winds', '34kt Winds']
    p_data = [pgdf, p1gdf, p2gdf]
    impacted_counties_final = pd.DataFrame()
    county_shapes_12.to_crs("EPSG:4269")
    for i, p in enumerate(p_data):
        impacted_counties = gpd.overlay(county_shapes_12, p, how='intersection', keep_geom_type=True)
        impacted_counties["Imp_area"] = impacted_counties.geometry.area
        impacted_counties_merged = pd.merge(county_shapes_12, impacted_counties[["COUNTYFP", "Imp_area"]], on='COUNTYFP', how='left')
        impacted_counties_merged["Imp_area"] = impacted_counties_merged["Imp_area"].fillna(0)
        impacted_counties_merged['Wind_Speed'] = wind_speeds[i]
        impacted_counties_final = impacted_counties_final.append(impacted_counties_merged)
    # Wind to Columns
    df_wind_col = (
        impacted_counties_final
        .pivot_table(values='Imp_area', index='COUNTYFP', columns='Wind_Speed', aggfunc='sum', fill_value='')
    )
    # Merged datasets
    df_merged_final = pd.merge(county_shapes_12, df_wind_col, how='left', on='COUNTYFP')
    df_merged_final['34kt_pct'] = (df_merged_final["34kt Winds"] / df_merged_final["full_area"])*100
    df_merged_final['50kt_pct'] = (df_merged_final["50kt Winds"] / df_merged_final["full_area"])*100
    df_merged_final['64kt_pct'] = (df_merged_final["64kt Winds"] / df_merged_final["full_area"])*100
    df_merged_final = df_merged_final.reindex(columns=['STATEFP', 'COUNTYFP', 'NAME', 'geometry','full_area', '34kt Winds', '34kt_pct', '50kt Winds', '50kt_pct', '64kt Winds', '64kt_pct'])
    df_merged = df_merged_final.drop(columns=["geometry","STATEFP"])
    df_merged ["id"] = h_id 
    df_merged["year"] = pd.to_datetime(sub_hurricane_df.date[0]).year
    return(df_merged)

dff = pd.DataFrame()
for i in hurricane_ids:
    dff = dff.append(get_percentage(i))
    

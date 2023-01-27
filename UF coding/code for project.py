#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 17:48:19 2022

@author: xiaomeihai
"""

import geopandas as gpd
from geopandas import GeoSeries, GeoDataFrame
import matplotlib.pyplot as plt 

shapes = gpd.read_file(r'AL112017_radii.shp')
shapes = shapes[shapes['RADII'] == 34]
shapes = shapes[shapes['SYNOPTIME'] == '2017083118']
fig, ax = plt.subplots(figsize=(10,8))
hurricane = gpd.GeoSeries(shapes['geometry'])
hurricane.plot(ax=ax, color='blue', edgecolor='black')
plt.title("Example Hurricane Irma Shapefile Made by NOAA")

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

hurricane_df = pd.read_csv(r'https://raw.githubusercontent.com/JackOgozaly/Hurricane_Crop_Acreage/main/Data/historical_hurricane_date.csv')

#Use our function for the different wind speeds
hurricane_irma_ts = hurricane_shapefile(hurricane_df, 'AL112017', wind_speed=34)
hurricane_irma_tc = hurricane_shapefile(hurricane_df, 'AL112017', wind_speed=50)
hurricane_irma_h = hurricane_shapefile(hurricane_df, 'AL112017', wind_speed=64)

#version of the wind speed plot
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches

county_shapes = gpd.read_file('https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip')


#plt.style.use('default')
fig, ax = plt.subplots(figsize=(10,8))

#Create our county shapefiles and plot them
counties = gpd.GeoSeries(county_shapes['geometry'])
counties.plot(ax=ax, color='white', edgecolor='black')

irma_shape = unary_union(hurricane_irma_ts)
p = gpd.GeoSeries(irma_shape)
p.plot(ax=ax, color = "yellow", alpha = .4)

irma_shape2 = unary_union(hurricane_irma_tc)
p2 = gpd.GeoSeries(irma_shape2)
p2.plot(ax=ax, color = "orange", alpha = .4)


irma_shape3 = unary_union(hurricane_irma_h)
p3 = gpd.GeoSeries(irma_shape3)
p3.plot(ax=ax, color = "red", alpha = .4)

#This is also just for decoration
plt.xlim([-90, -50])
plt.ylim([15, 40])
plt.axis('off')
plt.title("Hurricane Irma Wind Speeds Map")
fig.patch.set_facecolor('lightblue')
ax = plt.gca()
ax.patch.set_facecolor('white')

#Creating our legend manually
legend1 = mpatches.Patch(color='yellow', label='34kt Winds', alpha =.4)
legend2 = mpatches.Patch(color='orange', label='50kt Winds', alpha =.4)
legend3 = mpatches.Patch(color='red', label='64kt Winds', alpha =.4)
plt.legend(handles=[legend1, legend2, legend3])


##changing series to dataframe 
pgdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(p))
p2gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(p2))
p3gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(p3))


##intersecting and adding area to data frame
impacted_counties = gpd.overlay(county_shapes, pgdf, how = 'intersection')
impacted_counties.plot()
impacted_counties['Area']=impacted_counties.area
impacted_counties['Wind_Speed']= '34kt Winds'

impacted_counties2 = gpd.overlay(county_shapes, p2gdf, how = 'intersection')
impacted_counties2.plot()
impacted_counties2['Area']=impacted_counties2.area
impacted_counties2['Wind_Speed']= '50kt Winds' 


impacted_counties3 = gpd.overlay(county_shapes, p3gdf, how = 'intersection')
impacted_counties3.plot()
impacted_counties3['Area']=impacted_counties3.area
impacted_counties3['Wind_Speed']= '64kt Winds'

#Impacted counties dataframe appended together
impacted_counties_final = impacted_counties.append([impacted_counties2, impacted_counties3]) 
impacted_counties['ALAND'].value_counts()
#impacted_counties_final = impacted_counties_final.rename(columns={'COUNTYFP': 'county_code'})

#Cause of loss
cause_of_loss_4 = pd.read_parquet('https://github.com/JackOgozaly/Hurricane_Crop_Acreage/blob/main/Data/Cause_Of_Loss/crop_loss_data_4.parquet.gzip?raw=true')
cause_of_loss_4 = cause_of_loss_4.rename(columns = {'county_code':'COUNTYFP', 'state_code':'STATEFP'})


#idea for slicing df
df12 = county_shapes[['STATEFP','COUNTYFP','GEOID','NAME','geometry']]
df12 = df12[df12['STATEFP'] == '12']

df_wind_col = pd.merge(df12, impacted_counties_final, how = 'left', on = 'COUNTYFP')
df_wind_col = df_wind_col.groupby(['COUNTYFP', 'Wind_Speed'])['Area'].sum().reset_index()
df_wind_col = df_wind_col.set_index(['COUNTYFP', 'Wind_Speed']).Area.unstack(fill_value='')
df_wind_col

#Merged datasets
county_shapes_merge = county_shapes.copy()
county_shapes_merge = county_shapes_merge[county_shapes_merge['STATEFP'] == '12']
df_merged_final = pd.merge(county_shapes_merge, df_wind_col, how = 'left', on = 'COUNTYFP')

#Cause of loss clean
cause_of_loss_4_clean = cause_of_loss_4.drop(columns=["insurance_plan_name_abbreviation","subsidy","state/private_subsidy","additional_subsidy","efa_premium_discount","producer_paid_premium","insurance_plan_code","stage_code","net_endorsed_acres"])
cause_of_loss_4_clean = cause_of_loss_4_clean[cause_of_loss_4_clean['year_of_loss'] == 2017]
cause_of_loss_4_clean = cause_of_loss_4_clean[cause_of_loss_4_clean['cause_of_loss_description'] == 'Hurricane/Tropical Depression']
cause_of_loss_4_clean = cause_of_loss_4_clean[cause_of_loss_4_clean['STATEFP'] == '12'].reset_index(drop=True)
cause_of_loss_4_clean = cause_of_loss_4_clean.drop(columns = ['commodity_year_identifier','commodity_code','coverage_category','cause_of_loss_code','policies_earning_premium','policies_identified'])


#final df
final_df_cause_of_loss = cause_of_loss_4_clean.merge(df_merged_final, on='COUNTYFP', how='left')
































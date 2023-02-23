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
county_shapes_12 = county_shapes[county_shapes['STATEFP'] == '12'].reset_index(drop=True)
county_shapes_12 = county_shapes_12.set_crs("EPSG:3347", allow_override=True)
county_shapes_12["full_area"] = county_shapes_12.geometry.area
county_shapes_12 = county_shapes_12.drop(columns=["COUNTYNS","AFFGEOID","GEOID","LSAD","ALAND","AWATER"])

#plt.style.use('default')
fig, ax = plt.subplots(figsize=(10,8))

#Create our county shapefiles and plot them
counties = gpd.GeoSeries(county_shapes_12['geometry'])
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
pgdf  = pgdf.set_crs("EPSG:3347", allow_override=True)
p2gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(p2))
p2gdf = p2gdf.set_crs("EPSG:3347", allow_override=True)
p3gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(p3))
p3gdf = p3gdf.set_crs("EPSG:3347", allow_override=True)

##intersecting and adding area to data frame
impacted_counties = gpd.overlay(pgdf, county_shapes_12,  how = 'intersection', keep_geom_type=True)
impacted_counties.plot()
impacted_counties["Imp_area"] = impacted_counties.geometry.area
impacted_counties['OverlapIntersection']=(impacted_counties["Imp_area"] / impacted_counties["full_area"])*100
impacted_counties['Wind_Speed']= '34kt Winds'
impacted_counties = impacted_counties.reindex(columns=['STATEFP', 'COUNTYFP', 'NAME', 'geometry', 'full_area', 'Imp_area', 'OverlapIntersection','Wind_Speed'])

## need to add 0 to N/A county before calculate 
impacted_counties2 = gpd.overlay(county_shapes_12, p2gdf, how = 'intersection', keep_geom_type=True)
impacted_counties2["Imp_area"] = impacted_counties2.geometry.area
impacted_counties2.plot()

impacted_counties2_merged = pd.merge(county_shapes_12,impacted_counties2[["COUNTYFP","Imp_area"]],on='COUNTYFP', how='left')
impacted_counties2_merged["Imp_area"] = impacted_counties2_merged["Imp_area"].fillna(0)

impacted_counties2_merged['OverlapIntersection']= (impacted_counties2_merged["Imp_area"] / impacted_counties2_merged["full_area"])*100
impacted_counties2_merged['Wind_Speed']= '50kt Winds'

## need to add 0 to N/A county before calculate 
impacted_counties3 = gpd.overlay(county_shapes_12, p3gdf, how = 'intersection', keep_geom_type=True)
impacted_counties3["Imp_area"] = impacted_counties3.geometry.area
impacted_counties3.plot()

impacted_counties3_merged = pd.merge(county_shapes_12,impacted_counties3[["COUNTYFP","Imp_area"]],on='COUNTYFP', how='left')
impacted_counties3_merged["Imp_area"] = impacted_counties3_merged["Imp_area"].fillna(0)

impacted_counties3_merged['OverlapIntersection']= (impacted_counties3_merged["Imp_area"] / impacted_counties3_merged["full_area"])*100
impacted_counties3_merged['Wind_Speed']= '64kt Winds'

##Appending three dfs
impacted_counties_final = impacted_counties.append([impacted_counties2_merged, impacted_counties3_merged]) 

##Wind to Columns
df_wind_col = impacted_counties_final.groupby(['COUNTYFP', 'Wind_Speed'])[['Imp_area']].sum().reset_index()
df_wind_col = df_wind_col.set_index(['COUNTYFP', 'Wind_Speed']).Imp_area.unstack(fill_value='')

##Merged datasets
df_merged_final = pd.merge(county_shapes_12, df_wind_col, how = 'left', on = 'COUNTYFP')
df_merged_final['34kt_pct'] = (df_merged_final["34kt Winds"] / df_merged_final["full_area"])*100
df_merged_final['50kt_pct'] = (df_merged_final["50kt Winds"] / df_merged_final["full_area"])*100
df_merged_final['64kt_pct'] = (df_merged_final["64kt Winds"] / df_merged_final["full_area"])*100
df_merged_final = df_merged_final.reindex(columns=['STATEFP', 'COUNTYFP', 'NAME', 'geometry', 'full_area', '34kt Winds', '34kt_pct','50kt Winds', '50kt_pct', '64kt Winds', '64kt_pct'])
df_merged_final = df_merged_final.drop(columns=["geometry"])


#Cause of loss
cause_of_loss = pd.read_parquet('https://github.com/JackOgozaly/Hurricane_Crop_Acreage/blob/main/Data/Cause_Of_Loss/crop_loss_data_4.parquet.gzip?raw=true')
cause_of_loss = cause_of_loss.rename(columns = {'county_code':'COUNTYFP', 'state_code':'STATEFP'})

#Cause of loss clean
cause_of_loss_clean = cause_of_loss[cause_of_loss['year_of_loss'] == 2017]
cause_of_loss_clean = cause_of_loss_clean[(cause_of_loss_clean['month_of_loss'] == '8')|(cause_of_loss_clean['month_of_loss'] == '9')]
cause_of_loss_clean = cause_of_loss_clean[cause_of_loss_clean['cause_of_loss_description'] == 'Hurricane/Tropical Depression']
cause_of_loss_clean = cause_of_loss_clean[cause_of_loss_clean['STATEFP'] == '12'].reset_index(drop=True)
cause_of_loss_clean1 = cause_of_loss_clean.drop(cause_of_loss_clean.columns[[0,3,5,6,7,8,9,13,14,16,19,20,21,22,23,26,27]],axis = 1)

cause_of_loss_clean_group=cause_of_loss_clean1.groupby(['COUNTYFP'])[['total_premium','indemnity_amount','net_planted_quantity','liability','net_determined_quantity']].sum()

##Merge df_merged_final with cause of loss
final_data = pd.merge(df_merged_final,cause_of_loss_clean_group,on='COUNTYFP', how='left').fillna(0)

#Precipitation csv
precipitation = pd.read_excel('precipitation.xlsx')
PlantedAcres = pd.read_excel('PlantedAcres.xlsx')

#Merge cause of loss clean with preciptation dataset
final_data_pre = pd.merge(final_data,precipitation,on='NAME', how='left')
final_data_pre = pd.merge(final_data_pre,PlantedAcres,on='NAME', how='left')


##prepare the dataset for the model
## drop the columns not needed
data1 = final_data_pre.drop (columns=["NAME", "STATEFP", "COUNTYFP","total_premium","net_planted_quantity","net_determined_quantity","liability"])
##rearrange the columns as indemnity_amount first column
data1 = data1.iloc[:, [7] + list(range(7)) + list(range(8, len(data1.columns)))]


# library & dataset
import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm
import statsmodels.formula.api as sm

# check the correlation bewteen Y variable and the other variables in regression
sns.pairplot(data1 , kind="reg", diag_kws={"bins":"sqrt"})
plt.show()

##only limited to first 5 columns 
sns.pairplot(data1.iloc[:, :5], kind="reg", diag_kws={"bins":"sqrt"})
plt.show() 

##heatmaps
sns.heatmap(data1.corr(), annot=True)


data2 = data1.astype("double")
data2 = data2.rename(columns={"64kt_pct":"sixfour_pct"})
data2["sqrt_indemnity_amount"] = np.sqrt(data2["indemnity_amount"])
 
##model1
lm = sm.ols(formula = "sqrt_indemnity_amount ~ full_area + sixfour_pct + PlantedAcres + pcp", data=data2).fit()
print(lm.summary())


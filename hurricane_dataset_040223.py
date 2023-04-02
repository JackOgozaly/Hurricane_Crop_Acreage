
#Used for dealing with our dataframe
import pandas as pd
#Using for finding lat and long coords X nautical miles out
import geopy 
import geopy.distance
from geopy.distance import geodesic
#Used for drawing circles
from numpy import cos,sin
import numpy as np
#Used for polygon shenanigans
import shapely
from shapely.ops import unary_union
#Used for visualizing
import geopandas as gpd
#Used to join multiple dataframes together
from functools import reduce


#I lifted these two circle functions from this stackoverflow post
#https://stackoverflow.com/questions/11331854/how-can-i-generate-an-arc-in-numpy
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


def intersection_check(county_shapefiles, geoid_list, hurricane_shapefile_list,
                       wind_speed, hurricane_catgory = None):
    '''
    Parameters
    ----------
    county_shapefiles : Series
        column of county polygons
    geoid_list : list
        list of GOEIDs in same order as county shapefiles
    hurricane_shapefile_list : list
        Output from hurricane_shapefile f
    wind_speed : INT
        What wind speed your shapefile is. The default is 34.

    Returns
    -------
    Dataframe with counties that intersected + how many times they intersected, and wind speed

    '''
    counties = gpd.GeoSeries(county_shapefiles)
    county_fips = geoid_list
    #Create an empty list to store our intersect values
    intersects_list = []
    area_list = []
    for i in range(len(hurricane_shapefile_list)):
        hurricane_sub_shape = hurricane_shapefile_list[i]
        for i in range(len(county_shapes)):
            p1 = counties.iloc[i]
            #Store a list of counties that intersected with irma
            if p1.intersects(hurricane_sub_shape):
                intersects_list.append(county_fips[i])
                
                #Grab our area intersect value
                area_intersect = (p1.intersection(hurricane_sub_shape).area/p1.area)
                if hurricane_catgory is not None:    
                    area_list.append({'GEOID':county_fips[i], f'area_{hurricane_catgory}':area_intersect})
                else:
                    area_list.append({'GEOID':county_fips[i], f'area_{wind_speed}kt':area_intersect})
                
    if not intersects_list:
        return(None)
       
    #Create a dataframe with the count of times a storm was over a coutny
    count_df = (pd.DataFrame(intersects_list).groupby(0).size()).reset_index(drop=False, name='count').sort_values('count', ascending = False)
    
    #Create a dataframe with the max coverage a storm had over a county
    area_df = pd.DataFrame(area_list)
    
    if hurricane_catgory is not None:
        count_df.columns = ['GEOID', f'count_{hurricane_catgory}_intersects']
        area_df = area_df.groupby(['GEOID'])[f'area_{hurricane_catgory}'].max().reset_index(drop=False)
    else:
        count_df.columns = ['GEOID', f'count_{wind_speed}kt_intersects']
        area_df = area_df.groupby(['GEOID'])[f'area_{wind_speed}kt'].max().reset_index(drop=False)
    
    #Bring our count and max area dataframes together
    final_df = pd.merge(area_df, count_df, how = 'left', on = 'GEOID')
    
    #Return output
    return(final_df)


####____________________ Our code__________________________####

#Read in our county shapefile
county_shapes = gpd.read_file('https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip')


#Read in our historical hurricane dataframe
hurricane_df = pd.read_csv(r'https://raw.githubusercontent.com/JackOgozaly/Hurricane_Crop_Acreage/main/Data/historical_hurricane_date.csv')

#Remove rows where we have no shape data
hurricane_df = hurricane_df.dropna(subset=['34kt_wind_radii_NE_quad', '34kt_wind_radii_SE_quad',
                                           '34kt_wind_radii_SW_qud', '34kt_wind_radii_NW_qud' ],
                                   how='all')


#I'm sure there's a better way to do this, but I am lazy
hurricane_df['Category'] = np.where(hurricane_df['max_sustained_wind'].between(74, 95, inclusive ='both'), 'Cat1', np.NaN)
hurricane_df['Category'] = np.where(hurricane_df['max_sustained_wind'].between(76, 110, inclusive ='both'), 'Cat2', hurricane_df['Category'])
hurricane_df['Category'] = np.where(hurricane_df['max_sustained_wind'].between(111, 129, inclusive ='both'), 'Cat3', hurricane_df['Category'])
hurricane_df['Category'] = np.where(hurricane_df['max_sustained_wind'].between(130, 156, inclusive ='both'), 'Cat4', hurricane_df['Category'])
hurricane_df['Category'] = np.where(hurricane_df['max_sustained_wind'].between(157, 500, inclusive ='both'), 'Cat5', hurricane_df['Category'])

#Get a print out of how many of each category we have
print(hurricane_df.groupby(["Category"]).size())


#The cause of loss dataset only goes back to 1988, so need to summarize hurricanes before then
hurricane_df['date'] = pd.to_datetime(hurricane_df['date'])
hurricane_df = hurricane_df[hurricane_df['date'].dt.year > 1989]

#Get a list of our hurricane ids so we can loop through them
hurricane_ids = list(hurricane_df['id'].unique())
#Create an empty datframe to store our hurricane summaries
hurricane_summary_dfs = []


for hurricane_id in hurricane_ids:
    #Create our shapefiles for the different wind speeds
    hurricane_ts_shape = hurricane_shapefile(hurricane_df, hurricane_id, wind_speed=34)
    hurricane_tc_shape = hurricane_shapefile(hurricane_df, hurricane_id, wind_speed=50)
    hurricane_h_shape = hurricane_shapefile(hurricane_df, hurricane_id, wind_speed=64)


    #Again, I'm sure there's a better way, but I am lazy
    category_hurricanes = hurricane_df[hurricane_df['id'] == hurricane_id].copy()
    cat1_df = category_hurricanes[category_hurricanes['Category'] == 'Cat1'].copy()
    cat2_df = category_hurricanes[category_hurricanes['Category'] == 'Cat2'].copy()
    cat3_df = category_hurricanes[category_hurricanes['Category'] == 'Cat3'].copy()
    cat4_df = category_hurricanes[category_hurricanes['Category'] == 'Cat4'].copy()
    cat5_df = category_hurricanes[category_hurricanes['Category'] == 'Cat5'].copy()
    #Create our shapefiles
    cat1_shape = hurricane_shapefile(cat1_df, hurricane_id, wind_speed=64)
    cat2_shape = hurricane_shapefile(cat2_df, hurricane_id, wind_speed=64)
    cat3_shape = hurricane_shapefile(cat3_df, hurricane_id, wind_speed=64)
    cat4_shape = hurricane_shapefile(cat4_df, hurricane_id, wind_speed=64)
    cat5_shape = hurricane_shapefile(cat5_df, hurricane_id, wind_speed=64)


    #Perform our intersection analysis
    hurricane_ts_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 hurricane_ts_shape, 
                                                 wind_speed= 34)
    
    
    hurricane_tc_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 hurricane_tc_shape, 
                                                 wind_speed= 50)
    
    hurricane_h_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 hurricane_h_shape,
                                                 wind_speed = 64)


    cat1_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 cat1_shape,
                                                 wind_speed = 64,
                                                 hurricane_catgory = 'cat1')


    cat2_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 cat2_shape,
                                                 wind_speed = 64,
                                                 hurricane_catgory = 'cat2')
    
    cat3_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 cat3_shape,
                                                 wind_speed = 64,
                                                 hurricane_catgory = 'cat3')

    cat4_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 cat4_shape,
                                                 wind_speed = 64,
                                                 hurricane_catgory = 'cat4')


    cat5_intersects = intersection_check(county_shapes['geometry'],
                                                 list(county_shapes['GEOID']), 
                                                 cat5_shape,
                                                 wind_speed = 64,
                                                 hurricane_catgory = 'cat5')

    df_list = [hurricane_ts_intersects, hurricane_tc_intersects, hurricane_h_intersects,
               cat1_intersects, cat2_intersects, cat3_intersects, cat4_intersects, cat5_intersects]
    #Remove any null lists
    df_list = [df for df in df_list if df is not None]
    df_list = [df for df in df_list if not df.empty]

    if df_list:
        intersects_df = reduce(lambda x, y: pd.merge(x, y, on = 'GEOID', how = 'outer'), df_list)
        intersects_df['hurricane_id'] = hurricane_id
        hurricane_summary_dfs.append(intersects_df)

    else: 
        continue

hurricane_summary_dfs1 = [df for df in hurricane_summary_dfs if not df.empty]
#Combine all the dataframes together
hurricane_summary = pd.concat(hurricane_summary_dfs1)

#Reorganize columns for sanity check
final_df = hurricane_summary[['GEOID', 'hurricane_id', 'area_34kt', 'count_34kt_intersects', 
                              'area_50kt', 'count_50kt_intersects', 'area_64kt', 'count_64kt_intersects',
                              'area_cat1', 'count_cat1_intersects', 'area_cat2', 'count_cat2_intersects',
                              'area_cat3', 'count_cat3_intersects', 'area_cat4', 'count_cat4_intersects']].copy()

final_df['hurricane_category'] = np.where(final_df['count_cat1_intersects'].notnull(), 'cat1', 'TS')
final_df['hurricane_category'] = np.where(final_df['count_cat2_intersects'].notnull(), 'cat2', final_df['hurricane_category'])
final_df['hurricane_category'] = np.where(final_df['count_cat3_intersects'].notnull(), 'cat3', final_df['hurricane_category'])
final_df['hurricane_category'] = np.where(final_df['count_cat4_intersects'].notnull(), 'cat4', final_df['hurricane_category'])

final_df.to_csv(r'/Users/jackogozaly/Desktop/Python_Directory/Hurricane/Data/hurricane_final.csv', index=False)


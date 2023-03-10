#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:22:44 2023

@author: xiaomeihai
"""

from numpy import cos, sin, arccos
import geopy
from sklearn import metrics
from sklearn import linear_model as lm
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_score, train_test_split
from sklearn.linear_model import LinearRegression
from scipy.stats import boxcox
import statsmodels.formula.api as sm
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from shapely.ops import unary_union
from shapely.geometry import shape
from shapely.geometry import Polygon
import shapely
import numpy as np
from geopy.distance import geodesic
import geopy.distance
import pandas as pd
import geopandas as gpd
from geopandas import GeoSeries, GeoDataFrame
import matplotlib.pyplot as plt


## function 1 county shapefile 

county_shapes = gpd.read_file('https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip')
county_shapes_12 = county_shapes[county_shapes['STATEFP'] == '12'].reset_index(drop=True)
county_shapes_12 = county_shapes_12.set_crs("EPSG:3347", allow_override=True)
county_shapes_12["full_area"] = county_shapes_12.geometry.area
county_shapes_12 = county_shapes_12.drop(columns=["COUNTYNS", "AFFGEOID", "GEOID", "LSAD", "ALAND", "AWATER"])


## function 2- creating hurricane shapefile 

# Used for dealing with our dataframe
# Using for finding lat and long coords X nautical miles out
# Used for drawing circles
# Used for polygon shenanigans
# Used for visualizing
############################################################################################################################################
# Run Below

def parametric_circle(t, xc, yc, R):
    x = xc + R*cos(t)
    y = yc + R*sin(t)
    return x, y


def inv_parametric_circle(x, xc, R):
    t = np.arccos(((x-xc)/R))
    return t


def hurricane_shapefile(df, storm_id, wind_speed=34, quadrant_points=30):
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
    # Empty list to store our shapes
    storm_shapes = []
    # Bearing for our wind directions
    bearings = [0, 90, 180, 270]

    # Subset for our storm
    sub_df = df[df['id'] == storm_id]

    for i in range(len(sub_df)):
        # Get the center of our storm
        latitude = sub_df['latitude'].iloc[i]
        longitude = sub_df['longitude'].iloc[i]
        center = geopy.Point(latitude, longitude)

        wind_columns = [f'{wind_speed}kt_wind_radii_NE_quad',
                        f'{wind_speed}kt_wind_radii_SE_quad',
                        f'{wind_speed}kt_wind_radii_SW_qud',
                        f'{wind_speed}kt_wind_radii_NW_qud']

        # Empty list to store the quadrant shapefiles
        quadrant_shapes = []
        for j in range(4):
            wind_miles = sub_df[wind_columns[j]].iloc[i]

            if wind_miles == 0:
                continue

           # Get the two points that make up the ends of our arc
            p1 = geodesic(nautical=wind_miles).destination(
                point=geopy.Point(latitude, longitude), bearing=bearings[j])
            p2 = geodesic(nautical=wind_miles).destination(
                point=geopy.Point(latitude, longitude), bearing=(bearings[j] + 90))

            # Treat our radius different depending on the quadrant
            if j == 0:
                radius = (p1.latitude - center.latitude)
            if j == 1:
                radius = (p2.latitude - center.latitude) * -1
            elif j == 2:
                radius = (p1.latitude - center.latitude)
            elif j == 3:
                radius = (p2.latitude - center.latitude) * -1

            # Convert our start and end points into tuples
            start_point = tuple([p1.latitude, p1.longitude])
            end_point = tuple([p2.latitude, p2.longitude])
            # Draw our semi-circle
            start_t = inv_parametric_circle(start_point[0], latitude, radius)
            end_t = inv_parametric_circle(end_point[0], latitude, radius)
            arc_T = np.linspace(start_t, end_t, quadrant_points)
            # Get the points from our circle
            X, Y = parametric_circle(arc_T, latitude, longitude, radius)
            points = list(tuple(zip(Y, X)))
            # Insert the center of the storm as a first point so we can make the polygon
            points.insert(0, (longitude, latitude))
            # Create our polygon
            sub_poly = shapely.geometry.shape(
                {'type': "Polygon", "coordinates": [points]})
            # Append to our list
            quadrant_shapes.append(sub_poly)

        # Union together each quadrant to form a circle-ish
        storm_shapes.append(unary_union(quadrant_shapes))
    # Remove any records that had 0 for all wind radii
    storm_shapes = [i for i in storm_shapes if not i.is_empty]

    return(storm_shapes)


############################################################################################################################################

## function 3- overlap area for county shapefile and hurricane shapefile 

hurricane_df = pd.read_csv(r'https://raw.githubusercontent.com/JackOgozaly/Hurricane_Crop_Acreage/main/Data/historical_hurricane_date.csv')
h_list= hurricane_df.id.unique()


##'AL112017'= Irma, 'al052019'=DORIAN, 'al142016'=MATTHEW, 'al142018'=MICHAEL, 'al192020'=SALLY

h_list = ['AL112017', 'AL052019', 'AL142016', 'AL142018', 'AL192020']

h_full_df = gpd.GeoDataFrame()
for X in h_list:
    
    hurricane_ts = hurricane_shapefile(hurricane_df, X, wind_speed=34)
    hurricane_tc = hurricane_shapefile(hurricane_df, X, wind_speed=50)
    hurricane_h = hurricane_shapefile(hurricane_df, X, wind_speed=64)
    
    shape_ts = unary_union(hurricane_ts)
    p = gpd.GeoDataFrame(geometry=gpd.GeoSeries(shape_ts))
    shape_tc = unary_union(hurricane_tc)
    p2= gpd.GeoDataFrame(geometry=gpd.GeoSeries(shape_tc))
    shape_h = unary_union(hurricane_h)
    p3= gpd.GeoDataFrame(geometry=gpd.GeoSeries(shape_h))
    
    # changing series to dataframe
    pgdf = p.set_crs("EPSG:3347", allow_override=True)
    p2gdf = p2.set_crs("EPSG:3347", allow_override=True)
    p3gdf = p3.set_crs("EPSG:3347", allow_override=True)
    
    # intersecting and adding area to data frame
    impacted_counties1 = gpd.overlay(pgdf, county_shapes_12,  how='intersection', keep_geom_type=True)
    impacted_counties1["Imp_area"] = impacted_counties1.geometry.area
    impacted_counties1_merged = pd.merge(county_shapes_12, impacted_counties1[["COUNTYFP", "Imp_area"]], on='COUNTYFP', how='left')
    impacted_counties1_merged["Imp_area"] = impacted_counties1_merged["Imp_area"].fillna(0)
    impacted_counties1_merged['OverlapIntersection'] = (impacted_counties1_merged["Imp_area"] / impacted_counties1_merged["full_area"])*100
    impacted_counties1_merged['Wind_Speed'] = '34kt Winds'
    
    
    # need to add 0 to N/A county before calculate
    impacted_counties2 = gpd.overlay(county_shapes_12, p2gdf, how='intersection', keep_geom_type=True)
    impacted_counties2["Imp_area"] = impacted_counties2.geometry.area
    impacted_counties2_merged = pd.merge(county_shapes_12, impacted_counties2[["COUNTYFP", "Imp_area"]], on='COUNTYFP', how='left')
    impacted_counties2_merged["Imp_area"] = impacted_counties2_merged["Imp_area"].fillna(0)
    impacted_counties2_merged['OverlapIntersection'] = (impacted_counties2_merged["Imp_area"] / impacted_counties2_merged["full_area"])*100
    impacted_counties2_merged['Wind_Speed'] = '50kt Winds'
    
    # need to add 0 to N/A county before calculate
    impacted_counties3 = gpd.overlay(county_shapes_12, p3gdf, how='intersection', keep_geom_type=True)
    impacted_counties3["Imp_area"] = impacted_counties3.geometry.area
    impacted_counties3_merged = pd.merge(county_shapes_12, impacted_counties3[["COUNTYFP", "Imp_area"]], on='COUNTYFP', how='left')
    impacted_counties3_merged["Imp_area"] = impacted_counties3_merged["Imp_area"].fillna(0)
    impacted_counties3_merged['OverlapIntersection'] = (impacted_counties3_merged["Imp_area"] / impacted_counties3_merged["full_area"])*100
    impacted_counties3_merged['Wind_Speed'] = '64kt Winds'
    
    # Appending three dfs
    impacted_counties_final = impacted_counties1_merged.append([impacted_counties2_merged, impacted_counties3_merged])
    impacted_counties_final ["storm_id"] = X
    
    # Append to the h_full_shape
    h_full_df = pd.concat([h_full_df, impacted_counties_final], axis=0)
    
    
# After running the for loop    
df_wind_col = h_full_df.groupby(['COUNTYFP', 'storm_id','Wind_Speed'])[['Imp_area']].sum().reset_index()
df_wind_col2 = df_wind_col.set_index(['COUNTYFP', 'storm_id', 'Wind_Speed']).Imp_area.unstack(fill_value='')
    
df_merged_final1 = pd.merge(county_shapes_12, df_wind_col2, how='outer', on='COUNTYFP')

df_merged_final = pd.merge(county_shapes_12, df_wind_col2, how='left', on='COUNTYFP')
df_merged_final['34kt_pct'] = (df_merged_final["34kt Winds"] / df_merged_final["full_area"])*100
df_merged_final['50kt_pct'] = (df_merged_final["50kt Winds"] / df_merged_final["full_area"])*100
df_merged_final['64kt_pct'] = (df_merged_final["64kt Winds"] / df_merged_final["full_area"])*100
df_merged_final = df_merged_final.reindex(columns=['STATEFP', 'COUNTYFP', 'NAME', 'geometry','full_area', '34kt Winds', '34kt_pct', '50kt Winds', '50kt_pct', '64kt Winds', '64kt_pct'])
df_merged_final = df_merged_final.drop(columns=["geometry"])

############################################################################################################################################


# function 4- Cause of loss
cause_of_loss = pd.read_parquet('https://github.com/JackOgozaly/Hurricane_Crop_Acreage/blob/main/Data/Cause_Of_Loss/crop_loss_data_4.parquet.gzip?raw=true')
cause_of_loss = cause_of_loss.rename(columns={'county_code': 'COUNTYFP', 'state_code': 'STATEFP'})

Y_list= cause_of_loss.year_of_loss.unique()
M_list = cause_of_loss.month_of_loss.unique()
states = ['12', '34', '56']  # list of state codes to iterate over

for state in states:
    for year in Y_list:
        for month in M_list:
            cause_of_loss_clean = cause_of_loss[(cause_of_loss['year_of_loss'] == year) & (cause_of_loss['month_of_loss'] == month)]
            cause_of_loss_clean = cause_of_loss_clean[cause_of_loss_clean['cause_of_loss_description'] == 'Hurricane/Tropical Depression']
            cause_of_loss_clean = cause_of_loss_clean[cause_of_loss_clean['STATEFP'] == state].reset_index(drop=True)
            cause_of_loss_clean1 = cause_of_loss_clean.drop(cause_of_loss_clean.columns[[0, 3, 5, 6, 7, 8, 9, 13, 14, 16, 19, 20, 21, 22, 23, 26, 27]], axis=1)

            cause_of_loss_clean_group = cause_of_loss_clean1.groupby(['COUNTYFP'])[['total_premium', 'indemnity_amount', 'net_planted_quantity', 'liability', 'net_determined_quantity']].sum()

            # Do something with the cause_of_loss_clean_group DataFrame for this year, month, and state


# Merge df_merged_final with cause of loss
final_data = pd.merge(df_merged_final, cause_of_loss_clean_group,on='COUNTYFP', how='left').fillna(0)

# Precipitation csv https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/county/mapping/8/pcp/201709/2/value
precipitation = pd.read_excel('precipitation.xlsx')

# Crop Acreage Data: https://www.fsa.usda.gov/news-room/efoia/electronic-reading-room/frequently-requested-information/crop-acreage-data/index
PlantedAcres = pd.read_excel('PlantedAcres.xlsx')

# Merge cause of loss clean with preciptation dataset
final_data_pre = pd.merge(final_data, precipitation, on='NAME', how='left')
final_data_pre = pd.merge(final_data_pre, PlantedAcres, on='NAME', how='left')
final_data_pre.to_csv('/Users/xiaomeihai/Desktop//Project/final_data_pre.csv')


# prepare the dataset for the model
# drop the columns not needed
data1 = final_data_pre.drop(columns=["NAME", "STATEFP", "COUNTYFP", "total_premium","net_planted_quantity", "net_determined_quantity", "liability"])

# rearrange the columns as indemnity_amount first column
data1 = data1.iloc[:, [7] + list(range(7)) +list(range(8, len(data1.columns)))]


data1 = data1.astype("double")
data1 = data1.rename(columns={"64kt_pct": "sixfour_pct","34kt_pct": "threefour_pct", "50kt_pct": "fifity_pct"})

data1 = data1.rename(columns={"34kt Winds": "34kt_Winds", "50kt Winds": "50kt_Winds","64kt Winds": "64kt_Winds"})


# Transformations For Better Normal Distribution
# Log Transformation
data1["log_indemnity_amount"] = np.log1p(data1["indemnity_amount"])


# Square-Root Transformation
data1["sqrt_indemnity_amount"] = np.sqrt(data1["indemnity_amount"])

# Reciprocal Transformation- cannot be used as lot of inf values inside
data1["reciprocal_indemnity_amount"] = 1/data1["indemnity_amount"]


data1.to_csv('/Users/xiaomeihai/Desktop//Project/data1.csv')

# Box-Cox Transformation- cannot be used since we have 0 in indemnity_amount
#data1["bcx_indemnity_amount"] = bcx_indemnity_amount, lam = boxcox(data1["indemnity_amount"])
# lam is the best lambda for the distribution


datalog = data1.drop(columns=["indemnity_amount", "34kt_Winds", "50kt_Winds","reciprocal_indemnity_amount","threefour_pct", "fifity_pct", "64kt_Winds","sqrt_indemnity_amount"])
datalog = datalog.iloc[:, [4, 0, 1,2,3]]


datasqrt = data1.drop(columns=["indemnity_amount", "34kt_Winds", "50kt_Winds","reciprocal_indemnity_amount","threefour_pct", "fifity_pct", "64kt_Winds","log_indemnity_amount"])
datasqrt = datasqrt.iloc[:, [4, 0, 1,2,3]]


# check the correlation bewteen Y variable and the other variables in regressionof log data
sns.pairplot(datalog, kind="reg", diag_kws={"bins": "sqrt"})
plt.show()

# heatmaps
sns.heatmap(datalog.corr(), annot=True)


# check the correlation bewteen Y variable and the other variables in regressionof sqrt data
sns.pairplot(datasqrt, kind="reg", diag_kws={"bins": "sqrt"})
plt.show()

# heatmaps
sns.heatmap(datasqrt.corr(), annot=True)




# ******************************************************************************
## modeling part

import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as sm

# model1 with Log Transformation data
lm1 = sm.ols(formula="log_indemnity_amount ~ full_area + sixfour_pct + PlantedAcres + pcp", data=data1).fit()
print(lm1.summary())


## model2 with Square-Root Transformation data
lm2 = sm.ols(formula="sqrt_indemnity_amount ~ full_area + sixfour_pct + PlantedAcres + pcp", data=data1).fit()
print(lm2.summary())



# ******************************************************************************

# K-fold cross validation
# Import libraries

# Load and preprocess the dataset for Log Transformation data

X = datalog.iloc[:, [1, 2, 3, 4]].values
y = datalog.iloc[:, 0].values
# Perform data preprocessing, such as encoding categorical variables and scaling numerical variables using StandardScaler

# Split the dataset
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=0)

# Define the linear regression model
model = LinearRegression()

# Define the K-fold cross-validation
kfold = KFold(n_splits=10, shuffle=True, random_state=0)

# Train and evaluate the model
scores = cross_val_score(model, X_train, y_train,
                         cv=kfold, scoring='neg_mean_squared_error')
mse_scores = -scores  # convert negative MSE scores to positive
rmse_scores = np.sqrt(mse_scores)
print('RMSE scores:', rmse_scores)
print('Average RMSE:', np.mean(rmse_scores))

# Residuals plot
# Generate residuals plot
y_pred = model.fit(X_train, y_train).predict(X_test)
residuals = y_test - y_pred
fig, ax = plt.subplots(figsize=(10, 8))  # Set the figure size to 10x8 inches
ax.scatter(y_pred, residuals)
ax.axhline(y=0, color='r', linestyle='--')
ax.set_xlabel('Predicted values')
ax.set_ylabel('Residuals')
ax.set_title('Residuals plot')
plt.show()


# ******************************************************************************

# K-fold cross validation
# Import libraries

# Load and preprocess the dataset for sqrt Transformation data

X1 = datasqrt.iloc[:, [1, 2, 3, 4]].values
y1 = datasqrt.iloc[:, 0].values
# Perform data preprocessing, such as encoding categorical variables and scaling numerical variables using StandardScaler

# Split the dataset
X1_train, X1_test, y1_train, y1_test = train_test_split(
    X1, y1, test_size=0.2, random_state=0)

# Define the linear regression model
model = LinearRegression()

# Define the K-fold cross-validation
kfold = KFold(n_splits=10, shuffle=True, random_state=0)

# Train and evaluate the model
scores = cross_val_score(model, X1_train, y1_train,
                         cv=kfold, scoring='neg_mean_squared_error')
mse_scores = -scores  # convert negative MSE scores to positive
rmse_scores = np.sqrt(mse_scores)
print('RMSE scores:', rmse_scores)
print('Average RMSE:', np.mean(rmse_scores))

# Residuals plot
# Generate residuals plot
y1_pred = model.fit(X1_train, y1_train).predict(X1_test)
residuals = y1_test - y1_pred
fig, ax = plt.subplots(figsize=(10, 8))  # Set the figure size to 10x8 inches
ax.scatter(y1_pred, residuals)
ax.axhline(y=0, color='r', linestyle='--')
ax.set_xlabel('Predicted values')
ax.set_ylabel('Residuals')
ax.set_title('Residuals plot')
plt.show()



            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
      
            
      
            
      
            
      
            
      
            
      
            
      
            
      
            
      
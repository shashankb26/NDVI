#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 09:58:03 2022

@author: satyukt
"""

import geopandas as gpd
from shapely.geometry import Polygon
from matplotlib import colors
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
#from google.colab import drive
from glob import glob
import warnings
warnings.filterwarnings("ignore")
from gc import callbacks
import boto3
import os
import sys
import threading
from multiprocessing.pool import ThreadPool 
import time
import ast
import io
from PIL import Image
from datetime import date,datetime,timedelta
import boto3
import json
import botocore
from botocore.exceptions import ClientError
import pymysql
import json
import os
import datetime as dt
import pandas as pd
import ast
from natsort import natsorted
import ee

# Trigger the authentication flow.
#ee.Authenticate()

# Initialize the library.
ee.Initialize()


START_DATE = '2021-08-01'
END_DATE = '2022-05-30'
CLOUD_FILTER = 90
CLD_PRB_THRESH = 50
NIR_DRK_THRESH = 0.15
CLD_PRJ_DIST = 1
BUFFER = 50

def ndviS2image(image) :
  nir = image.select('B8');
  red = image.select('B4');

  ndvi =  nir.subtract(red).divide(nir.add(red)).rename('NDVI');
  return image.addBands(ndvi)

#def clipS2image(image) :
#  return image.clip(AOI)


# SCL band cloud mask
def sclCloudMask(image):
  ndvi = image.select('NDVI')
  crs = (ndvi.projection()).crs()
  scl = image.select('SCL')
  reScl = scl.resample('bilinear').reproject(crs = crs, scale = 10)
  mask = reScl.gt(3) and reScl.lte(7)
  # mask = mask.lte(7) # values which are not cloud
  maskedNdvi = ndvi.updateMask(mask)
  
  return ee.Image(maskedNdvi)

#def ndviS2image(image) :
#  nir = image.select('B8');
#  red = image.select('B4');
#
#  ndvi =  nir.subtract(red).divide(nir.add(red)).rename('NDVI');
#  return image.addBands(ndvi)

# create sample rectangles
def mapRectangle(image) :
  return ee.Image(image).sampleRectangle(region = AOI, defaultValue = float(-9999))

# reducing resolution function
def reduce_resolution(image) :
  crs = (image.select('NDVI').projection()).crs()
  return image.reproject(crs = crs, scale = scale_res)

'''dbname="sat2farmdb"
user="remote"
password="satelite"
host="34.125.75.120"
port=3306
conn = pymysql.connect(host=host, user=user,port=port,
                           passwd=password, db=dbname)

active_client_list=pd.read_csv("/home/satyukt/Downloads/ActiveUsers_MonthlyData.csv",usecols=['ClientID']).values.tolist()
active_list=[f[0] for f in active_client_list]

fill_df = pd.DataFrame()
for cl in active_list: 
#    print(cl)
    conn = pymysql.connect(host=host, user=user,port=port,passwd=password, db=dbname)
    cur = conn.cursor()
    cur.execute(f"select id,polyinfo from polygonStore where active=1 and clientID={cl}")
    
    
    
    res = cur.fetchall()
#    client_id.append(cl)
    for r in res:
#        print(r,cl)
        farm_id=r[0]
#        print(farm_id,cl)
        dict=ast.literal_eval(r[1])
        coordinates=dict['geo_json']['geometry']['coordinates'][0]
        temp_df = pd.DataFrame([{"farm_id" : farm_id, "client_id" : cl, "coordinates":coordinates}])
        fill_df = fill_df.append(temp_df)
#        fill_df[i,"farmid_list"] = (farm_id)
#        coordinate_list.append(coordinates)
#        fill_df[i,"client_id"] = (cl)
    cur.close()
    conn.close()

"""Establishing connection to S3"""

s3 = boto3.resource(
    's3',
    region_name='us-east-1',
    aws_access_key_id= "AKIAXMFIWZ7AO7DFOCVK",
    aws_secret_access_key="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"
)

s3_client = boto3.client(
                's3',
                region_name = "ap-south-1",
                aws_access_key_id="AKIAXMFIWZ7AO7DFOCVK",
                aws_secret_access_key="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"
            )

def genurl(paths):
            url = s3_client.generate_presigned_url(
                ClientMethod='get_object',
                Params={'Bucket': "data.satyukt.com", 'Key': paths, },
                ExpiresIn=180,
            )
            # url = url.replace("https://s3.ap-south-1.amazonaws.com/","http://")
            return url



class ProgressPercentage(object):

    def __init__(self, filename):
        self._filename = filename
        self._size = float(os.path.getsize(filename))
        self._seen_so_far = 0
        self._lock = threading.Lock()

    def __call__(self, bytes_amount):
        # To simplify, assume this is hooked up to a single filename
        with self._lock:
            self._seen_so_far += bytes_amount
            percentage = (self._seen_so_far / self._size) * 100
            sys.stdout.write(
                "\r%s  %s / %s  (%.2f%%)" % (
                    self._filename, self._seen_so_far, self._size,
                    percentage))
            sys.stdout.flush()

def writefile(s3path,content,Type):
    
    s3.Object('data.satyukt.com',s3path).put(Body=content,ACL = "public-read", ContentType = Type)

def uploadfile(localpath,s3path):
    s3.meta.client.upload_file( localpath, "data.satyukt.com", s3path,ExtraArgs={"ContentType":"image/png","ACL":"public-read"},Callback=ProgressPercentage(localpath))'''


# generating custom colormap
levels = [0, 0.2, 0.4, 0.6, 0.8, 1]
clrs = [(0, "#ff0000"), (0.13, "#e9331f"), (0.26, "#fe591e"), (0.39, "#dcde4c"), (0.52, "#ffff00"), (0.653846, "#94f25d"), (0.78, "#3eca6a"), (0.9, "#006d2c"), (1, "#055005")]
legend_elements = [Patch(facecolor="#055005", edgecolor='black',label='Very High'), Patch(facecolor="#3eca6a", edgecolor='black',label='High'), Patch(facecolor="#ffff00", edgecolor='black',label='Medium'), Patch(facecolor="#fe591e", edgecolor='black',label='Low'), Patch(facecolor="#ff0000", edgecolor='black',label='Very Low')]
cmap = colors.LinearSegmentedColormap.from_list('rsm', clrs, N=256)

# getting farm WKTs
#farms = glob(drive_root+'/*.csv')



def img_create(image, datee, img, fig) :
  ''' function to create images and write them in s3'''
#  client_id=254
#  url = f"micro/{farm_id}/satelite_data/SM/PNG/SM_{date}.png"
  img.set_data(image)

  # saving image in bytes to store in memory
  buf = io.BytesIO()
  fig.savefig(buf)
  buf.seek(0)
#  plt.savefig(f"/home/satyukt/Downloads/{client_id}/{farm_id}/NDVI_{datee}.png")
  
  plt.savefig(f"/home/satyukt/Projects/1259/Sitapur/Aalampur/ndvi/NDVI_{datee}")
#  print(genurl(filepath))
  

  # writing data to S3
#  writefile(filepath,buf.read(),"image/png")


# main function to compute rsm and call plotting function
def main(coordinates) :
    

  try : 
    global AOI, polygon_geom, polygon, min_ras, sensitivity, extent_gdf, polygon_bound, buff, x_axis, y_axis,scale_res,arrayPlot

    # farm_id = info['polygon_id'].loc[i]
    #farm_id=i
    # print("farm_id_is"+str(farm_id))
    # coordinates = info['coordinates'].loc[i]
    AOI = ee.Geometry.Polygon(coordinates)

    # generating farm polygon
    lon_list = [x[0] for x in coordinates]
    lat_list = [x[1] for x in coordinates]
    polygon_geom = Polygon(zip(lon_list, lat_list))
    crs = {'init': 'epsg:4326'}
    polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])  

    extent_gdf = polygon.bounds.values[0]
    
    buff = 0.00008
    x_axis = [extent_gdf[0]-buff, extent_gdf[2]+buff]
    y_axis = [extent_gdf[1]-buff, extent_gdf[3]+buff]
    
    # generating MBR for farm polygon
    poly_bound = Polygon([[extent_gdf[0],extent_gdf[1]],
                          [extent_gdf[2],extent_gdf[1]],
                          [extent_gdf[2],extent_gdf[3]],
                          [extent_gdf[0],extent_gdf[3]]])
    polygon_bound = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[poly_bound]) 
    
    # dataset collection
    dataset = (ee.ImageCollection('COPERNICUS/S2_SR')
      #.filterDate((datetime.now() - timedelta(days = 30)).strftime('%Y-%m-%d'), (datetime.now()).strftime('%Y-%m-%d'))
      .filterDate('2021-06-01', '2022-06-20')
      .filterBounds(AOI))
    
    # image collection for relative soil moisture
    #rsm = dataset.filterDate('2021-08-01', '2022-04-30') # segregate data for 6 months
    eedates = [x.split('/')[2][:8] for x in dataset.aggregate_array('system:id').getInfo()]
    # s3dates = [(object_summary.key).split('/')[-1].split('.')[0].split('_')[-1] for object_summary in myBucket.objects.filter(Prefix=f"micro/{farm_id}/satelite_data/SM/PNG/")]
    s3dates = []
    limit = [eedates.index(x) for x in eedates if x not in s3dates]

    # calculate RSM for dates not present in S3
    if len(limit) > 0 :
      scale_res = 10

      ndvi_time_series = (dataset.map(ndviS2image).map(sclCloudMask).select('NDVI').map(mapRectangle)).toList(dataset.size())
      # print(ndvi_time_series.getInfo()[0])
      try:
          arrayPlot = [np.asarray(x['properties']['NDVI']) for x in ndvi_time_series.getInfo()] # throws error for bigger farms
          
      except ee.EEException as eex:
#          print(eex)
          if "Too many pixels in sample" in str(eex) or "Response size exceeds limit " in str(eex):
              print(str(eex))
              choice = True
              while choice:
                 try:
                      scale_res += 10
                      data=dataset.map(ndviS2image).map(sclCloudMask).select('NDVI').map(reduce_resolution)
                      ndvi_time_series_re = (data.map(mapRectangle)).toList(data.size()) # sample rectangle calculation
                      arrayPlot = [np.asarray(x['properties']['NDVI']) for x in ndvi_time_series_re.getInfo()] # for bigger farms
                      len_for_axs  = scale_res*0.0001
                      buff = 0.01
                      x_axis = [extent_gdf[0]-len_for_axs, extent_gdf[2]+len_for_axs]
                      y_axis = [extent_gdf[1]-len_for_axs, extent_gdf[3]+len_for_axs]
                      choice = False
                 except ee.EEException as eex:
                     if "Too many pixels in sample" in str(eex) :
                        print(scale_res, str(eex))
                        continue
                     else:
                         break

      # creating plot handle for each farm
      fig, ax = plt.subplots(figsize=(5,5), facecolor = 'white')
      im = ax.imshow(arrayPlot[limit[0]], cmap = cmap, extent = [extent_gdf[0], extent_gdf[2], extent_gdf[1], extent_gdf[3]], vmin = 0, vmax = 1)
      (polygon.overlay(gpd.GeoDataFrame(geometry=polygon_bound.geometry.buffer(buff)), how='symmetric_difference')).plot(ax=ax, facecolor = 'white')
      polygon.plot(ax=ax, facecolor = 'none', edgecolor = 'black')
      ax.axis('off')
      cb1 = plt.colorbar(mappable = im, orientation="vertical", ticks = [round(x,1) for x in levels], fraction=0.03)
      ax.set_xlim(x_axis)
      ax.set_ylim(y_axis)
      fig.suptitle(f"Normalized Difference Vegetation Index", ha = 'center', fontsize=12, fontweight = 'bold')
      plt.legend(handles=legend_elements, loc='lower center', ncol = 5, bbox_to_anchor=(0.5,0), bbox_transform=plt.gcf().transFigure,  prop={'size': 8})
      plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

      # calling plotting function  and passing the image handle
      for j in limit :
        in_array = arrayPlot[j].astype(float)
        in_array[in_array == -9999] = np.nan
        if ~np.isnan(in_array).all():
          # print('\n',in_array)
          img_create(in_array, eedates[j], im, fig)
  
  except ee.EEException as eex :
    print(f"earth engine error for -> {eex}")
    pass
  
  except Exception as e :
    print(f"error for -> {e}")
    pass

"""# **Calling main function**"""

shp_path = f"/home/satyukt/Projects/1001/shaas/Sitapur/Aalampur.shp"
gdf_ext_shp = (gpd.read_file(shp_path)).to_crs(4326)
geom = gdf_ext_shp["geometry"] 
jsonDict = eval(geom.to_json())
coordinates = jsonDict['features'][0]['geometry']['coordinates'][0]
main(coordinates)
        
  # except ee.EEException as eex :
  #   print(f"earth engine error for {farm_id} -> {eex}, reducing resolution and trying again")
    # if "Too many pixels in sample" in str(eex) :
    #   rsm_re = rsm.map(reduce_resolution)
    #   rsm_time_series_re = (rsm_re.map(mapRectangle)).toList(rsm_re.size()) # sample rectangle calculation
    #   arrayPlot = [np.asarray(x['properties']['VV']) for x in rsm_time_series_re.getInfo()] # for bigger farms
    #   for i in arrayPlot :
    #     print(i)
  
  # except ValueError as ve :
  #   print(f"value error for {farm_id} -> {ve}")
  #   pass

# farm_index = list(range(len(info)))

#for i in fill_df['farm_id']:                              
#  print(i)  
#  coordinates=fill_df.loc[fill_df['farm_id'] == i, 'coordinates'].iloc[0]
#  client_id=fill_df.loc[fill_df['farm_id'] == i, 'client_id'].iloc[0]
#  main(i,coordinates,client_id)
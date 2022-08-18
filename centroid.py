#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 10:53:45 2022

@author: satyukt
"""

import pandas as pd
import ast
import json
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

in_csv = pd.read_csv("/media/edrive1/Shashank/Shashank_Projects/UP_classification/csv/NDVI_Sitapur_Kharif.csv", usecols = ['Position'])
in_csv.head()

coordinates = ast.literal_eval(in_csv['Position'][2])['coordinates'][0]

poly = Polygon(coordinates)
center = [poly.centroid.x, poly.centroid.y]
print(center)
plt.plot(*poly.exterior.xy)





final_list = []


for i in range(len(in_csv)):
    try:
        coordinates = ast.literal_eval(in_csv['Position'][i])['coordinates'][0]
        poly = Polygon(coordinates)
        center = [poly.centroid.x, poly.centroid.y]
        final_list.append('{ "coordinates" : '+ str(center) + '}')
    except:
        final_list.append('{ "coordinates" : '+ str([0, 0]) + '}')
#        print(i)


df = pd.DataFrame(final_list)
df.to_csv("/home/satyukt/Downloads/pragya.csv")
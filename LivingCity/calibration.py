import time
import os, zipfile, requests, pandas as pd
import geopandas as gpd
import osmnx as ox, networkx as nx
import ast
import statistics
import numpy as np
from skopt import gp_minimize
from skopt.plots import plot_convergence
#from sklearn.neighbors import BallTree
from shapely.geometry import Point
import random
from matplotlib import pyplot as plt
import matplotlib.path as mpltPath
import matplotlib.cm as cm
import matplotlib.colors as colors
import json 
from ast import literal_eval
import itertools
from collections import Counter
import descartes
import matplotlib
from shapely.geometry import LineString # To create line geometries that can be used in a GeoDataFrame
import math
import folium # To generate a Leaflet-based map of my data throughout my analysis
import branca.colormap as cm
import time
from time import gmtime, strftime
from selenium import webdriver
import imageio
import imgkit
import requests
import calendar
import seaborn as sns
import sys
#sys.path.insert(1, '../google_data/baytraffic_share/credentials')
#from google_key import GOOGLE_MAPS_API_KEY
import subprocess

from keplergl import KeplerGl 
from mapboxgl.utils import *
from mapboxgl.viz import *

# Load an empty map
from geopandas import GeoDataFrame # To create a GeoDataFrame from a DataFrame
from shapely.geometry import LineString # To create line geometries that can be used in a GeoDataFrame
pd.options.display.max_colwidth = 100



print('ox {}\nnx {}'.format(ox.__version__, nx.__version__))
start_time = time.time()


def load_network(nodes_file, edges_file, path):
    nodes = pd.read_csv("{}/{}".format(path, nodes_file))
    edges = pd.read_csv("{}/{}".format(path, edges_file))
    return nodes, edges


def create_per_edge_speeds_df(num_printouts, path, edges_df, nodes_df, remove_error=False, save=False):
    edge_vels_df = pd.read_csv('{}/all_edges_vel_0.txt'.format(path), names=['time_0'], header=None)
    for ind in range(1,num_printouts):
        edge_vels_df['time_{}'.format(ind)] = pd.read_csv('{}/all_edges_vel_{}.txt'.format(path, ind))
    
    if remove_error:
        edges_error = pd.read_csv('{}/edges_error.txt'.format(path), names=['index'], header=None)
        edges_error_vals = edges_error['index'].tolist()
        filtered_edges = edges_df.drop(edges_df.index[edges_error_vals])
    
    #merge osmids, lats, longs from     
    edge_vels_df = edge_vels_df.join(filtered_edges, lsuffix='_caller', rsuffix='_other')
    edge_vels_df_node1 = pd.merge(edge_vels_df, nodes_df[['osmid', 'x', 'y']], left_on='u', right_on='osmid', how='left')
    edge_vels_df_node1_node2 = pd.merge(edge_vels_df_node1, nodes_df[['osmid', 'x', 'y']], left_on='v', right_on='osmid', how='left', suffixes=['_u', '_v'])
    
    df = edge_vels_df_node1_node2
    if save:
        df['geometry'] = df.apply(lambda row: 'LINESTRING ({} {}, {} {})'.format(row['x_u'], row['y_u'], row['x_v'], row['y_v']), axis=1)
        df.to_csv('{}/edges_over_time.csv'.format(path), index=False)
    else:
        df['geometry'] = df.apply(lambda row: LineString([(row['x_u'], row['y_u']), (row['x_v'], row['y_v'])]), axis=1)
        
    return df

def merge_edges_with_uber_edges(edges_df, uber_df):
    new_df = pd.merge(edges_df, uber_df,  how='left', left_on=['u','v'], right_on = ['osm_start_node_id','osm_end_node_id'])
    return new_df

def merge_edges_with_times_with_uber_edges(edge_vels_df, same_edges_df):
    new_df = pd.merge(edge_vels_df, same_edges_df,  how='left', on=['u','v'])
    return new_df

def write_options_file(a_prev, b_prev, T_prev, s_0_prev, a, b, T, s_0):
    fin = open("command_line_options.ini", "r+")

    # Read in the file
    with open('command_line_options.ini', 'r') as file :
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('A={}'.format(a_prev), 'A={}'.format(a)).replace('B={}'.format(b_prev), 'B={}'.format(b)).replace('T={}'.format(T_prev), 'T={}'.format(T)).replace('s_0={}'.format(s_0_prev), 's_0={}'.format(s_0))

    # Write the file out again
    with open('command_line_options.ini', 'w') as file:
        file.write(filedata)

    #for line in fin:
    #    fin.write(line.replace('A={}'.format(a_prev), '{}'.format(a)))
    #    fin.write(line.replace('B={}'.format(b_prev), '{}'.format(b)))
    #    fin.write(line.replace('T={}'.format(T_prev), '{}'.format(T)))
    #fin.close()

def load_edges_u_v(path):
    edges_u = pd.read_csv('{}/edges_u.txt'.format(path), names=['index'], header=None)
    edges_v = pd.read_csv('{}/edges_v.txt'.format(path), names=['index'], header=None)
    edges_u = edges_u.astype(int)
    edges_v = edges_v.astype(int)
    edges_u = edges_u['index'].tolist()
    edges_v = edges_v['index'].tolist()
        
    return edges_u, edges_v

def create_network_from_edges_node_ids(edges_u, edges_v, nodes, edges, path, save=True):
    num_printouts = 7
    new_df = pd.DataFrame({'u': edges_u, 'v': edges_v})
    route_edge_node1 = pd.merge(new_df, nodes[['osmid', 'x', 'y']], left_on='u', right_on='osmid', how='left')
    route_edge_node1_node2 = pd.merge(route_edge_node1, nodes[['osmid', 'x', 'y']], left_on='v', right_on='osmid', how='left', suffixes=['_u', '_v'])
    #route_edge_node1_node2 = pd.merge(route_edge_node1_node2, edges, left_on='u', right_on='osmid_u', how='left')
    #route_edge_node1_node2 = pd.merge(route_edge_node1_node2, edges, left_on='v', right_on='osmid_v', how='left')

    df = route_edge_node1_node2
    if save:
        df['geometry'] = df.apply(lambda row: 'LINESTRING ({} {}, {} {})'.format(row['x_u'], row['y_u'], row['x_v'], row['y_v']), axis=1)
        df.to_csv('{}/edges_over_time.csv'.format(path), index=False)
    else:
        df['geometry'] = df.apply(lambda row: LineString([(row['x_u'], row['y_u']), (row['x_v'], row['y_v'])]), axis=1)
        
    #new_df = pd.merge(df, edges,  how='left', left_on=['osmid_u','osmid_v'], right_on = ['u','v'])
    
    new_df['time_0'] = pd.read_csv('{}/all_edges_vel_0.txt'.format(path), names=['time_0'], header=None)
    for ind in range(1,num_printouts):
        new_df['time_{}'.format(ind)] = pd.read_csv('{}/all_edges_vel_{}.txt'.format(path, ind)) * 2.23694
        
    new_df['microsim_avg'] = new_df[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)
    
    #if remove_error:
    #    edges_error = pd.read_csv('{}/edges_error.txt'.format(path), names=['index'], header=None)
    #    edges_error_vals = edges_error['index'].tolist()
    #    filtered_edges = edges_df.drop(edges_df.index[edges_error_vals])

    
    print(new_df.shape[0])
    new_df = new_df.replace(0, pd.np.nan)
    new_df = new_df.dropna()
    return new_df

def compare_to_uber_new(path, start_time=5, end_time=12, one_time_slice=False):

    nodes, edges = load_network("nodes.csv", "edges.csv", "../new_full_network_full_demand")
    #edge_vels_df = create_per_edge_speeds_df(7, "{}".format(path), edges, nodes, remove_error=False, save=True)
    edges_u, edges_v = load_edges_u_v(path)
    edge_vels_df = create_network_from_edges_node_ids(edges_u, edges_v, nodes, edges, path)
    #edge_vels_df = edge_vels_df.replace(0, pd.np.nan)
    #print(edge_vels_df.shape[0])
    #edge_vels_df = edge_vels_df.dropna()
    #print(edge_vels_df.shape[0])
    #edge_vels_df = edge_vels_df.loc[~(edge_vels_df==0).all(axis=1)]
    #print(edge_vels_df.head(50))
    #load Uber speeds
    uber_speeds_df = pd.read_csv("berkeley_2018/uber_data/movement-speeds-quarterly-by-hod-san-francisco-2019-Q2.csv")

    #filter Uber data to 5am to 12pm
    uber_speeds_5_to_12 = uber_speeds_df[(uber_speeds_df['hour_of_day'] >= start_time) & (uber_speeds_df['hour_of_day'] < end_time)]
    print("number of trips from {} to {} = {} in actual Uber data".format(start_time, end_time, uber_speeds_5_to_12.shape[0]))

    # Make new df that calculates the avg. speed across the 7 hours for every edge (as one row as opposed to 7 different rows for the same edge)
    new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_mean']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
    #new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_p85']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
    new_uber = new_uber.reset_index(level=['osm_start_node_id', 'osm_end_node_id'])

    # Merge our network edges with Uber edges
    merged_edges = merge_edges_with_uber_edges(edges, new_uber)
    merged_edges = merged_edges.dropna()

    # Merge edge vels (with x timesteps) with the merged Uber data
    merged_edges = merge_edges_with_times_with_uber_edges(edge_vels_df, merged_edges)
    merged_edges = merged_edges.dropna()
    print("number of trips from {} to {} = {}".format(start_time, end_time, merged_edges.shape[0]))

    # Average the speeds over all the timesteps from the microsim and add column to the df
    merged_time_edges = merged_edges.copy()
    if one_time_slice:
        merged_time_edges['microsim_avg'] = merged_time_edges['time_{}'.format(start_time - 5)]
    else: 
        merged_time_edges['microsim_avg'] = merged_time_edges[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)
    #merged_time_edges['microsim_avg'] = merged_time_edges[['time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)

    # Drop unnecessary columns
    """
    merged_time_edges = merged_time_edges.drop(['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6',
                                                    'uniqueid_x' , 'osmid_u',
                                                    'x_u', 'y_u', 'osmid_v', 'x_v', 'y_v', 'uniqueid_y',
                                                    'length_y', 'lanes_y', 'speed_mph_y', 'osm_start_node_id', 'osm_end_node_id'], axis=1)
    """
    # Add *difference* column to show delta between our microsim and Uber speeds
    merged_time_edges['diff'] = merged_time_edges['speed_mph_mean'] - merged_time_edges['microsim_avg']
    #merged_time_edges['diff'] = abs(merged_time_edges['speed_mph_p85'] - merged_time_edges['microsim_avg'])
    merged_time_edges['microsim_avg'] = merged_time_edges[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)
    merged_time_edges['speed_limit_vs_microsim'] = merged_time_edges['speed_mph'] - merged_time_edges['microsim_avg']
    print("merged with uber edges length = {}".format(merged_time_edges.shape[0]))
    print("uber avg = {}".format(merged_time_edges['speed_mph_mean'].mean()))
    print("microsim avg = {}".format(merged_time_edges['microsim_avg'].mean()))

    merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_mean'] / merged_time_edges['microsim_avg']
    #merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_p85'] / merged_time_edges['speed_mph']
    #merged_time_edges['uber_speed_limit_ratio'] = merged_time_edges['speed_mph_p85'] / merged_time_edges['speed_mph_x']
    print("average uber / microsim ratio = {}".format(merged_time_edges['uber_microsim_ratio'].mean()))
    #print("average uber / speed_limit ratio = {}".format(merged_time_edges['uber_speed_limit_ratio'].mean()))


    #print(merged_time_edges.head())
    print("diff mean = {}".format(merged_time_edges['diff'].mean()))

    #calculate RMSE
    rmse = np.sqrt(sum(merged_time_edges['diff']**2) / merged_time_edges.shape[0])
    print("RMSE = {}".format(rmse))
    return merged_time_edges


def calibrate_new(a_b_T_s_0_vec, path, start_time=5, end_time=12, one_time_slice=False):
    a_prev = 0.8
    b_prev = 0.8
    T_prev = 1.5
    s_0_prev = 7.0

    diff_vec = []
    



    for x in a_b_T_s_0_vec:
        a = x[0]
        b = x[1]
        T = x[2]
        s_0 = x[3]

        write_options_file(a_prev, b_prev, T_prev, s_0_prev, a, b, T, s_0)
        print(a, b, T, s_0)

        subprocess.run(["./LivingCity", "&"])

        #load street network and edge/node data
        nodes, edges = load_network("nodes.csv", "edges.csv", "berkeley_2018/new_full_network")
        edges_u, edges_v = load_edges_u_v(path)
        #get edge speed data from microsim
        edge_vels_df = create_network_from_edges_node_ids(edges_u, edges_v, nodes, edges, path)
    
        #filter Uber data to 5am to 12pm
        uber_speeds_df = pd.read_csv("berkeley_2018/uber_data/movement-speeds-quarterly-by-hod-san-francisco-2019-Q2.csv")
        uber_speeds_5_to_12 = uber_speeds_df[(uber_speeds_df['hour_of_day'] >= start_time) & (uber_speeds_df['hour_of_day'] < end_time)]
        new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_mean']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
        #new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_p85']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
        new_uber = new_uber.reset_index(level=['osm_start_node_id', 'osm_end_node_id'])
        print("number of trips from {} to {} = {} in actual Uber data".format(start_time, end_time, uber_speeds_5_to_12.shape[0]))


        # Merge our network edges with Uber edges
        merged_edges = merge_edges_with_uber_edges(edges, new_uber)
        merged_edges = merged_edges.dropna()

        # Merge edge vels (with x timesteps) with the merged Uber data
        merged_edges = merge_edges_with_times_with_uber_edges(edge_vels_df, merged_edges)
        merged_edges = merged_edges.dropna()
        print("number of trips from {} to {} = {}".format(start_time, end_time, merged_edges.shape[0]))

        # Average the speeds over all the timesteps from the microsim and add column to the df
        merged_time_edges = merged_edges.copy()
        if one_time_slice:
            merged_time_edges['microsim_avg'] = merged_time_edges['time_{}'.format(start_time - 5)]
        else: 
            merged_time_edges['microsim_avg'] = merged_time_edges[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)
        #merged_time_edges['microsim_avg'] = merged_time_edges[['time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)

        # Drop unnecessary columns
        """
        merged_time_edges = merged_time_edges.drop(['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6',
                                                    'uniqueid_x' , 'osmid_u',
                                                    'x_u', 'y_u', 'osmid_v', 'x_v', 'y_v', 'uniqueid_y',
                                                    'length_y', 'lanes_y', 'speed_mph_y', 'osm_start_node_id', 'osm_end_node_id'], axis=1)
        """
        # Add *difference* column to show delta between our microsim and Uber speeds
        merged_time_edges['diff'] = merged_time_edges['speed_mph_mean'] - merged_time_edges['microsim_avg']
        #merged_time_edges['diff'] = abs(merged_time_edges['speed_mph_p85'] - merged_time_edges['microsim_avg'])
        merged_time_edges['microsim_avg'] = merged_time_edges[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)
        merged_time_edges['speed_limit_vs_microsim'] = merged_time_edges['speed_mph'] - merged_time_edges['microsim_avg']
        print("merged with uber edges length = {}".format(merged_time_edges.shape[0]))
        print("uber avg = {}".format(merged_time_edges['speed_mph_mean'].mean()))
        print("microsim avg = {}".format(merged_time_edges['microsim_avg'].mean()))

        merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_mean'] / merged_time_edges['microsim_avg']
        #merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_p85'] / merged_time_edges['speed_mph']
        #merged_time_edges['uber_speed_limit_ratio'] = merged_time_edges['speed_mph_p85'] / merged_time_edges['speed_mph_x']
        print("average uber / microsim ratio = {}".format(merged_time_edges['uber_microsim_ratio'].mean()))
        #print("average uber / speed_limit ratio = {}".format(merged_time_edges['uber_speed_limit_ratio'].mean()))


        #calculate RMSE
        #rmse = np.sqrt(sum(merged_time_edges['diff']**2) / merged_time_edges.shape[0])
        diff = abs(merged_time_edges['diff'].mean())
        #print(diff)
        diff_vec += [diff,]
        print("diff mean = {}".format(diff))

        #calculate RMSE
        rmse = np.sqrt(sum(merged_time_edges['diff']**2) / merged_time_edges.shape[0])
        print("RMSE = {}".format(rmse))

        T_prev = T
        b_prev = b
        a_prev = a
        s_0_prev = s_0
        print("diff vector = {}".format(diff_vec))

    min_diff = min(diff_vec)
    index = diff_vec.index(min_diff)
    return diff_vec, min_diff, index, a_prev, b_prev, T_prev, s_0_prev




def calibrate(a_b_T_vec):
    
    a_prev = 0.8
    b_prev = 0.8
    T_prev = 1.5

    diff_vec = []
    for x in a_b_T_vec:
        a = x[0]
        b = x[1]
        T = x[2]

        write_options_file(a_prev, b_prev, T_prev, a, b, T)
        print(a, b, T)

        subprocess.run(["./LivingCity", "&"])

        #load street network
        nodes, edges = load_network("nodes.csv", "edges.csv", "berkeley_2018/new_full_network")

        #create dataframe with every edge's speed for each timestep
        edge_vels_df = create_per_edge_speeds_df(7, ".", edges, nodes, remove_error=True, save=False)

        #load Uber speeds
        uber_speeds_df = pd.read_csv("berkeley_2018/uber_data/movement-speeds-quarterly-by-hod-san-francisco-2019-Q2.csv")

        #filter Uber data to 5am to 12pm
        uber_speeds_5_to_12 = uber_speeds_df[(uber_speeds_df['hour_of_day'] >= 5) & (uber_speeds_df['hour_of_day'] < 12)]

        # Make new df that calculates the avg. speed across the 7 hours for every edge (as one row as opposed to 7 different rows for the same edge)
        new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_mean']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
        new_uber = new_uber.reset_index(level=['osm_start_node_id', 'osm_end_node_id'])

        # Merge our network edges with Uber edges
        merged_edges = merge_edges_with_uber_edges(edges, new_uber)
        merged_edges = merged_edges.dropna()

        # Merge edge vels (with x timesteps) with the merged Uber data
        merged_edges = merge_edges_with_times_with_uber_edges(edge_vels_df, merged_edges)
        merged_edges = merged_edges.dropna()

        # Average the speeds over all the timesteps from the microsim and add column to the df
        merged_time_edges = merged_edges.copy()
        merged_time_edges['microsim_avg'] = merged_time_edges[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)

        # Drop unnecessary columns
        merged_time_edges = merged_time_edges.drop(['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6',
                                                    'uniqueid_x' , 'osmid_u',
                                                    'x_u', 'y_u', 'osmid_v', 'x_v', 'y_v', 'uniqueid_y',
                                                    'length_y', 'lanes_y', 'speed_mph_y', 'osm_start_node_id', 'osm_end_node_id'], axis=1)

        # Add *difference* column to show delta between our microsim and Uber speeds
        merged_time_edges['diff'] = merged_time_edges['speed_mph_mean'] - merged_time_edges['microsim_avg']


        #print(merged_time_edges.head())

        #calculate RMSE
        #rmse = np.sqrt(sum(merged_time_edges['diff']**2) / merged_time_edges.shape[0])
        diff = merged_time_edges['diff'].mean()
        #print(diff)
        diff_vec += [diff,]
                
        T_prev = T
        b_prev = b
        a_prev = a
        print("diff vector = {}".format(diff_vec))

    min_diff = min(diff_vec)
    index = diff_vec.index(min_diff)
    return diff_vec, min_diff, index


def gradient_descent(epsilon=0.02):
    # epsilon = max delta (in mph) that we're willing to have our RMSE
    flag = True
    iteration = 0
    diff_vec = []
    min_diff = 0
    min_index = 0
    while flag:

        if iteration == 0:
            #a_vec = [1 + np.random.uniform(0,1), 6 + np.random.uniform(0,1)]
            #b_vec = [1 + np.random.uniform(0,1), 6 + np.random.uniform(0,1)]
            #T_vec = [.8 + np.random.uniform(0,1), 2 + np.random.uniform(0,1)]
            a_vec = np.random.uniform(low=1, high=10, size=(5,)) 
            b_vec = np.random.uniform(low=1, high=10, size=(5,)) 
            T_vec = np.random.uniform(low=0.1, high=2.0, size=(5,)) 
            s_0_vec = np.random.uniform(low=1, high=5.0, size=(5,)) 
        else:
            write_options_file(a_prev, b_prev, T_prev, s_0_prev, 0.8, 0.8, 1.5, 7.0)
            a_vec = a_b_T_s_0_vec[min_index][0] + np.random.uniform(-1, 1, size=(5,))
            b_vec = a_b_T_s_0_vec[min_index][1] + np.random.uniform(-1, 1, size=(5,))
            T_vec = a_b_T_s_0_vec[min_index][2] + np.random.uniform(-1, 1, size=(5,))
            s_0_vec = a_b_T_s_0_vec[min_index][3] + np.random.uniform(-1, 1, size=(5,))

        a_b_T_s_0_vec = [x for x in zip(a_vec, b_vec, T_vec, s_0_vec)]
        print("ITERATION #{}\n".format(iteration))
        print("a_b_T_s_0_vec = {}".format(a_b_T_s_0_vec))

        #diff_vec, min_diff, min_index = calibrate(a_b_T_vec)
        diff_vec, min_diff, min_index, a_prev, b_prev, T_prev, s_0_prev = calibrate_new(a_b_T_s_0_vec, ".")

        

        if abs(min_diff) < epsilon:
            flag = False
            print("Found the best model parameters a, b, T, and s_0!")
            print("a = {}, b = {}, T = {}, s_0 = {}".format(a_b_T_s_0_vec[min_index][0], a_b_T_s_0_vec[min_index][1], a_b_T_s_0_vec[min_index][2], a_b_T_s_0_vec[min_index][3]))
            return a_b_T_s_0_vec[min_index][0], a_b_T_s_0_vec[min_index][1], a_b_T_s_0_vec[min_index][2], a_b_T_s_0_vec[min_index][3]
        iteration += 1


def compare_to_uber(path):

    nodes, edges = load_network("nodes.csv", "edges.csv", "berkeley_2018/new_full_network")
    edge_vels_df = create_per_edge_speeds_df(7, "{}".format(path), edges, nodes, remove_error=True, save=False)
    edge_vels_df = edge_vels_df.replace(0, pd.np.nan)
    print(edge_vels_df.shape[0])
    edge_vels_df = edge_vels_df.dropna()
    print(edge_vels_df.shape[0])
    #edge_vels_df = edge_vels_df.loc[~(edge_vels_df==0).all(axis=1)]
    #print(edge_vels_df.head(50))
    #load Uber speeds
    uber_speeds_df = pd.read_csv("berkeley_2018/uber_data/movement-speeds-quarterly-by-hod-san-francisco-2019-Q2.csv")

    #filter Uber data to 5am to 12pm
    uber_speeds_5_to_12 = uber_speeds_df[(uber_speeds_df['hour_of_day'] >= 5) & (uber_speeds_df['hour_of_day'] < 12)]

    # Make new df that calculates the avg. speed across the 7 hours for every edge (as one row as opposed to 7 different rows for the same edge)
    new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_mean']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
    #new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_p85']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
    new_uber = new_uber.reset_index(level=['osm_start_node_id', 'osm_end_node_id'])

    # Merge our network edges with Uber edges
    merged_edges = merge_edges_with_uber_edges(edges, new_uber)
    merged_edges = merged_edges.dropna()

    # Merge edge vels (with x timesteps) with the merged Uber data
    merged_edges = merge_edges_with_times_with_uber_edges(edge_vels_df, merged_edges)
    merged_edges = merged_edges.dropna()

    # Average the speeds over all the timesteps from the microsim and add column to the df
    merged_time_edges = merged_edges.copy()
    merged_time_edges['microsim_avg'] = merged_time_edges[['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6']].mean(axis=1)

    # Drop unnecessary columns
    merged_time_edges = merged_time_edges.drop(['time_0', 'time_1', 'time_2', 'time_3', 'time_4', 'time_5', 'time_6',
                                                    'uniqueid_x' , 'osmid_u',
                                                    'x_u', 'y_u', 'osmid_v', 'x_v', 'y_v', 'uniqueid_y',
                                                    'length_y', 'lanes_y', 'speed_mph_y', 'osm_start_node_id', 'osm_end_node_id'], axis=1)

    # Add *difference* column to show delta between our microsim and Uber speeds
    merged_time_edges['diff'] = merged_time_edges['speed_mph_mean'] - merged_time_edges['microsim_avg']
    #merged_time_edges['diff'] = abs(merged_time_edges['speed_mph_p85'] - merged_time_edges['microsim_avg'])

    merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_mean'] / merged_time_edges['microsim_avg']
    #merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_p85'] / merged_time_edges['speed_mph']
    #merged_time_edges['uber_speed_limit_ratio'] = merged_time_edges['speed_mph_p85'] / merged_time_edges['speed_mph_x']
    print("average uber / microsim ratio = {}".format(merged_time_edges['uber_microsim_ratio'].mean()))
    #print("average uber / speed_limit ratio = {}".format(merged_time_edges['uber_speed_limit_ratio'].mean()))


    print(merged_time_edges.head())
    print("diff mean = {}".format(merged_time_edges['diff'].mean()))

    #calculate RMSE
    rmse = np.sqrt(sum(merged_time_edges['diff']**2) / merged_time_edges.shape[0])
    print("RMSE = {}".format(rmse))

    

if __name__== "__main__":

    gradient_descent()


    #compare_to_uber("high_acceleration_with_1_point_19_multiplier")
    #compare_to_uber("high_acceleration")



    





                






    """
    # Make a new column for every time step as opposed to each row having its own time step. 
    ts_uber = uber_speeds_5_to_12.copy()
    pivoted_ts_uber = ts_uber.pivot_table(index=['osm_start_node_id', 'osm_end_node_id'], columns='hour_of_day', values='speed_mph_mean')
    pivoted_ts_uber = pivoted_ts_uber.reset_index(level=['osm_start_node_id', 'osm_end_node_id'])
    pivoted_ts_uber = pivoted_ts_uber.dropna()

    # Merge our network edges with Uber timestep edges
    pivoted_ts_merged_edges = merge_edges_with_uber_edges(edges, pivoted_ts_uber)
    pivoted_ts_merged_edges = pivoted_ts_merged_edges.dropna()

    # Create Uber average
    pivoted_ts_merged_edges['uber_avg'] = pivoted_ts_merged_edges[[5, 6, 7, 8, 9, 10, 11]].mean(axis=1)
    """
    
  

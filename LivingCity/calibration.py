import os, sys
import time
import math
import pandas as pd
import numpy as np
#from sklearn.neighbors import BallTree
from shapely.geometry import LineString # To create line geometries that can be used in a GeoDataFrame
import time
import json
#sys.path.insert(1, '../google_data/baytraffic_share/credentials')
#from google_key import GOOGLE_MAPS_API_KEY
import subprocess

from shapely.geometry import LineString # To create line geometries that can be used in a GeoDataFrame
pd.options.display.max_colwidth = 100

from datetime import datetime

import re
from termcolor import colored
from pdb import set_trace as st


#start_time = time.time()


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

def write_options_file(params):
    fin = open("command_line_options.ini", "r+")

    # Read in the file
    with open('command_line_options.ini', 'r') as file :
        filedata = file.read()

    # Replace the target string
    for (parameter_name, parameter_value) in params.items():
        filedata = re.sub('{}=(-|(0-9)|.)*\n'.format(parameter_name), \
                          "{}={}\n".format(parameter_name,parameter_value), \
                          filedata)


    # Write the file out again
    with open('command_line_options.ini', 'w') as file:
        file.write(filedata)

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


class benchmark():
    def __init__(self, task_name):
        self.task_name = task_name
        self.start_time_benchmark = time.time()
        print(colored("Starting task '{}'".format(task_name), 'green'))
    
    def end(self):
        elapsed_time = time.time() - self.start_time_benchmark
        print(colored("Ended task '{}'. Time elapsed: ~{} mins"\
            .format(self.task_name, round(elapsed_time/60,4)),'green'))


def calibrate_new(param_list, path, new_uber, start_time, end_time, one_time_slice=False):
    min_diff = None
    min_params = None

    for params in param_list:
        a = params["a"]
        b = params["b"]
        T = params["T"]
        s_0 = params["s_0"]

        write_options_file(params)

        print("Running LivingCity with params [a:{}, b:{}, T:{}, s_0:{}]".\
                        format(params["a"], params["b"], params["T"], params["s_0"]))
        subprocess.run(["./LivingCity", "&"])


        edges_u, edges_v = load_edges_u_v(path)
        
        #load street network and edge/node data
        nodes, edges = load_network("nodes.csv", "edges.csv", "berkeley_2018/new_full_network")

        #get edge speed data from microsim
        edge_vels_df = create_network_from_edges_node_ids(edges_u, edges_v, nodes, edges, path)

        # Merge our network edges with Uber edges
        merged_edges = merge_edges_with_uber_edges(edges, new_uber)

        merged_edges = merged_edges.dropna()

        # Merge edge vels (with x timesteps) with the merged Uber data
        merged_edges = merge_edges_with_times_with_uber_edges(edge_vels_df, merged_edges)

        merged_edges = merged_edges.dropna()
        #print("number of trips from {} to {} = {}".format(start_time, end_time, merged_edges.shape[0]))

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
        #print("merged with uber edges length = {}".format(merged_time_edges.shape[0]))
        print("uber avg = {}".format(merged_time_edges['speed_mph_mean'].mean()))
        print("microsim avg = {}".format(merged_time_edges['microsim_avg'].mean()))

        merged_time_edges['uber_microsim_ratio'] = merged_time_edges['speed_mph_mean'] / merged_time_edges['microsim_avg']
        #print("average uber / microsim ratio = {}".format(merged_time_edges['uber_microsim_ratio'].mean()))

        #calculate RMSE
        diff = abs(merged_time_edges['diff'].mean())

        #calculate RMSE
        rmse = np.sqrt(sum(merged_time_edges['diff']**2) / merged_time_edges.shape[0])

        print("For [a:{}, b:{}, T:{}, s_0:{}] we have diff:{} and RMSE:{}".format(params["a"], \
                    params["b"], params["T"], params["s_0"], diff, rmse))
        
        if (min_diff == None or diff < min_diff):
            print("Found new minimum.")
            min_diff = diff
            min_params = params.copy()

        sys.stdout.flush()
    
    return min_diff, min_params


def generate_next_step_parameters(params, decay, learning_rate_params, lower_bound_params = None):
    # Based on the current parameters, we return a possibility for the next step's parameters
    
    upper_bound_params = {"a":10, "b":10, "T":math.inf, "s_0":math.inf}
    next_step_params = {}

    for param in params.keys():
        if not param in lower_bound_params or lower_bound_params[param] < 0:
            lower_bound_params[param] = 0

    for param in params.keys():
        found_next_step_in_bounds = False
        while not found_next_step_in_bounds:
            step = np.random.uniform(-1*learning_rate_params[param]*decay, learning_rate_params[param]*decay)
            next_step_params[param] = params[param] + step
            
            if next_step_params[param] > lower_bound_params[param] and next_step_params[param] < upper_bound_params[param]:
                found_next_step_in_bounds = True
    
    
    return next_step_params.copy()

def first_step_parameters():
    first_step_params = {"a": np.random.uniform(low=1, high=10), \
                         "b": np.random.uniform(low=1, high=10), \
                         "T": np.random.uniform(low=0.1, high=2.0), \
                         "s_0": np.random.uniform(low=1, high=5.0)}

    return first_step_params.copy()

def gradient_descent(epsilon=0.04, learning_rate_params = None, progress_filename = None, lower_bound_params = None, start_time=5, end_time=12):
    print("Using progress file {}".format(progress_filename))
    print("Using lower bound for parameters {}".format(lower_bound_params))
    # epsilon = max delta (in mph) that we're willing to have our RMSE

    #filter Uber data to 5am to 12pm
    uber_speeds_df = pd.read_csv("berkeley_2018/uber_data/movement-speeds-quarterly-by-hod-san-francisco-2019-Q2.csv")
    uber_speeds_5_to_12 = uber_speeds_df[(uber_speeds_df['hour_of_day'] >= start_time) & (uber_speeds_df['hour_of_day'] < end_time)]
    new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_mean']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
    #new_uber = uber_speeds_5_to_12[['osm_start_node_id','osm_end_node_id','speed_mph_p85']].groupby(['osm_start_node_id','osm_end_node_id']).mean()
    new_uber = new_uber.reset_index(level=['osm_start_node_id', 'osm_end_node_id'])
    print("number of trips from {} to {} = {} in actual Uber data".format(start_time, end_time, uber_speeds_5_to_12.shape[0]))

    progress = []

    iteration = 0
    prev_params = {"a": None, "b": None, "T": None, "s_0": None}
    current_params = {"a": None, "b": None, "T": None, "s_0": None}
    prev_diff = None

    default_learning_rate_params = {"a":1, "b":1, "T": 0.5, "s_0": 0.1}
    for param in default_learning_rate_params.keys():
        if not param in learning_rate_params:
            learning_rate_params[param] = default_learning_rate_params[param]
    print("Learning rate for parameters: {}".format(learning_rate_params))
    
    if progress_filename and os.path.isfile('{}.csv'.format(progress_filename)):
        print("Loading saved progress from {}...".format(progress_filename))
        with open('{}.json'.format(progress_filename), 'r') as fp:
            loaded_progress = json.load(fp)
        iteration = loaded_progress["number_of_iterations"]+1
        progress = loaded_progress["progress"]
        current_params = progress[-1]["params"].copy()
        current_diff = progress[-1]["current_diff"]
        prev_diff = current_diff
    

    while True:
        if iteration == 0 and not os.path.isfile('{}.csv'.format(progress_filename)):
            possible_next_params = [first_step_parameters() for _ in range(5)]
            current_diff = None
        else:
            # we store the previous values
            prev_params = current_params.copy()
            if len(progress) <= 7:
                decay = 1 # no decay
            elif len(progress) > 7:
                decay = 0.1
            elif len(progress) > 12:
                decay = 0.01
            elif len(progress) > 15:
                decay = 0.001
            possible_next_params = [generate_next_step_parameters(current_params, decay, learning_rate_params=learning_rate_params, lower_bound_params=lower_bound_params) for _ in range(5)]

        print(colored("----------------- ITERATION #{}\n".format(iteration), "cyan"))
        print("Current progress: {}" .format(progress))
        print("Current parameters: {}".format(current_params))
        print("Current diff: {}".format(current_diff))
        
        print("Next steps that will be tested on this iteration:")
        for one_possible_step in possible_next_params:
            print(one_possible_step)
        prev_params = current_params.copy()
        current_diff, current_params = calibrate_new(possible_next_params, ".", new_uber, start_time, end_time)

        if (prev_diff != None and current_diff >= prev_diff):
            print("Not moving in that direction since current_diff={} and prev_diff={}".format(current_diff, prev_diff))
            current_params = prev_params.copy()
            current_diff = prev_diff
        else:
            prev_diff = current_diff
            progress.append({"params":current_params, "current_diff":current_diff})

        if progress_filename:
            print("Saving progress at {}...".format(progress_filename))
            saved_progress = {}
            saved_progress["number_of_iterations"] = iteration
            saved_progress["progress"] = progress
            with open('{}.json'.format(progress_filename), 'w') as fp:
                json.dump(saved_progress, fp)

        

        if abs(current_diff) < epsilon:
            print("Found the best model parameters a, b, T, and s_0!")
            print(current_parms)
            return current_parms
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
    now = datetime.now()
    print("Running at {}hs PST on the day {} (d/m/Y).".format(now.strftime("%H:%M:%S"),now.strftime("%d/%m/%Y")))

    gradient_descent(learning_rate_params = {"a":3, "b":3}, progress_filename = "calibration_progress_6oct", lower_bound_params={"s_0":1, "T":0.1})



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
    
  

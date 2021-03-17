import os
import subprocess
import random
import filecmp
import re
import pytest
import pandas as pd
from termcolor import colored
import csv
import numpy as np
from tqdm import tqdm
from pdb import set_trace as st  # used for debugging


# ========================== Aux ==========================
def write_options_file(params):
    filedata = "\n".join(["[General]",\
                        "GUI=false",\
                        "USE_CPU=false",\
                        "NETWORK_PATH=berkeley_2018/new_full_network/",\
                        "USE_JOHNSON_ROUTING=false",\
                        "USE_SP_ROUTING=true",\
                        "USE_PREV_PATHS=false",\
                        "LIMIT_NUM_PEOPLE=256000",\
                        "ADD_RANDOM_PEOPLE=false",\
                        "NUM_PASSES=1",\
                        "TIME_STEP=0.5",\
                        "START_HR=5",\
                        "END_HR=12",\
                        "SHOW_BENCHMARKS=false",\
                        "OD_DEMAND_FILENAME=od_demand_5to12.csv",\
                        "REROUTE_INCREMENT=15",\
                        " "])

    for (parameter_name, parameter_value) in params.items():
        filedata = re.sub('{}=(-|(0-9)|.|_)*\n'.format(parameter_name), "{}={}\n".format(parameter_name, parameter_value), filedata)

    with open('command_line_options.ini', 'w') as file:
        file.write(filedata)


def log(text, color = 'cyan'):
    print(colored(text, color))


def length_of_csv(csv_file, delimiter):
    with open(csv_file) as file_object:
        return sum(1 for row in file_object)

def route_csv_string_to_list(route_csv_string):
    route_csv_string = route_csv_string.replace("[", "").replace("]", "")
    route_list = route_csv_string.split(",")
    route_list = route_list[:-1]  # delete last extra comma
    return route_list
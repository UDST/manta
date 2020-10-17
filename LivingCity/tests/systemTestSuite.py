""" System tests that consider LivingCity as a black box,
    where inputs are given and outputs are validated to be correct

    REQUIREMENTS: pytest (Install with pip3 install pytest)
"""

import subprocess
import filecmp
import re
import pytest
import pandas as pd
from pdb import set_trace as st  # used for debugging


def write_options_file(params):
    filedata = """[General]
                GUI=false
                USE_CPU=false
                NETWORK_PATH=berkeley_2018/new_full_network/
                USE_JOHNSON_ROUTING=false
                USE_SP_ROUTING=true
                USE_PREV_PATHS=false
                LIMIT_NUM_PEOPLE=256000
                ADD_RANDOM_PEOPLE=false
                NUM_PASSES=1
                TIME_STEP=0.5
                START_HR=5
                END_HR=12
                """

    for (parameter_name, parameter_value) in params.items():
        filedata = re.sub('{}=(-|(0-9)|.)*\n'.format(parameter_name),
                          "{}={}\n".format(parameter_name, parameter_value),
                          filedata)

    with open('command_line_options.ini', 'w') as file:
        file.write(filedata)


""" 
    * Runs LivingCity with prev_paths = false and then with prev_paths = true
    * Produces 0_people5to12_first_run.csv, 0_indexPathVec5to12_first_run.csv, 0_route5to12_first_run.csv on the first run
    * Produces 0_people5to12_second_run.csv, 0_indexPathVec5to12_second_run.csv, 0_route5to12_second_run.csv on the second run
    * Should be no difference between route and indexPathVec files
    * Should be no difference in certain columns of the people files:
        p,init_intersection,end_intersection,time_departure,distance,a,b,T
"""
# ========================== Aux ==========================
def aux_prev_paths_preserve_output_files_in_given_network(network_path):
    output_files = ["route", "people", "indexPathVec"]

    # Previous files are deleted
    for f in output_files:
        subprocess.call(
            "rm ./0_{0}5to12.csv ./0_people5to12_first_run.csv ./0_{0}5to12_{0}_run.csv".format(f), shell=True)

    # First run
    write_options_file({"NETWORK_PATH":"{}".format(network_path), "USE_PREV_PATHS":"false"})
    subprocess.run(["./LivingCity", "&"])
    for f in output_files:
        subprocess.call(
            "mv ./0_{0}5to12.csv ./0_{0}5to12_first_run.csv".format(f), shell=True)

    # Second run
    write_options_file({"NETWORK_PATH":"{}".format(network_path), "USE_PREV_PATHS":"true"})
    subprocess.run(["./LivingCity", "&"])
    for f in output_files:
        subprocess.call(
            "mv ./0_{0}5to12.csv ./0_{0}5to12_second_run.csv".format(f), shell=True)

    # Compare both runs, routes and indexPathVec should be equal
    for f in ["route", "indexPathVec"]:
        assert filecmp.cmp("0_{0}5to12_first_run.csv".format(f), "0_{0}5to12_second_run.csv".format(
            f), shallow=False), "{} file are not equal between runs".format(f)

    # Compare both people files, they should be equal in certain columns:
    for (row_first_run, row_second_run) in pd.read_csv("0_people5to12_first_run.csv", "0_people5to12_second_run.csv"):
        parameters_that_should_be_equal = [
            "p", "init_intersection", "end_intersection", "time_departure", "distance", "a", "b", "T"]
        for param in parameters_that_should_be_equal:
            assert row_first_run[param] == row_second_run, "{} parameter is not equal in the people file between runs".format(
                param)


# ========================== Tests ==========================

@pytest.mark.skip()
def test01_prev_paths_preserve_output_files_in_small_network_1():
    aux_prev_paths_preserve_output_files_in_given_network("test/small_network_1")

@pytest.mark.skip()
def test02_prev_paths_preserve_output_files_in_small_network_2():
    aux_prev_paths_preserve_output_files_in_given_network("test/small_network_2")

def test03_prev_paths_preserve_output_files_in_new_full_network():
    aux_prev_paths_preserve_output_files_in_given_network("berkeley_2018/new_full_network/")

test03_prev_paths_preserve_output_files_in_new_full_network()
import subprocess
import random
import re
import pandas as pd
import numpy as np
import time
from pdb import set_trace as st


def initialize_combination():
    return {'A': random.uniform(1, 10),
            'B': random.uniform(1, 10),
            'T': random.uniform(0.1, 2.0),
            's_0': random.uniform(1, 5.0)}


def write_options_file(a_combination):
    # Read in the file
    with open('command_line_options.ini', 'r') as file:
        filedata = file.read()

    # Replace the parameter values
    for parameter_name in a_combination.keys():
        parameter_value = a_combination[parameter_name]
        filedata = re.sub('\n{}=*\n'.format(parameter_name),
                          '\n{}={}\n'.format(parameter_name, parameter_value), filedata)

    # to do: make sure prev paths are true or false?
    #USE_PREV_PATHS=false

    # Write the file out again
    with open('command_line_options.ini', 'w') as file:
        file.write(filedata)


def perform_step(a_combination, avg_v_mean_current_combination, avg_v_mean_goal, gradient_epsilon, learning_rate):

    gradient = {}

    # we iterate through all parameters
    for parameter in a_combination.keys():
        # we calculate gradient in the direction of the current parameter
        parameter_prev_value = a_combination[parameter]
        a_combination[parameter] += gradient_epsilon
        write_options_file(a_combination)

        subprocess.run(["./LivingCity", "&"],check=True)

        df_people = pd.read_csv("0_people5to12.csv")
        avg_v_mean_new_combination = df_people["avg_v(mph)"].mean()
        a_combination[parameter] = parameter_prev_value

        # (f(x+e) - f(x)) / e
        gradient[parameter] = (abs(avg_v_mean_current_combination - avg_v_mean_goal) -
                               abs(avg_v_mean_new_combination - avg_v_mean_goal))/gradient_epsilon

    # we update the combination according to the gradient and the learning rate
    for parameter in a_combination.keys():
        a_combination[parameter] += gradient[parameter] * learning_rate

    return a_combination


"""
avg_v_goal = avg mph that we want to approximate to
avg_v_goal_margin = max delta (in mph) that we're willing to have in our RMSE
gradient_epsilon = e in (f(x+e)-f(x))/e when calculating the gradient
"""


def gradient_descent(avg_v_mean_goal=30, avg_v_goal_margin=0.02, gradient_epsilon=0.0001, learning_rate=0.01, max_steps=10):
    iteration = 0
    combination = initialize_combination()

    print("Running gradient descent with the following hyperparameters:")
    print("> Average velocity mean goal: {}")
    print("> Goal margin: {}")
    print("> epsilon when calculating gradient: {}")
    print("> learning rate: {}")
    print("> Maximum number of steps: {}")
    print("> Initial combination given at random: {}".format(combination))

    for step in range(max_steps):
        print("\n Step #{}".format(step))
        # load street network and edge/node data
        df_people = pd.read_csv("0_people5to12.csv")
        avg_v_mean_current_combination = df_people["avg_v(mph)"].mean()

        print("> Current combination:")
        print("> {} | average velocity mean: {}".format(
            combination, avg_v_mean_current_combination))

        if (abs(avg_v_mean_goal - avg_v_mean_current_combination) < avg_v_goal_margin):
            print("Current combination reached goal. Stopping.")
            return

        # calculate gradient
        combination = perform_step(
            combination, avg_v_mean_current_combination, avg_v_mean_goal, gradient_epsilon, learning_rate)

    print("Goal not reached at {} steps".format(max_steps))
    print("Final combination:")
    df_people = pd.read_csv("0_people5to12.csv")
    avg_v_mean_current_combination = df_people["avg_v(mph)"].mean()
    print("> {} | average velocity mean: {}".format(
        combination, avg_v_mean_current_combination))


if (__name__ == "__main__"):
    gradient_descent()

""" This file contains benchmarking scripts.
    Requirements under linux:
    - Python 3.6.5
"""
import os
import sys
import time
import json
import subprocess
import re
import statistics
import pandas as pd
from termcolor import colored
import psutil
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import xlwt
import argparse
from pdb import set_trace as st
import numpy as np
from datetime import date
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

# ================= Aux =================

def write_options_file(params = None):
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
                            "START_HR=0",\
                            "END_HR=24",\
                            "SHOW_BENCHMARKS=true",\
                            "OD_DEMAND_FILENAME=activity_od_demand_0to24_new.csv", \
                            " "])

    if params is not None:
        for (parameter_name, parameter_value) in params.items():
            filedata = re.sub('{}=(-|(0-9)|.)*\n'.format(parameter_name),
                            "{}={}\n".format(parameter_name, parameter_value),
                            filedata)

    with open('command_line_options.ini', 'w') as file:
        file.write(filedata)

def log(text):
    print(colored(text, 'cyan'))

def obtain_gpu_memory_used():
    # we execute nvidia-smi asking for the gpu memory used right now
    # it should output it in MiB
    command = ("nvidia-smi --query-gpu=memory.used --format=csv")
    assert subprocess.check_output(command.split()).decode('ascii').split('\n')[:-1][-1].split()[1] == "MiB"
    return int(subprocess.check_output(command.split()).decode('ascii').split('\n')[:-1][-1].split()[0])

def time_since_epoch_miliseconds():
    return int(time.time()*1000)

def convertMiBtoMB(value):
    return value * 1.049

# ================= Components and resources =================

all_component_names  =  ["Load_network", \
                "Load_OD_demand_data", \
                "Routing_CH", \
                "CH_output_nodes_to_edges_conversion", \
                "Convert_routes_into_GPU_data_structure_format", \
                "File_output", \
                "Lane_Map_creation", \
                "Microsimulation_in_GPU"]

all_resource_names = ["cpu_used", "mem_used", "gpu_memory_used", "Elapsed_time_(ms)"]


# ================= Benchmarking =================

"""
    In order to benchmark how many resources the simulation consumes, this function:
    1. Runs the simulation while polling the resources being consumed for each half second
    2. Parses the simulation output, obtaining the timestamp for each component (when it started and ended)
    3. Combine the components' timestamps with the resources' poll
    4. Save these combinations at benchmarks.csv
"""
def benchmark_one_run(number_of_run, benchmark_name, params):
    network_path = params['NETWORK_PATH']
    # ============== Run the simulation while polling the resources
    resources_timestamps_idle = { \
        "cpu_used": psutil.cpu_percent(), \
        "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
        "gpu_memory_used": convertMiBtoMB(obtain_gpu_memory_used()), \
        "Elapsed_time_(ms)": "NaN"
        }

    write_options_file(params)

    args = "./LivingCity &".format(number_of_run).split()
    start_time_simulation = time.time()
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # While the simulation is running we check the resources and write them down per second
    resources_timestamps = []
    while process.poll() is None:
        resources_timestamps.append({ \
            "time_since_epoch_ms": int(time_since_epoch_miliseconds()), \
            "cpu_used": psutil.cpu_percent(), \
            "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
            "gpu_memory_used": convertMiBtoMB(obtain_gpu_memory_used())
        })
        time.sleep(.5)
    
    full_output = str(process.stdout.read())
    elapsed_time_simulation = time.time() - start_time_simulation
    log("Total simulation duration: {} secs")

    assert process.returncode == 0
    
    # ============== Parse the simulation benchmarks =================
    log("Parsing the simulation benchmarks...")
    all_components = { "Load_network":{}, "Load_OD_demand_data":{}, "Routing_CH":{}, "CH_output_nodes_to_edges_conversion":{}, "Convert_routes_into_GPU_data_structure_format":{}, "File_output":{}, "Lane_Map_creation":{}, "Microsimulation_in_GPU":{}}

    log("Showing full output...")
    full_output_printable = full_output.replace("\\n", "\n")
    print(full_output_printable)
    for output_line in full_output[2:].split("\\n"):
        for component_name,component_timestamp in all_components.items():
            if re.search('<{}>'.format(component_name), output_line, flags=0):
                time_since_epoch = int(output_line.split()[-1]) # the last word is the time since epoch
                if re.search('Started', output_line, flags=0):
                    component_timestamp["Started"] = time_since_epoch
                if re.search('Ended', output_line, flags=0):
                    component_timestamp["Ended"] = time_since_epoch

    # ================= Combine and process resource polls and parsed benchmarks =================
    log("Combining and parsing proces resource polls and parsed benchmarks...")
    for component_name, component_resources in all_components.items():
        component_resources["cpu_used"] = []
        component_resources["mem_used"] = []
        component_resources["gpu_memory_used"] = []
        for one_resource_timestamp in resources_timestamps:
            if one_resource_timestamp["time_since_epoch_ms"] > component_resources["Started"] and \
                one_resource_timestamp["time_since_epoch_ms"] < component_resources["Ended"]:

                component_resources["cpu_used"].append(one_resource_timestamp["cpu_used"])
                component_resources["mem_used"].append(one_resource_timestamp["mem_used"])
                component_resources["gpu_memory_used"].append(one_resource_timestamp["gpu_memory_used"])

    for component_name, component_resources in all_components.items():
        component_resources["Elapsed_time_(ms)"] = component_resources["Ended"] - \
                                                    component_resources["Started"]
        assert "Started" in component_resources.keys() and \
                "Ended" in component_resources.keys()

        del component_resources["Ended"]
        del component_resources["Started"]

    # ================= Write it down to a csv =================
    log("(CPU data is in %, memory is in GB, GPU memory is in MB, elapsed time in ms)")
    with open("benchmarking/{}.csv".format(benchmark_name),"a") as benchmarks_file:
        for component_name in all_component_names:
            line_to_write = [str(number_of_run), component_name]
            for one_resource in all_resource_names:
                if one_resource == "Elapsed_time_(ms)":
                    line_to_write.append(str(all_components[component_name][one_resource]))
                else:
                    if all_components[component_name][one_resource] == []:
                        line_to_write.append("0")
                    else:
                        line_to_write.append(str(statistics.mean(all_components[component_name][one_resource])))

            line_to_write_string = ','.join(line_to_write)
            benchmarks_file.write(line_to_write_string + "\n")
        
        # === idle ===
        line_to_write_idle = [str(number_of_run), "idle"]
        line_to_write_idle += [str(resources_timestamps_idle[resource]) for resource in all_resource_names]
        idle_line_to_write_string = ','.join(line_to_write_idle)
        benchmarks_file.write(idle_line_to_write_string + "\n")

    # ================ Metadata ================
    log("Saving metadata...")
    metadata = {'elapsed_time_simulation': elapsed_time_simulation}
    for output_line in full_output[2:].split("\\n"):
        if re.search('Number of blocks: '.format(component_name), output_line, flags=0):
            metadata['number_of_blocks'] = int(output_line.split()[-1])
        if re.search('Number of threads per block: '.format(component_name), output_line, flags=0):
            metadata['number_of_threads_per_block'] = int(output_line.split()[-1])
    
    df_nodes = pd.read_csv(os.path.join(network_path, "nodes.csv"))
    metadata['number_of_nodes'] = len(df_nodes)

    df_edges = pd.read_csv(os.path.join(network_path, "edges.csv"))
    metadata['number_of_edges'] = len(df_edges)

    with open("benchmarking/metadata_{}.json".format(benchmark_name),"w") as metadata_file:
        json.dump(metadata, metadata_file)

    log("Done")


def benchmark_multiple_runs_and_save_csv(number_of_runs, benchmark_name = "benchmarks", params = None):
    log("Running system benchmarks. Benchmark name: {}. Number of runs: {}. Custom parameters: {}".format(benchmark_name, number_of_runs, params))
    log("Please do not run anything else on this PC until it is finished, since that may alter the benchmarking results.")

    with open("benchmarking/{}.csv".format(benchmark_name),"w") as benchmarks_file:
        head = ["number_of_run,Component", "cpu_used", "mem_used", "gpu_memory_used", "Elapsed_time_(ms)"]
        benchmarks_file.write(",".join(head)+"\n")

    for i in range(number_of_runs):
        log("Running benchmark for run {}...".format(i))
        benchmark_one_run(i, benchmark_name, params)


def load_csv_and_generate_benchmark_report(benchmark_name = "benchmarks", all_in_one_plot = True):
    # ======================= Read metadata =======================
    with open('benchmarking/metadata_{}.json'.format(benchmark_name)) as metadata_file:
        metadata = json.load(metadata_file)

    df_people = pd.read_csv('0_people5to12.csv')

    # ======================= Plot benchmarks =======================
    log("Plotting benchmarks")
    readable_components = { \
            "Load_network":"Load network", \
            "Load_OD_demand_data":"Load OD demand", \
            "Routing_CH":"Routing CH", \
            "CH_output_nodes_to_edges_conversion":"CH output conv.", \
            "Convert_routes_into_GPU_data_structure_format":"Routes conv. to GPU", \
            "File_output":"File output", \
            "Lane_Map_creation":"Lane Map creation", \
            "Microsimulation_in_GPU":"Microsimulation in GPU"
            }

    sns.set_theme(style="whitegrid")

    df = pd.read_csv('benchmarking/{}.csv'.format(benchmark_name))

    if all_in_one_plot:
        fig_benchmarks = plt.figure(figsize=(18,9))

    # change from ms to secs
    df["Elapsed_time_(ms)"] = df["Elapsed_time_(ms)"].apply(lambda x: x / 1000)

    resources = [("RAM usage (GB)", "mem_used", "GB"), ("CPU usage (%)", "cpu_used", "%"), \
        ("GPU memory usage (MB)", "gpu_memory_used", "MB"), ("Elapsed time (secs)", "Elapsed_time_(ms)", "secs")]
    for i, (title, resource_name, resource_axis) in enumerate(resources):
        plt.subplots_adjust(wspace=0.8, hspace=0.8)
        if all_in_one_plot:
            ax_benchmarks = fig_benchmarks.add_subplot(2, 2, i+1)
        else:
            fig_benchmarks = plt.figure(figsize=(16,8))
            ax_benchmarks = matplotlib.axes.Axes(fig = fig_benchmarks, rect = [0,0, 18, 9])
        sns.set(style="whitegrid")

        components_in_order = ["Load_network", \
                            "Load_OD_demand_data", \
                            "Routing_CH", \
                            "CH_output_nodes_to_edges_conversion", \
                            "Convert_routes_into_GPU_data_structure_format", \
                            "Lane_Map_creation", \
                            "Microsimulation_in_GPU", \
                            "File_output"]

        ticks = [readable_components[comp] for comp in components_in_order]
        if resource_name == "Elapsed_time_(ms)":
            df = df.drop(df[df["Component"] == "idle"].index)
            current_components = components_in_order
        else:
            ticks += ["idle"]
            current_components = components_in_order + ["idle"]
        
        ax_benchmarks = sns.barplot(
            data=df,
            x=resource_name, y="Component",
            order=components_in_order,
            capsize=.2
        )
        ax_benchmarks.set_yticks(range(len(ticks)))
        ax_benchmarks.set_yticklabels(ticks)
        plt.ylabel("")
        plt.xlabel(resource_axis)
        if all_in_one_plot:
            plt.title(title)
        else:
            plt.title("MANTA benchmarking  |  {}".format(title))
        component_label_values = [df[df["Component"] == component_name][resource_name].mean() for component_name in current_components]


        def add_custom_legend(ax, txt, fontsize = 12, loc = 1, *args, **kwargs):
            at = AnchoredText(txt,
                            prop=dict(size=fontsize), 
                            frameon=True,
                            loc=loc)
            at.patch.set_boxstyle("round,pad=0.2,rounding_size=0.2")
            ax.add_artist(at)
            return at

        for i, value in enumerate(component_label_values):
            plt.text(value + 0.05*max(component_label_values), i, str(round(value, 3)))
    
        if not all_in_one_plot:
            add_custom_legend(ax_benchmarks, "# Nodes: {}\n# Edges: {}\n# Trips: {}".format(
                metadata["number_of_nodes"], metadata["number_of_edges"], len(df_people)))
            if resource_name == "Elapsed_time_(ms)":
                plt.text(0.73, 0.02, 'Total: {} secs ({} mins)'.format(
                        round(sum(component_label_values),2), round(sum(component_label_values)/60,2)),
                        fontsize=14, transform=plt.gcf().transFigure, weight = 'bold')
            output_filepath = "benchmarking/{}_{}.png".format(benchmark_name, resource_name)
            log("Saving figure at {}".format(output_filepath))
            fig_benchmarks.savefig(output_filepath, dpi = 400)
            plt.cla()
            plt.clf()
            plt.close()
        elif resource_name == "Elapsed_time_(ms)":
            plt.text(0.80, 0.055, 'Total: {} secs ({} mins)'.format(
                    round(sum(component_label_values),2), round(sum(component_label_values)/60,2)),
                    fontsize=12, transform=plt.gcf().transFigure, weight = 'bold')

    number_of_runs = df['number_of_run'].nunique()
    if all_in_one_plot:
        if number_of_runs == 1:
            fig_benchmarks.suptitle("MANTA benchmark for 1 run".format(number_of_runs), fontsize=16)    
        else:
            fig_benchmarks.suptitle("MANTA benchmark for {} runs".format(number_of_runs), fontsize=16)
        
        output_filepath = "benchmarking/{}.png".format(benchmark_name)
        log("Saving figure at {}".format(output_filepath))
        fig_benchmarks.savefig(output_filepath, dpi = 400)
        plt.close()
        plt.clf()

    # ======================= Plot metadata (run time in minutes) =======================
    log("Plotting num steps distribution in minutes...")
    df_num_steps_mins = df_people['num_steps'] / 60
    
    fig_metadata = plt.figure(2, figsize=(6,2.4))
    ax_metadata = df_num_steps_mins.plot.hist(bins=50, alpha=1)
    fig_metadata = ax_metadata.get_figure()
    run_time_minutes_output_filepath = 'benchmarking/run_time_in_minutes_{}.png'.format(benchmark_name)
    log("Saving figure at {}".format(run_time_minutes_output_filepath))
    plt.title("Trip run time distribution (in simulated mins)")
    plt.xlabel("Minutes")
    plt.axvline(df_num_steps_mins.mean(), color='red', linestyle='dashed', linewidth=2)
    ax_metadata.legend(['Mean: {}\nStandard dev: {}'.format(round(df_num_steps_mins.mean(),2), round(df_num_steps_mins.std(),2))])    
    fig_metadata.savefig(run_time_minutes_output_filepath)
    plt.close()
    plt.clf()


    # ======================= Save xls file =======================
    log("Saving xls file at {}.xls".format(benchmark_name))
    
    book = xlwt.Workbook(encoding="utf-8")
    sheet_benchmarking = book.add_sheet("Benchmarking")
    
    title_font = xlwt.Font()
    title_font.bold = True
    title_style = xlwt.XFStyle()
    title_style.title_font = title_font
    sheet_benchmarking.write(0, 0, "Component profiling", style = title_style)
    vertical_offset = 1
    for component_index, component_name in enumerate(components_in_order):
        sheet_benchmarking.write(vertical_offset + component_index+1, 0, component_name)
    
    for resource_index, (resource_name, _, _) in enumerate(resources):
        sheet_benchmarking.write(vertical_offset, resource_index+1, resource_name)

    for component_index, component_name in enumerate(components_in_order):
        for resource_index, (_, resource_name, _) in enumerate(resources):
            sheet_benchmarking.write(vertical_offset + component_index + 1,
                resource_index+1, round(df[df["Component"] == component_name][resource_name].mean(), 2))

    vertical_offset += len(components_in_order) + 1
    
    metadata_to_write = [("Total simulation time (in secs)", round(df['Elapsed_time_(ms)'].sum(), 2)), \
                        ("Total simulation time (in mins)", round(df['Elapsed_time_(ms)'].sum() / 60, 2)), \
                        ("Number of trips simulated", len(df_num_steps_mins)), \
                        ("Number of blocks", round(metadata["number_of_blocks"], 2)), \
                        ("Number of threads per block", round(metadata["number_of_threads_per_block"], 2))]

    trip_duration_simulated_mins = [("Trip duration mean (in simulated mins)", round(df_num_steps_mins.mean(), 2)), \
                                    ("Trip duration std dev (in simulated mins)", round(df_num_steps_mins.std(), 2)), \
                                    ("Total trip duration (in simulated mins)", round(df_num_steps_mins.sum(), 2))]
    
    def write_block_to_xls(block, vertical_position, horizontal_position):
        for (i, (elem_name, elem_value)) in enumerate(block):
            sheet_benchmarking.write(vertical_position + i, horizontal_position, elem_name)
            sheet_benchmarking.write(vertical_position + i, horizontal_position + 1, elem_value)

    vertical_offset += 1
    sheet_benchmarking.write(vertical_offset, 0, "Simulation metadata", style = title_style)

    vertical_offset += 1
    write_block_to_xls(metadata_to_write, vertical_offset, 0)

    vertical_offset += len(metadata_to_write) + 1
    sheet_benchmarking.write(vertical_offset, 0, "Trips metadata", style = title_style)

    vertical_offset += 1
    write_block_to_xls(trip_duration_simulated_mins, vertical_offset, 0)

    sheet_benchmarking.write(0, 7,"Date (Y/M/D)")
    sheet_benchmarking.write(0, 8, date.today().strftime("%Y/%m/%d"))
    
    book.save(os.path.join("benchmarking", "{}.xls".format(benchmark_name)))


def parse_input_arguments():
    argument_values = {}

    parser = argparse.ArgumentParser(description='Runs multiple benchmarks over simulation and outputs a csv, a plot'\
                                                ' and a xls file with the results.')
    parser.add_argument('-i', "--runs", type=int, help='Number of runs to run (min: 1).')
    parser.add_argument('-n', "--name", type=str, help='Name of the benchmark for the output files')

    args = parser.parse_args()

    if not args.runs:
        log("Setting number of runs by default to: 1")
        argument_values['number_of_runs'] = 1
    else:
        argument_values['number_of_runs'] = args.runs
    
    if not args.name:
        log('Setting benchmark name by default to: "benchmark"')
        argument_values['benchmark_name'] = "benchmark"
    else:
        argument_values['benchmark_name'] = args.name
    
    return argument_values

def experiment_increasing_demand_trips():
    runs = 1
    # create all csvs from the 0to24 one
    full_trips = pd.read_csv("../berkeley_2018/new_full_network/activity_od_demand_0to24_new.csv")
    for end in [3, 9, 15, 21]:
        print("running for {} ...".format(end))
        full_trips[full_trips["dep_time"] < end * 3600].to_csv("activity_od_demand_0to{}_new.csv".format(end))

    # run the simulations
    for end in range(6, 24, 3):
        od_demand = "activity_od_demand_0to{}_new.csv".format(end)
        name = "benchmarking_mainbranch_3090_0to{}".format(end)
        benchmark_multiple_runs_and_save_csv(runs, name, {"NETWORK_PATH":"berkeley_2018/new_full_network/", "OD_DEMAND_FILENAME": od_demand, "START_HR":"0", "END_HR": str(end)})
        load_csv_and_generate_benchmark_report(name, all_in_one_plot = True)


if __name__ == "__main__":
    experiment_increasing_demand_trips()
    sys.exit()
    argument_values = parse_input_arguments()
    benchmark_multiple_runs_and_save_csv(argument_values['number_of_runs'],
                                        benchmark_name = argument_values['benchmark_name'])
    load_csv_and_generate_benchmark_report(benchmark_name = argument_values['benchmark_name'],
                                            all_in_one_plot = True)

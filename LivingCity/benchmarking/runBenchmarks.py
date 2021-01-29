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
                            "START_HR=5",\
                            "END_HR=12",\
                            "SHOW_BENCHMARKS=true",\
                            "OD_DEMAND_FILENAME=od_demand_5to12.csv", \
                            "REROUTE_INCREMENT=30"
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

def all_values_are_equal(a_dictionary):
    return all(one_value == list(a_dictionary.values())[0] for one_value in a_dictionary.values())

def all_elems_in_list_are_equal(a_list):
    return all(elem == a_list[0] for elem in a_list)

def count_number_of_reroute_increments_in_full_output(full_output_lines):
    repeated_components_count = {"Routing_CH": 0,
        "CH_output_nodes_to_edges_conversion": 0,
        "Convert_routes_into_GPU_data_structure_format": 0}

    for output_line in full_output_lines:
        for component_name in repeated_components_count.keys():
            if re.search('<{}'.format(component_name), output_line, flags=0):
                if re.search('Started', output_line, flags=0):
                    repeated_components_count[component_name] += 1

    assert all_values_are_equal(repeated_components_count), "Discrepancy in number of repetitions between components"
    
    return list(repeated_components_count.values())[0]

def merge_dictionaries(a, b):
    return {**a, **b}

# ================= Components and resources =================

all_resource_names = ["cpu_used", "mem_used", "gpu_memory_used", "Elapsed_time_(ms)"]


# ================= Benchmarking =================

"""
    In order to benchmark how many resources the simulation consumes, this function:
    1. Runs the simulation while polling the resources being consumed for each half second
    2. Parses the simulation output, obtaining the timestamp for each component (when it started and ended)
    3. Combine the components' timestamps with the resources' poll
    4. Save these combinations at benchmarks.csv
"""
def benchmark_one_run(number_of_run, benchmark_name, network_path):
    # ============== Run the simulation while polling the resources
    resources_timestamps_idle = { \
        "cpu_used": psutil.cpu_percent(), \
        "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
        "gpu_memory_used": convertMiBtoMB(obtain_gpu_memory_used()), \
        "Elapsed_time_(ms)": "NaN"
        }

    write_options_file({"NETWORK_PATH": network_path})

    with open('LivingCity.out', 'w') as outputFile:
        args = "./LivingCity &".split()
        start_time_simulation = time.time()
        process = subprocess.Popen(args, stdout=outputFile, stderr=outputFile, shell=True)
        
        # While the simulation is running we check the resources every half second
        full_output = ""
        resources_timestamps = []
        while process.poll() is None:
            resources_timestamps.append({ \
                "time_since_epoch_ms": int(time_since_epoch_miliseconds()), \
                "cpu_used": psutil.cpu_percent(), \
                "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
                "gpu_memory_used": convertMiBtoMB(obtain_gpu_memory_used())
            })
            time.sleep(.5)

    elapsed_time_simulation = time.time() - start_time_simulation
    log("Total simulation duration: {} secs")
    assert process.returncode == 0

    with open('LivingCity.out', 'r') as outputFile:
        full_output = outputFile.read()
        print(full_output)

    # ============== Parse the simulation benchmarks =================
    log("Parsing the simulation benchmarks...")
    
    full_output_lines = full_output.split("\n")
    number_of_reroute_increments = count_number_of_reroute_increments_in_full_output(full_output_lines)

    all_components_keys = ["{}{}".format(component, i)
        for i in range(number_of_reroute_increments)
        for component in ["Routing_CH_batch_",
                          "CH_output_nodes_to_edges_conversion_batch_",
                          "Convert_routes_into_GPU_data_structure_format_batch_",
                          "Microsimulation_in_GPU_batch_"]
        ] + ["Load_network",
            "Load_OD_demand_data",
            "File_output",
            "Lane_Map_creation"]

    all_components = {one_component_key:{} for one_component_key in all_components_keys}
    
    for output_line in full_output_lines:
        for component_name, component_timestamp in all_components.items():
            if re.search('<{}>'.format(component_name), output_line, flags=0):
                time_since_epoch = int(output_line.split()[-1]) # the last word is the time since epoch
                try:
                    if re.search('Started', output_line, flags=0):
                        component_timestamp["Started"] = time_since_epoch
                    if re.search('Ended', output_line, flags=0):
                        component_timestamp["Ended"] = time_since_epoch
                except:
                    st()

    # ================= Combine resource polling and parsed benchmarks =================
    log("Combining process resource polling and parsed benchmarks...")
    for component_name, component_resources in all_components.items():
        component_resources["cpu_used"] = []
        component_resources["mem_used"] = []
        component_resources["gpu_memory_used"] = []
        for one_resource_timestamp in resources_timestamps:
            try:
                if one_resource_timestamp["time_since_epoch_ms"] > component_resources["Started"] and \
                    one_resource_timestamp["time_since_epoch_ms"] < component_resources["Ended"]:

                    component_resources["cpu_used"].append(one_resource_timestamp["cpu_used"])
                    component_resources["mem_used"].append(one_resource_timestamp["mem_used"])
                    component_resources["gpu_memory_used"].append(one_resource_timestamp["gpu_memory_used"])
            except:
                st()

    try:
        for component_name, component_resources in all_components.items():
            component_resources["Elapsed_time_(ms)"] = component_resources["Ended"] - \
                                                        component_resources["Started"]
            assert "Started" in component_resources.keys() and \
                    "Ended" in component_resources.keys()
    except:
        st()
    st()


    # ================= Write it down to a csv =================
    log("(CPU data is in %, memory is in GB, GPU memory is in MB, elapsed time in ms)")
    with open("benchmarking/{}.csv".format(benchmark_name),"a") as benchmarks_file:
        for component_name in all_components:
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


def benchmark_multiple_runs_and_save_csv(number_of_runs,
                                        benchmark_name,
                                        network_path):
    log("Running system benchmarks for network {}. Benchmark name: {}. Number of runs: {}".format(network_path, benchmark_name, number_of_runs))
    log("Please do not run anything else on this PC until it is finished, since that may alter the benchmarking results.")

    with open("benchmarking/{}.csv".format(benchmark_name),"w") as benchmarks_file:
        head = ["number_of_run,Component", "cpu_used", "mem_used", "gpu_memory_used", "Elapsed_time_(ms)"]
        benchmarks_file.write(",".join(head)+"\n")

    for i in range(number_of_runs):
        log("Running benchmark for run {}...".format(i))
        benchmark_one_run(i, benchmark_name, network_path)


def load_csv_and_generate_benchmark_report(benchmark_name = "benchmarks", all_in_one_plot = True):
    # ======================= Read metadata =======================
    with open('benchmarking/metadata_{}.json'.format(benchmark_name)) as metadata_file:
        metadata = json.load(metadata_file)

    df_people = pd.read_csv('0_people5to12.csv')

    # ======================= Plot benchmarks =======================
    log("Plotting benchmarks")

    df = pd.read_csv('benchmarking/{}.csv'.format(benchmark_name))
    repeated_components =   ["Routing_CH",
                            "CH_output_nodes_to_edges_conversion",
                            "Convert_routes_into_GPU_data_structure_format"]

    def get_rows_with_component_prefix(df, prefix):
        return df[df['Component'].str[:len(prefix)].str.contains(prefix)]
    

    repeated_components_count = [len(get_rows_with_component_prefix(df, prefix))
                            for prefix in repeated_components]

    
    assert all_elems_in_list_are_equal(repeated_components_count), \
        "Discrepancy between repeated components in benchmark file."

    for component in repeated_components:
        for batch_number in range(repeated_components_count[0]):
            replace_dictionary = {'Component':'{}_batch_{}'.format(component, batch_number)}
            df = df.replace(replace_dictionary, component)

    df = df.groupby(['number_of_run', 'Component']).agg({'cpu_used': ['mean'],'mem_used': ['mean'],'gpu_memory_used': ['mean'],'Elapsed_time_(ms)': ['mean']})

    df = df.reset_index()
    df.columns = ['number_of_run', 'Component', 'cpu_used', 'mem_used', 'gpu_memory_used', 'Elapsed_time_(ms)']
    

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


    if all_in_one_plot:
        fig_benchmarks = plt.figure(figsize=(18,9))

    # change from ms to secs
    df["Elapsed_time_(ms)"] = df["Elapsed_time_(ms)"].apply(lambda x: x / 1000)
    resources = [("RAM usage (GB)", "mem_used", "GB"), ("CPU usage (%)", "cpu_used", "%"), \
        ("GPU memory usage (MB)", "gpu_memory_used", "MB"), ("Elapsed time (secs)", "Elapsed_time_(ms)", "secs")]
    for i, (title, resource_name, resource_axis) in enumerate(resources):
        df_current_resource = df.copy()
        st()
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
        df_current_resource = df_current_resource.drop(df_current_resource[df_current_resource[resource_name] == 0].index)
        if resource_name == "Elapsed_time_(ms)":
            df_current_resource = df_current_resource.drop(df_current_resource[df_current_resource["Component"] == "idle"].index)
            current_components = components_in_order
        else:
            ticks += ["idle"]
            current_components = components_in_order + ["idle"]

        ax_benchmarks = sns.barplot(
            data=df_current_resource,
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
        component_label_values = [df_current_resource[df_current_resource["Component"] == component_name][resource_name].mean() for component_name in current_components]


        def add_custom_legend(ax, txt, fontsize = 12, loc = 1, *args, **kwargs):
            at = AnchoredText(txt,
                            prop=dict(size=fontsize), 
                            frameon=True,
                            loc=loc)
            at.patch.set_boxstyle("round,pad=0.2,rounding_size=0.2")
            ax.add_artist(at)
            return at

        for i, value in enumerate(component_label_values):
            if resource_name == "Elapsed_time_(ms)":
                plt.text(value + 0.05*max(component_label_values), i, format_seconds_to_mins_and_seconds(value))
            else:
                plt.text(value + 0.05*max(component_label_values), i, str(round(value, 3)))
    
        if not all_in_one_plot:
            add_custom_legend(ax_benchmarks, "# Nodes: {}\n# Edges: {}\n# Trips: {}".format(
                metadata["number_of_nodes"], metadata["number_of_edges"], len(df_people)))
            if resource_name == "Elapsed_time_(ms)":
                plt.text(0.73, 0.02, 'Total: {} secs ({} mins)'.format(
                        round(metadata['elapsed_time_simulation'],2), round(metadata['elapsed_time_simulation']/60,2)),
                        fontsize=14, transform=plt.gcf().transFigure, weight = 'bold')
            output_filepath = "benchmarking/{}_{}.png".format(benchmark_name, resource_name)
            log("Saving figure at {}".format(output_filepath))
            fig_benchmarks.savefig(output_filepath, dpi = 400)
            plt.cla()
            plt.clf()
            plt.close()
        elif resource_name == "Elapsed_time_(ms)":
            plt.text(0.80, 0.055, 'Total: {} secs ({} mins)'.format(
                    round(metadata['elapsed_time_simulation'],2), round(metadata['elapsed_time_simulation']/60,2)),
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


def benchmark_runtime_varying_reroute_increment_and_save_csv(network_path = "berkeley_2018/new_full_network/"):
    for increment in [5, 10] + list(range(15, 181, 15)):
        log("Running the simulation with increment {}".format(increment))
        write_options_file({"NETWORK_PATH": network_path, "REROUTE_INCREMENT": str(increment)})
        start_time_simulation = time.time()
        os.system("./LivingCity")
        #process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elapsed_time_simulation = time.time() - start_time_simulation
        log("Finished. Elapsed time: {}".format(elapsed_time_simulation))

        with open("benchmarking/routingTimes.csv","a") as reroute_runtime_file:
            line_to_write = ["dynamic_routing_{}_mins".format(increment), str(elapsed_time_simulation)]
            line_to_write_string = ','.join(line_to_write)
            reroute_runtime_file.write(line_to_write_string + "\n")

def format_seconds_to_mins_and_seconds(total_seconds):
    total_seconds = int(total_seconds)
    
    seconds = int(total_seconds % 60)
    minutes = int(total_seconds / 60)

    if minutes == 0:
        return "{}s".format(seconds)
    else:
        return "{}m{}s".format(minutes, seconds)

def load_routing_times_csv_and_plot():
    df = pd.read_csv("benchmarking/routingTimes.csv")
    sns.color_palette("pastel", 8)
    fig, ax = matplotlib.pyplot.subplots(figsize=(12,8))

    plt.title("Runtime varying dynamic routing increment")

    ax = sns.barplot(
        data=df,
        x='runtime_in_secs',y="experiment",
        order=['dynamic_routing_{}_mins'.format(mins) for mins in [5, 10] + list(range(15, 181, 15))] + ['static'],
        palette='Blues_d'
    )
    plt.xlabel("")
    plt.ylabel("")

    for bar in ax.patches: 
        ax.annotate(format_seconds_to_mins_and_seconds(bar.get_width()),
                   (bar.get_width() + 50, bar.get_y() + bar.get_height() / 2),
                   ha='center', va='center', 
                   size=12, xytext=(0, 8), 
                   textcoords='offset points') 
  
    plt.subplots_adjust(left=0.2, bottom=0.01, right=0.92, top=0.9)

    fig.savefig("benchmarking/rerouting_times.png")


def parse_input_arguments():
    argument_values = {}

    parser = argparse.ArgumentParser(description='Runs multiple benchmarks over simulation and outputs a csv, a plot'\
                                                ' and a xls file with the results.')
    parser.add_argument('-i', "--runs", type=int, help='Number of runs to run (min: 1).')
    parser.add_argument('-name', "--name", type=str, help='Name of the benchmark for the output files.')
    parser.add_argument('-network', "--networkpath", type=str, help='Relative address of network path folder.')
    parser.add_argument('-p', "--justplot", action="store_true", help='Just plots from the csv without running it previously.')

    args = parser.parse_args()

    if args.runs:
        argument_values['number_of_runs'] = args.runs
    else:
        log("Setting number of runs by default to: 1")
        argument_values['number_of_runs'] = 1
    
    if args.name:
        argument_values['benchmark_name'] = args.name
    else:
        log('Setting benchmark name by default to: "benchmark"')
        argument_values['benchmark_name'] = "benchmark"

    if args.justplot:
        argument_values['just_plot'] = args.justplot
    else:
        argument_values['just_plot'] = False

    if args.networkpath:
        argument_values['network_path'] = args.networkpath
    else:
        argument_values['network_path'] = "berkeley_2018/new_full_network/"
    return argument_values

if __name__ == "__main__":
    argument_values = parse_input_arguments()

    if not argument_values['just_plot']:
        benchmark_multiple_runs_and_save_csv(argument_values['number_of_runs'],
                                            benchmark_name = argument_values['benchmark_name'],
                                            network_path = argument_values['network_path'])

    load_csv_and_generate_benchmark_report(benchmark_name = argument_values['benchmark_name'],
                                            all_in_one_plot = True)

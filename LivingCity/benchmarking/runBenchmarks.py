""" This file contains benchmarking scripts.
    Requirements under linux:
    - Python 3.6.5
"""
import time
import subprocess
import re
import statistics
import pandas as pd
from termcolor import colored
import psutil
import seaborn as sns
import matplotlib.pyplot as plt
from pdb import set_trace as st

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
                            "SHOW_BENCHMARKS=true"])

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
    In order to benchmark how many resources LivingCity consumes, this function:
    1. Runs the simulation while polling the resources being consumed for each half second
    2. Parses the simulation output, obtaining the timestamp for each component (when it started and ended)
    3. Combine the components' timestamps with the resources' poll
    4. Save these combinations at benchmarks.csv
"""
def benchmark_one_run(number_of_run, network_path):
    # ============== Run the simulation while polling the resources
    resources_timestamps_idle = { \
        "cpu_used": psutil.cpu_percent(), \
        "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
        "gpu_memory_used": convertMiBtoMB(obtain_gpu_memory_used()), \
        "Elapsed_time_(ms)": "NaN"
        }

    write_options_file()

    args = "./LivingCity > benchmarking/runsOutput/run_{}.out &".format(number_of_run).split()
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # While LivingCity is running we check the resources and write them down per second
    resources_timestamps = []
    while process.poll() is None:
        resources_timestamps.append({ \
            "time_since_epoch_ms": int(time_since_epoch_miliseconds()), \
            "cpu_used": psutil.cpu_percent(), \
            "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
            "gpu_memory_used": convertMiBtoMB(obtain_gpu_memory_used())
        })
        time.sleep(.5)
    
    with open('resources_mock.log', 'w') as f:
        for item in resources_timestamps:
            f.write("%s," % item)

    assert process.returncode == 0
    
    # ============== Parse the simulation benchmarks =================
    log("Parsing the simulation benchmarks...")
    all_components = { \
            "Load_network":{}, \
            "Load_OD_demand_data":{}, \
            "Routing_CH":{}, \
            "CH_output_nodes_to_edges_conversion":{}, \
            "Convert_routes_into_GPU_data_structure_format":{}, \
            "File_output":{}, \
            "Lane_Map_creation":{}, \
            "Microsimulation_in_GPU":{} 
            }

    full_output = str(process.stdout.read())
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
    with open("benchmarking/benchmarks.csv","a") as benchmarks_file:
        for component_name in all_component_names:
            line_to_write = [str(number_of_run), component_name]
            for one_resource in all_resource_names:
                if one_resource == "Elapsed_time_(ms)":
                    line_to_write.append(str(all_components[component_name][one_resource]))
                else:
                    line_to_write.append(str(statistics.mean(all_components[component_name][one_resource])))

            line_to_write_string = ','.join(line_to_write)
            benchmarks_file.write(line_to_write_string + "\n")
        
        # === idle ===
        line_to_write_idle = [str(number_of_run), "idle"]
        line_to_write_idle += [str(resources_timestamps_idle[resource]) for resource in all_resource_names]
        idle_line_to_write_string = ','.join(line_to_write_idle)
        benchmarks_file.write(idle_line_to_write_string + "\n")

    log("Done")


def benchmark_multiple_runs(number_of_runs, network_path = "berkeley_2018/new_full_network/"):
    log("Running system benchmarks for network {}. Number of runs: {}".format(network_path, number_of_runs))
    log("Please do not run anything else on this PC until it is finished, since that may alter the benchmarking results.")
    subprocess.call("rm -r benchmarking/benchmarks.csv", shell=True)

    with open("benchmarking/benchmarks.csv","w") as benchmarks_file:
        head = ["number_of_run,Component", "cpu_used", "mem_used", "gpu_memory_used", "Elapsed_time_(ms)"]
        benchmarks_file.write(",".join(head)+"\n")

    for i in range(number_of_runs):
        log("Running benchmark for run {}...".format(i))
        benchmark_one_run(i, network_path)


def plot_benchmarks():
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
    number_of_run = 0

    df = pd.read_csv('benchmarking/benchmarks.csv'.format(number_of_run))

    fig = plt.figure(figsize=(19.20,10.80))

    # change from ms to secs
    df["Elapsed_time_(ms)"] = df["Elapsed_time_(ms)"].apply(lambda x: x / 1000)

    resources = [("RAM usage (GB)", "mem_used", "GB"), ("CPU usage (%)", "cpu_used", "%"), \
        ("GPU memory usage (MB)", "gpu_memory_used", "MB"), ("Elapsed time (secs)", "Elapsed_time_(ms)", "secs")]
    for i, (title, resource_name, resource_axis) in enumerate(resources):
        plt.subplots_adjust(wspace=0.8, hspace=0.8)
        ax = fig.add_subplot(2, 2, i+1)
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
            current_components = components_in_order+["idle"]
        
        ax = sns.barplot(
            data=df,
            x=resource_name, y="Component",
            order = components_in_order,
            capsize=.2
        )
        ax.set_yticks(range(len(ticks)))
        ax.set_yticklabels(ticks)
        plt.ylabel("")
        plt.xlabel(resource_axis)
        plt.title(title)
        
        lista = [df[df["Component"] == component_name][resource_name].mean() for component_name in current_components]
    
        for i, value in enumerate(lista):
            plt.text(value + 0.05*max(lista), i, str(round(value, 3)))

    fig.suptitle("LivingCity benchmarks for {} runs".format(df['number_of_run'].nunique()), fontsize=16)
    fig.show()
    output_filepath = "benchmarking/benchmarks.png"
    log("Saving figure at {}".format("output_filepath"))
    fig.savefig(output_filepath,dpi = 400)


if __name__ == "__main__":
    benchmark_multiple_runs(40)
    plot_benchmarks()

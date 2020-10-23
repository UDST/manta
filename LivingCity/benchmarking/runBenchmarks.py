""" This file contains benchmarking scripts.
    Requirements under linux:
    - Python 3.6.5
"""
import time
import subprocess
import random
import filecmp
import re
import pandas as pd
from termcolor import colored
import csv
import psutil
from pdb import set_trace as st


mock_output = """[TIME] Started: <Load_network>. Time since epoch (ms): 1603407603301
berkeley_2018/new_full_network/edges.csv as edges file
berkeley_2018/new_full_network/nodes.csv as nodes file
berkeley_2018/new_full_network/od_demand_5to12.csv as OD file
total edges = 547697
# of edges: 540827
# of vertices: 223328
[TIME] Ended <Load_network>. Elapsed: 8227 ms. Time since epoch (ms): 1603407611529
[TIME] Started: <Load_OD_demand_data>. Time since epoch (ms): 1603407611529
[TIME] Ended <Load_OD_demand_data>. Elapsed: 7400 ms. Time since epoch (ms): 1603407618929
# of OD pairs = 3441953
Loading berkeley_2018/new_full_network/all_paths_ch.txt as paths file
Person 0: init 1723768799 end 65319536 time 35570.496094 a 6.377270 b 3.590473 T 0.888143
loadB18TrafficPeople: People 3441953
[TIME] Started: <Lane_Map_creation>. Time since epoch (ms): 1603407638732
  >> createLaneMap
max_edge_id 547696
Edge count is exactly like number of edges: 540827
Num edges 540827 Num Lanes 610933 Num Lanes Width 627199 Max Leng 26716.257812 Max num lanes 11
Total Memory 1284503552
Start intersection info
[TIME] Ended <Lane_Map_creation>. Elapsed: 8342 ms. Time since epoch (ms): 1603407647074
[TIME] Started: <Microsimulation_in_GPU>. Time since epoch (ms): 1603407647074
traffic person vec size = 3441953
index path vec size = 222515167
edgesData size = 627199
laneMap size = 1284503552
intersections size = 223328
trafficPersonVec size = 3441953
indexPathVec size = 222515167
edgesData size = 627199
laneMap size = 1284503552
intersections size = 223328
<<b18InitCUDA
GPU memory usage: used = 2951, free = 5166 MB, total = 8118 MB
First time_departure 18000.000000
Running main loop from 5 to 12 with 3441953 person...
Total # iterations = 504000%) with 4.67 ms per simulation step (average over 1800)
Total num steps 12203238400.0 Avg 29.55 min Avg CO 727.69
Simulation time = 266444 ms
[TIME] Started: <File output>. Time since epoch (ms): 1603407921590
Save indexPathVec 3441953
Save route 3441953
Save people 3441953
[TIME] Ended <File output>. Elapsed: 67260 ms. Time since epoch (ms): 1603407988850
[TIME] Ended <Microsimulation_in_GPU>. Elapsed: 341808 ms. Time since epoch (ms): 1603407988882

>>Simulation Ended"""
def write_options_file(params):
    filedata = """[General]
                mem_used = GUI=false
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

def log(text):
    print(colored(text, 'cyan'))

def length_of_csv(csv_file, delimiter):
    with open(csv_file) as file_object:
        return sum(1 for row in file_object)


def obtain_gpu_memory_used():
    command = ("nvidia-smi --query-gpu=memory.used --format=csv")
    return int(subprocess.check_output(command.split()).decode('ascii').split('\n')[:-1][-1].split()[0])

def time_since_epoch_miliseconds():
    return int(time.time()*100)

csv_columns  =  ["Load_network", \
                "Load_OD_demand data", \
                "Routing_CH", \
                "CH_output_nodes_to_edges_conversion", \
                "Convert_routes_into_GPU_data_structure_format", \
                "File_output", \
                "Lane_Map_creation", \
                "Microsimulation_in_GPU"]

def benchmark_one_run(number_of_run, network_path = "berkeley_2018/new_full_network/"):
    # ============== Run the simulation while polling the resources
    log("Running system benchmarks for network {}".format(network_path))
    log("Please do not run anything else on this PC until it is finished, since that may alter the benchmarking results.")

    #subprocess.call("./LivingCity > benchmarking/runsOutput/run_{}.out &".format(number_of_run), shell=True)
    args = "./LivingCity > benchmarking/runsOutput/run_{}.out &".format(number_of_run).split()
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    
    resources_timestamps_idle = { \
        "cpu_used": psutil.cpu_percent(), \
        "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
        "gpu_memory_used": obtain_gpu_memory_used()
        }
    
    # While LivingCity is running we check the resources and write them down per second
    resources_timestamps = []
    while process.poll() is None:
        resources_timestamps.append({ \
            "time_since_epoch_ms": int(time_since_epoch_miliseconds()), \
            "cpu_used": psutil.cpu_percent(), \
            "mem_used": int(psutil.virtual_memory()._asdict()["used"]) >> 30, \
            "gpu_memory_used": obtain_gpu_memory_used()
        })
        time.sleep(.5)
    
    with open('resources_mock.log', 'w') as f:
        for item in resources_timestamps:
            f.write("%s," % item)



    assert process.returncode == 0
    
    # ============== Parse the simulation benchmarks
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


    # ================= Combine and process resource polls and parsed benchmarks
    log("Combining and parsing proces resource polls and parsed benchmarks...")
    for component_name, component_resources in all_components.items():
        for one_resource_timestamp in resources_timestamps:
            try:
                if one_resource_timestamp["time_since_epoch_ms"] > component_resources["Started"] and \
                    one_resource_timestamp["time_since_epoch_ms"] < component_resources["Ended"]:

                    component_resources["cpu_used"].append(one_resource_timestamp["cpu_used"])
                    component_resources["mem_used"].append(one_resource_timestamp["mem_used"])
                    component_resources["gpu_memory_used"].append(one_resource_timestamp["gpu_memory_used"])
            except:
                st()

    for component_name, component_resources in all_components.items():
        component_resources["Elapsed_time_(ms)"] = component_resources["Ended"] - \
                                                    component_resources["Started"]
        assert "Started" in component_resources.keys() and \
                "Ended" in component_resources.keys()

        del component_resources["Ended"]
        del component_resources["Started"]

    # ================= Write it down to a csv
    for one_resource in ["cpu_used", "mem_used", "gpu_memory_used", "Elapsed_time(ms)"]:
        log("Saving csv file for resource {}...".format(one_resource))
        with open("benchmarking/{}.csv".format(one_resource),"a") as benchmarks_file:
            line_to_write = []
            for component_name in csv_columns:
                line_to_write.append(all_components[component_name][one_resource])
            benchmarks_file.write(','.join(line_to_write)+"\n")    

    log("Done")

                
    
    


def benchmark_multiple_runs(number_of_runs):
    subprocess.call("rm -r benchmarking/runsOutput/*", shell=True)
    with open("benchmarking/benchmarks.csv","w") as benchmarks_file:
        benchmarks_file.write(','.join(csv_columns)+"\n")

    for i in range(number_of_runs):
        benchmark_one_run(i)



if __name__ == "__main__":
    try:
        benchmark_multiple_runs(2)
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        last_frame = lambda tb=tb: last_frame(tb.tb_next) if tb.tb_next else tb
        frame = last_frame().tb_frame
        ns = dict(frame.f_globals)
        ns.update(frame.f_locals)
        code.interact(local=ns)
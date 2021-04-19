""" This file contains system tests that consider LivingCity as a black box,
    where inputs are given and outputs are validated.
    
    Tests are divided in Shallow and Deep tests. Shallow tests are fast (a couple of minutes at most),
    while deep tests are more exhaustive and take longer. Shallow tests are executed for all simulated
    people while deep tests are executed for a small set of people given by `test_setup.deep_testing_people`.

    Pytest ensures that each test runs automatically for all Setups given by the global variable `pytest.test_setups`.
    In addition, it ensures its simulation runs before each test, if it was not run before.
    This allows the user to run a specific test, i.e.
    `pytest -x -s tests/outputFilesTestSuite.py::test06_last_edge_of_trip_should_go_to_end_intersection`

    In the context of each test, the variable `test_setup` has all input and output variables regarding its setup.
    Therefore, no additional reading from the disk is necessary, since all output files are already loaded in memory.
    Since these variables are shared between tests, they should not be modified in the scope of the tests, only read.
"""

from utils import *
import threading

class TestSetup:
    def __init__(self, name = "test", parameters = {},
                 network_path="berkeley_2018/new_full_network/", distance_epsilon = 1,
                 deep_testing_people = 1000, run_simulation = True,
                 verbose = False, number_of_threads = 20):
        self.name = name
        self.parameters = parameters
        self.network_path = network_path
        self.distance_epsilon = distance_epsilon
        self.deep_testing_people = deep_testing_people
        self.run_simulation = run_simulation
        self.verbose = verbose
        self.has_run = False
        self.df_route = None
        self.df_people = None
        self.df_edges = None
        self.df_indexPathInit = None
        self.start_hr=parameters["START_HR"]
        self.end_hr=parameters["END_HR"]
        self.number_of_threads = number_of_threads

# ============= Setups on which the tests will run =============
setup_static_5to7 = TestSetup("_static_5to7", parameters={"REROUTE_INCREMENT": 0, "PREV_PATHS": "False",
                                           "START_HR": 5, "END_HR": 7, "OD_DEMAND_FILENAME": "od_demand_5to12.csv"},
              run_simulation=True)

setup_dyn_60m_5to7 = TestSetup("_dyn_60m_5to7", parameters={"REROUTE_INCREMENT": 60, "PREV_PATHS": "False",
                                           "START_HR": 5, "END_HR": 7, "OD_DEMAND_FILENAME": "od_demand_5to12.csv"},
              run_simulation=True)

setup_static_5to12 = TestSetup("_static_5to12", parameters={"REROUTE_INCREMENT": 0, "PREV_PATHS": "False",
                                           "START_HR": 5, "END_HR": 12, "OD_DEMAND_FILENAME": "od_demand_5to12.csv"},
              run_simulation=True)

setup_dyn_60m_5to12 = TestSetup("_dyn_60m_5to12", parameters={"REROUTE_INCREMENT": 60, "PREV_PATHS": "False",
                                           "START_HR": 5, "END_HR": 12, "OD_DEMAND_FILENAME": "od_demand_5to12.csv"},
              run_simulation=True)

setup_dyn_15m_5to12 = TestSetup("_dyn_15m_5to12", parameters={"REROUTE_INCREMENT": 15, "PREV_PATHS": "False",
                                           "START_HR": 5, "END_HR": 12, "OD_DEMAND_FILENAME": "od_demand_5to12.csv"},
              run_simulation=True)

# Premade subsets you can test on
all_test_setups = [setup_static_5to7, setup_dyn_60m_5to7, setup_static_5to12, setup_dyn_60m_5to12, setup_dyn_15m_5to12]
dynamic_test_setups = [setup_dyn_60m_5to7, setup_dyn_60m_5to12, setup_dyn_15m_5to12]
static_test_setups = [setup_static_5to7, setup_static_5to12]

# Initialize test_setups as a list with the setups you want to test
pytest.test_setups = [all_test_setups]

# ======================================

@pytest.fixture()
def run_all_test_setups(request):
    for current_test_setup in pytest.test_setups:
        if current_test_setup.has_run:
            continue
        log("Running setup {}...".format(current_test_setup.name))
        current_test_setup.has_run = True

        if current_test_setup.run_simulation:
            log("Running system test setup on output files for network {}".format(current_test_setup.network_path))

            log("Removing previous run files...")
            for output_file_name in ["route", "people", 'indexPathInit']:
                filename = "./0_{}{}to{}_{}.csv".format(output_file_name, current_test_setup.start_hr, \
                                                        current_test_setup.end_hr, current_test_setup.name)
                if os.path.exists(filename):
                    log("Removing old {}...".format(filename))
                    os.remove(filename)

            log("Running simulation...")
            write_options_file(current_test_setup.parameters)
            subprocess.run(["./LivingCity", "&"], check=True)

            for output_file_name in ["route", "people", 'indexPathInit']:
                filename = "./0_{}{}to{}_{}.csv".format(output_file_name, current_test_setup.start_hr, \
                                                        current_test_setup.end_hr, current_test_setup.name)
                os.rename("./0_{}{}to{}.csv".format(output_file_name, current_test_setup.start_hr, \
                                                    current_test_setup.end_hr), filename)
        else:
            log("Skipping the simulation")

        log("Loading output csv files...")
        current_test_setup.df_people = pd.read_csv("0_people{}to{}_{}.csv". \
            format(current_test_setup.start_hr, current_test_setup.end_hr, current_test_setup.name))

        current_test_setup.df_edges = pd.read_csv(current_test_setup.network_path + "edges.csv")

        current_test_setup.df_route = pd.read_csv("0_route{}to{}_{}.csv". \
            format(current_test_setup.start_hr, current_test_setup.end_hr, current_test_setup.name), delimiter=":")

        current_test_setup.df_indexPathInit = pd.read_csv("0_indexPathInit{}to{}_{}.csv". \
            format(current_test_setup.start_hr, current_test_setup.end_hr, current_test_setup.name))

        # ======== deep people to test ========
        if current_test_setup.deep_testing_people == 'all':
            log("Deep testing all people.")
            current_test_setup.deep_testing_people = list(range(len(current_test_setup.df_people)))
        elif type(current_test_setup.deep_testing_people) == int:
            log("Deep testing on {} random people.".\
                format(current_test_setup.deep_testing_people))
            current_test_setup.deep_testing_people = random.sample(range(len(current_test_setup.df_people)),
                                                                current_test_setup.deep_testing_people)
        else:
            assert type(current_test_setup.deep_testing_people) == list, "deep_testing_people should either be 'all', a number or a list."
        deep_testing_people_indexes = pd.DataFrame(current_test_setup.deep_testing_people)
        deep_testing_people_indexes.columns = ['p']
        current_test_setup.df_deep_testing_people = pd.merge(current_test_setup.df_people, deep_testing_people_indexes)
        log("Finished setup {}.".format(current_test_setup.name))
        log("_____________________________")

    return

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_01_all_output_files_should_have_the_same_length(run_all_test_setups, test_setup):
    log("Test 01 (shallow) | Setup {}. All output files should have the same length.".format(test_setup.name))
    assert len(test_setup.df_people) == len(test_setup.df_route)
    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_02_number_of_edges_in_people_file_should_match_for_cpu_and_gpu(run_all_test_setups, test_setup):
    log("Test 02 (shallow) | Setup {}. The number in the people file should match between CPU and GPU".format(test_setup.name))

    mismatch_people = test_setup.df_people[(test_setup.df_people['path_length_cpu'] !=
                                            test_setup.df_people['path_length_gpu']) & (test_setup.df_people['active'] != 1)]

    assert mismatch_people.empty

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_03_all_unfinished_trips_should_have_same_last_time_simulated(run_all_test_setups, test_setup):
    log("Test 03 (shallow) | Setup {}. Unfinished trips should have as last_time_simulated the same value because they were interrupted"\
    .format(test_setup.name))

    unfinished_trips = test_setup.df_people[test_setup.df_people['active'] == 1]

    assert len(unfinished_trips) == 0 or len(unfinished_trips['last_time_simulated'].value_counts()) == 1

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_04a_distance_in_people_file_should_match_the_distance_in_the_route_file(run_all_test_setups, test_setup):
    log("Test 04a (shallow) | Setup {}. The distance in the people file should match the distance".format(test_setup.name) + \
        " outputted in the route file with a margin of {}.".format(test_setup.distance_epsilon))

    joined_distances = test_setup.df_people.merge(test_setup.df_deep_testing_people).join(test_setup.df_route, lsuffix='_people', rsuffix='_route')
    mismatch_people = joined_distances[abs(joined_distances['distance_people'] - joined_distances['distance_people']) > test_setup.distance_epsilon]
    error_message = str(len(mismatch_people)) + " people have a mismatch in their distance."

    assert len(mismatch_people) == 0, error_message
    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_04b_distance_in_people_file_should_match_the_sum_of_the_edges_in_the_route_file(run_all_test_setups, test_setup):
    log("Test 04b (deep) | Setup {}. The distance in the people file should match the sum of the edges".format(test_setup.name) + \
        " in the route file with a margin of {}.".format(test_setup.distance_epsilon))

    mismatch_people = {}
    df_people_to_iterate = pd.merge(test_setup.df_people, test_setup.df_deep_testing_people)

    def test04b_thread(start, end):
        df_people_chunk = df_people_to_iterate[start:end]
        for _, row in df_people_chunk.iterrows():
            person_id = int(row["p"])

            active = int(float(df_people_chunk[df_people_chunk['p'] == int(float(person_id))]['active'].iloc[0]))
            last_time_simulated = int(df_people_chunk[df_people_chunk['p'] == int(float(person_id))]['last_time_simulated'].iloc[0])

            if not int(float(person_id)) in test_setup.deep_testing_people:
                continue

            route = str(test_setup.df_route[test_setup.df_route['p'] == person_id]['route'].iloc[0])
            route = route.replace("[", "").replace("]", "")
            route = route.split(",")

            duplicated_edges = test_setup.df_edges[test_setup.df_edges.duplicated(subset=['osmid_u', 'osmid_v'])]

            distance_traffic_simulator = int(test_setup.df_route[test_setup.df_route['p'] == person_id]['distance'].iloc[0])

            distance_sum_of_edges = 0
            if route == ['']:
                distance_sum_of_edges = 0
            else:
                for edge_id in route:
                    try:
                        distance_sum_of_edges += test_setup.df_edges[test_setup.df_edges['uniqueid'] == int(float(edge_id))]["length"].iloc[0]
                    except:
                        st()


            distance_people_info = df_people_chunk[df_people_chunk['p'] == int(float(person_id))]['distance'].iloc[0]

            if test_setup.verbose:
                print("Person_id: {} distance_sum_of_edges: {}, distance_people_info: {}, active: {}".format(person_id, distance_sum_of_edges, distance_people_info, active))
            
            if (active == 2 and abs(distance_people_info - distance_sum_of_edges) > test_setup.distance_epsilon):
                mismatch_people[person_id] = {"distance_people_info": distance_people_info,
                                            "distance_sum_of_edges": distance_sum_of_edges,
                                            "active": active,
                                            "last_time_simulated": last_time_simulated,
                                            "distance_traffic_simulator": distance_traffic_simulator}

        error_message = str(len(mismatch_people)) + " people have a mismatch in their distance."

        assert len(mismatch_people) == 0, error_message

    thread_list = []
    people_count = 0
    chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads)
    last_chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads + \
                      len(test_setup.df_people) % test_setup.number_of_threads)
    for i in range(test_setup.number_of_threads):
        end = int(people_count + last_chunk_size if (i == test_setup.number_of_threads - 1) else people_count + chunk_size)
        thread = threading.Thread(target=test04b_thread, args=(people_count, end))
        thread_list.append(thread)
        thread_list[i].start()
        people_count += chunk_size

    for thread in thread_list:
        thread.join()

    error_message = str(len(mismatch_people)) + " people have a mismatch in their distance."
    assert len(mismatch_people) == 0, error_message

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_05_first_edge_of_trip_should_come_from_init_intersection(run_all_test_setups, test_setup):
    log("Test 05 (deep) | Setup {}. The first edge of each trip should come from the init intersection.".format(test_setup.name))
    mismatch_people = {}
    
    def test05_thread(start, end):
        df_people_chunk = test_setup.df_deep_testing_people[start:end]
        for _, row in df_people_chunk.iterrows():
            person_id = int(row["p"])

            init_intersection_according_to_people_file = \
                int(df_people_chunk[df_people_chunk['p'] == person_id]['init_intersection'].iloc[0])
            route_list = route_csv_string_to_list(str(test_setup.df_route[test_setup.df_route['p'] == person_id]['route'].iloc[0]))
            if route_list == []:
                continue
            first_edge = int(route_list[0])
            init_intersection_according_to_route_and_edges_csv = \
                int(test_setup.df_edges[test_setup.df_edges['uniqueid'] == first_edge]['osmid_u'])

            if init_intersection_according_to_people_file != init_intersection_according_to_route_and_edges_csv:
                mismatch_people[person_id] = {"init_intersection_according_to_people_file": init_intersection_according_to_people_file,
                                            "init_intersection_according_to_route_and_edges_csv" : init_intersection_according_to_route_and_edges_csv}

    thread_list = []
    people_count = 0
    chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads)
    last_chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads + \
                      len(test_setup.df_people) % test_setup.number_of_threads)
    for i in range(test_setup.number_of_threads):
        end = int(people_count + last_chunk_size if (i == test_setup.number_of_threads - 1) else people_count + chunk_size)
        thread = threading.Thread(target=test05_thread, args=(people_count, end))
        thread_list.append(thread)
        thread_list[i].start()
        people_count += chunk_size

    for thread in thread_list:
        thread.join()

    error_message = str(len(mismatch_people)) + " people have a mismatch in their initial edge and intersection."
    assert len(mismatch_people) == 0, error_message

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_06_last_edge_of_trip_should_go_to_end_intersection(run_all_test_setups, test_setup):
    log("Test 06 (deep) | Setup {}. The last edge of each trip should finish at the end intersection.".format(test_setup.name))
    mismatch_people = {}
        
    def test06_thread(start, end):
        df_people_chunk = test_setup.df_deep_testing_people[start:end]
        for _, row in df_people_chunk.iterrows():
            person_id = int(row["p"])

            end_intersection_according_to_people_file = \
                int(test_setup.df_people[test_setup.df_people['p'] == person_id]['end_intersection'].iloc[0])
            route_list = route_csv_string_to_list(str(test_setup.df_route[test_setup.df_route['p'] == person_id]['route'].iloc[0]))
            if route_list == []:
                continue

            last_edge = int(route_list[-1])
            end_intersection_according_to_route_and_edges_csv = \
                int(test_setup.df_edges[test_setup.df_edges['uniqueid'] == last_edge]['osmid_v'])

            if end_intersection_according_to_people_file != end_intersection_according_to_route_and_edges_csv:
                mismatch_people[person_id] = {"end_intersection_according_to_people_file": end_intersection_according_to_people_file,
                                            "end_intersection_according_to_route_and_edges_csv" : end_intersection_according_to_route_and_edges_csv}


    thread_list = []
    people_count = 0
    chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads)
    last_chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads + \
                      len(test_setup.df_people) % test_setup.number_of_threads)
    for i in range(test_setup.number_of_threads):
        end = int(people_count + last_chunk_size if (i == test_setup.number_of_threads - 1) else people_count + chunk_size)
        thread = threading.Thread(target=test06_thread, args=(people_count, end))
        thread_list.append(thread)
        thread_list[i].start()
        people_count += chunk_size

    for thread in thread_list:
        thread.join()

    error_message = str(len(mismatch_people)) + " people have a mismatch in their end edge and intersection."
    assert len(mismatch_people) == 0, error_message

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_07_each_edge_of_trip_should_connect_to_the_following_edge(run_all_test_setups, test_setup):
    log("Test 07 (deep) | Setup {}. Each edge of each trip should connect to the following edge.".format(test_setup.name))
    mismatch_people = {}
    
    def test07_thread(start, end):
        df_people_chunk = test_setup.df_deep_testing_people[start:end]
        for _, row in df_people_chunk.iterrows():
            person_id = int(row["p"])

            route_list = route_csv_string_to_list(str(test_setup.df_route[test_setup.df_route['p'] == person_id]['route'].iloc[0]))
            if route_list == []:
                continue

            for first_edge_index in range(len(route_list)-1):
                first_edge = int(route_list[first_edge_index])
                first_edge_end_intersection = int(test_setup.df_edges.loc[test_setup.df_edges['uniqueid'] == first_edge]['osmid_v'])

                second_edge = int(route_list[first_edge_index + 1])
                second_edge_init_intersection = int(test_setup.df_edges.loc[test_setup.df_edges['uniqueid'] == second_edge]['osmid_u'])
                
                error_message = "route for person {} does is not a valid path. Route is {} but edge {} and {} are not" + \
                    "connected in the network. Edge {} finishes at intersection {} but edge {} starts at intersection {}" \
                    .format(person_id, route_list, first_edge, second_edge, first_edge,
                    first_edge_end_intersection, second_edge, second_edge_init_intersection)
                    
                if first_edge_end_intersection != second_edge_init_intersection:
                    mismatch_people[person_id] = error_message

    thread_list = []
    people_count = 0
    chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads)
    last_chunk_size = int(len(test_setup.df_people) / test_setup.number_of_threads + \
                      len(test_setup.df_people) % test_setup.number_of_threads)
    for i in range(test_setup.number_of_threads):
        end = int(people_count + last_chunk_size if (i == test_setup.number_of_threads - 1) else people_count + chunk_size)
        thread = threading.Thread(target=test07_thread, args=(people_count, end))
        thread_list.append(thread)
        thread_list[i].start()
        people_count += chunk_size

    for thread in thread_list:
        thread.join()

    assert len(mismatch_people) == 0
    log("Passed")


@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_08_indexPathInit_should_have_all_different_values(run_all_test_setups, test_setup):
    log("Test 08 (shallow) | Setup {}. Testing if indexPathInit has all different values...".format(test_setup.name))

    assert not test_setup.df_indexPathInit['indexPathInit'].duplicated().any(), \
        "All indexPathInit values should be different"

    log("Passed")


@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_09_distances_shouldnt_have_too_many_repeated_values(run_all_test_setups, test_setup):
    log("Test 09 (shallow) | Setup {}. Testing if the travel distances are different...".format(test_setup.name))
    
    distance_values = test_setup.df_people['distance'].value_counts().tolist()[1:10]
    
    assert all([dist < len(test_setup.df_people) / 10 for dist in distance_values]), \
        "Too many repeated values" + test_setup.df_people['distance'].value_counts()

    log("Passed")

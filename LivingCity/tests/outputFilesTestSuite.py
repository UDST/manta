""" This file contains system tests that consider LivingCity as a black box,
    where inputs are given and outputs are validated.
    Requirements under linux:
    - Python 3.6.5
    - pytest 6.1.1
    - pytest-cov 2.10.1
    - pytest-remotedata 0.3.2
"""

from utils import *

# Global variables across tests
pytest.default_network_setup_has_run = False
pytest.network_path = "berkeley_2018/new_full_network/"
pytest.distance_margin_between_route_and_people_file = 100
pytest.edges_path = "berkeley_2018/new_full_network/edges.csv"
pytest.number_of_people_to_test = 1000
pytest.test_all_people = False
pytest.people_to_test = None
pytest.run_simulation = False

def route_csv_string_to_list(route_csv_string):
    route_csv_string = route_csv_string.replace("[", "").replace("]", "")
    route_list = route_csv_string.split(",")
    route_list = route_list[:-1]  # delete last extra comma
    return route_list

"""
default_network_setup runs automatically before the tests that take default_network_setup as a parameter
"""
@pytest.fixture(params=["berkeley_2018/new_full_network/"], ids=["new_full_network"])
def default_network_setup(request):
    log("Running default network setup...")
    if pytest.default_network_setup_has_run:
        return pytest.network_path
    pytest.default_network_setup_has_run = True

    if pytest.run_simulation:
        log("Running system test setup on output files for network {}".format(pytest.network_path))
        output_files = ["route", "people", "indexPathVec"]

        log("Removing previous run files...")
        for output_file_name in output_files:
            for run in ["", "_first_run", "_second_run"]:
                filename = "./0_{}5to12{}.csv".format(output_file_name, run)
                if os.path.exists(filename):
                    log("Removing old {}...".format(filename))
                    os.remove(filename)

        log("Running simulation...")
        write_options_file({"NETWORK_PATH": "{}".format(
            pytest.network_path), "USE_PREV_PATHS": "false", "REROUTE_INCREMENT": 0})
        subprocess.run(["./LivingCity", "&"], check=True)
    else:
        log("Not running the simulation")

    log("Loading output csv files for testing...")
    pytest.pd_people = pd.read_csv("0_people5to12.csv")
    pytest.pd_edges = pd.read_csv(pytest.edges_path)
    pytest.pd_route = pd.read_csv("0_route5to12.csv", delimiter=':')

    if (pytest.test_all_people):
        log("Testing all people...")
        pytest.people_to_test = list(range(1, len(pytest.pd_people)))
    else:
        log("(Since testing for > 2 million people takes several hours, it is tested on {} random people).".format(pytest.number_of_people_to_test))
        pytest.people_to_test = random.sample(range(1, len(pytest.pd_people)), pytest.number_of_people_to_test)

    log("Finished setup.")

    return pytest.network_path

def test01_all_output_files_should_have_the_same_length(default_network_setup):
    lengths = []
    log("Checking that the length of people and route files are the same...")
    for output_file in ["route", "people"]:
        delimiter = ":" if output_file == "route" else ","
        lengths.append(length_of_csv(
            '0_{}5to12.csv'.format(output_file), delimiter))
    try:
        assert all([lengths[0] == elem for elem in lengths])
    except:
        st()
    log("Passed")

def test02_first_edge_of_trip_should_come_from_init_intersection(default_network_setup):
    log("Testing if the first edge of each trip comes from the int intersection")
    for person_id in pytest.people_to_test:
        init_intersection_according_to_people_file = int(pytest.pd_people.iloc[person_id]['init_intersection'])
        route_list = route_csv_string_to_list(str(pytest.pd_route.iloc[person_id][1]))
        if route_list == []:
            log("Skipping empty route for person {}".format(person_id))
            continue
        first_edge = int(route_list[0])
        init_intersection_according_to_route_and_edges_csv = int(pytest.pd_edges.iloc[first_edge]['osmid_u'])

        assert init_intersection_according_to_people_file == init_intersection_according_to_route_and_edges_csv, \
            "Initial intersections do not match! According to the peoples' file, init intersection is " + \
            init_intersection_according_to_people_file + \
            " but according to the routes' and edges files, the init intersection is " + \
            init_intersection_according_to_route_and_edges_csv

    log("Passed")

def test03_last_edge_of_trip_should_go_to_end_intersection(default_network_setup):
    log("Testing if the last edge of each trip goes to the end intersection")
    for person_id in pytest.people_to_test:
        end_intersection_according_to_people_file = int(pytest.pd_people.iloc[person_id]['end_intersection'])
        route_list = route_csv_string_to_list(str(pytest.pd_route.iloc[person_id][1]))
        if route_list == []:
            log("Skipping empty route for person {}".format(person_id))
            continue

        last_edge = int(route_list[-1])
        end_intersection_according_to_route_and_edges_csv = int(pytest.pd_edges.iloc[last_edge]['osmid_v'])

        assert end_intersection_according_to_people_file == end_intersection_according_to_route_and_edges_csv, \
            "End intersections do not match! According to the peoples' file, end intersection is " + \
            end_intersection_according_to_people_file + \
            " but according to the routes' and edges files, the end intersection is " + \
            end_intersection_according_to_route_and_edges_csv
    log("Passed")

def test_04_each_edge_of_trip_should_connect_to_the_following_edge(default_network_setup):
    log("Testing if each edge of each trip connects to the following edge")
    for person_id in pytest.people_to_test:
        route_list = route_csv_string_to_list(str(pytest.pd_route.iloc[person_id][1]))
        if route_list == []:
            log("Skipping empty route for person {}".format(person_id))
            continue

        for first_edge_index in range(len(route_list)-1):
            first_edge = int(route_list[first_edge_index])
            first_edge_end_intersection = int(pytest.pd_edges.loc[pytest.pd_edges['uniqueid'] == first_edge]['osmid_v'])

            second_edge = int(route_list[first_edge_index + 1])
            second_edge_init_intersection = int(pytest.pd_edges.loc[pytest.pd_edges['uniqueid'] == second_edge]['osmid_u'])
            
            error_message = "route for person {} does is not a valid path. Route is {} but edge {} and {} are not" + \
                "connected in the network. Edge {} finishes at intersection {} but edge {} starts at intersection {}" \
                .format(person_id, route_list, first_edge, second_edge, first_edge,
                first_edge_end_intersection, second_edge, second_edge_init_intersection)
            assert first_edge_end_intersection == second_edge_init_intersection, error_message
                
    log("Passed")

def test_05_distance_in_people_file_should_match_the_sum_of_the_edges_in_the_route_file(default_network_setup):
    log("Comparing that the distance in the people file matches the sum of the edges in the route file with a margin of {}...".format(
        pytest.distance_margin_between_route_and_people_file))
    pd_people = pd.read_csv("0_people5to12.csv")
    pd_edges = pd.read_csv(pytest.edges_path)


    for chunk_route in tqdm(pd.read_csv("0_route5to12.csv", sep=":", chunksize=1000),
                            total=len(pytest.pd_people)/1000):
        for _, row in chunk_route.iterrows():
            person_id = str(row["p"])

            if not pytest.test_all_people and not person_id in pytest.people_to_test:
                continue

            route = str(row["route"])
            route = route.replace("[", "").replace("]", "")
            route = route.split(",")
            route = route[:-1]  # delete last extra comma

            distance_sum_of_edges = 0
            for edge_id in route:
                distance_sum_of_edges += pd_edges.loc[int(edge_id)]["length"]

            distance_people_info = pd_people.loc[int(person_id)]['distance']

            error_message = "Discrepancy has been found for person {}. \
                            Distance according to people info: {}. \
                            Distance according to sum of edges: {}.\n Stopping.". \
                            format(person_id, distance_people_info,
                                   distance_sum_of_edges)
            try:
                assert abs(distance_people_info -
                       distance_sum_of_edges) < pytest.distance_margin_between_route_and_people_file, error_message
            except:
                st()

    log("Passed")
""" This file contains system tests that consider LivingCity as a black box,
    where inputs are given and outputs are validated.
    Requirements under linux:
    - Python 3.6.5
    - pytest 6.1.1
    - pytest-cov 2.10.1
    - pytest-remotedata 0.3.2
"""
from utils import *
from utils import route_csv_string_to_list


class TestSetup:
    def __init__(self, name = "test", parameters = {}, edges_path="berkeley_2018/new_full_network/edges.csv",
                 network_path="berkeley_2018/new_full_network/", distance_epsilon = 500,
                 people_to_test = 1000, run_simulation = True, maximum_percentage_distance_mismatch_allowed = 0.1,
                 verbose = False):
        self.name = name
        self.parameters = parameters
        self.edges_path = edges_path
        self.network_path = network_path
        self.distance_epsilon = distance_epsilon
        self.people_to_test = people_to_test
        self.run_simulation = run_simulation
        assert(maximum_percentage_distance_mismatch_allowed > 0 and maximum_percentage_distance_mismatch_allowed < 1)
        self.maximum_percentage_distance_mismatch_allowed = maximum_percentage_distance_mismatch_allowed
        self.verbose = verbose
        self.has_run = False
        self.pd_route = None
        self.pd_people = None
        self.pd_edges = None
        self.pd_indexPathInit = None
        self.start_hr=parameters["START_HR"]
        self.end_hr=parameters["END_HR"]

# ============= All setups on which the tests will run =============

pytest.test_setups = [
    #TestSetup("_dyn_60m", parameters={"REROUTE_INCREMENT": 60, "PREV_PATHS": "False", "START_HR": 5, "END_HR":7, \
    #        "OD_DEMAND_FILENAME": "od_demand_5to7_no_invalid_trips_no_edge_cases.csv"},
    #        run_simulation=False, people_to_test=list(range(100)))#,
#    TestSetup("_dyn_60m_5to7", parameters={"REROUTE_INCREMENT": 60, "PREV_PATHS": "False", "START_HR": 5, "END_HR":7, \
#            "OD_DEMAND_FILENAME": "od_demand_5to7_no_invalid_trips_no_edge_cases.csv"}, run_simulation=True,
#            people_to_test=list(range(1000)))
    TestSetup("_dyn_60m_5to7_with_edge_Cases", parameters={"REROUTE_INCREMENT": 60, "PREV_PATHS": "False", "START_HR": 5, "END_HR":7, \
            "OD_DEMAND_FILENAME": "od_demand_5to7_no_invalid_trips.csv"}, run_simulation=True,
            people_to_test=list(range(100)))
]

# to do: assert all names and all configurations are differnet

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

        log("Loading output csv files for testing...")
        current_test_setup.pd_people = pd.read_csv("0_people{}to{}_{}.csv". \
            format(current_test_setup.start_hr, current_test_setup.end_hr, current_test_setup.name))

        current_test_setup.pd_edges = pd.read_csv(current_test_setup.network_path + "edges.csv")

        current_test_setup.pd_route = pd.read_csv("0_route{}to{}_{}.csv". \
            format(current_test_setup.start_hr, current_test_setup.end_hr, current_test_setup.name), delimiter=':')

        current_test_setup.pd_indexPathInit = pd.read_csv("0_indexPathInit{}to{}_{}.csv". \
            format(current_test_setup.start_hr, current_test_setup.end_hr, current_test_setup.name))

        if current_test_setup.people_to_test == 'all':
            log("Testing all people...")
            current_test_setup.people_to_test = list(range(1, len(current_test_setup.pd_people)))
        elif type(current_test_setup.people_to_test) == int:
            log("(Since testing for > 2 million people takes several hours, it is tested on {} random people).".\
                format(current_test_setup.people_to_test))
            current_test_setup.people_to_test = random.sample(range(1, len(current_test_setup.pd_people)),
                                                                current_test_setup.people_to_test)
        else:
            assert type(current_test_setup.people_to_test) == list, "people_to_test should either be 'all', a number or a list."

    log("Finished all setups")

    return

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test01_all_output_files_should_have_the_same_length(run_all_test_setups, test_setup):
    log("Test 01: All output files should have the same length for setup {}".format(test_setup.name))
    assert len(test_setup.pd_people) == len(test_setup.pd_route)
    log("Passed test 01")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_02_distance_in_people_file_should_match_the_sum_of_the_edges_in_the_route_file(run_all_test_setups, test_setup):
    log("Test 02: The distance in the people file should match the sum of the edges" + \
        " in the route file with a certain margin for setup {}".format(test_setup.name))
    df_people_to_test = pd.DataFrame(test_setup.people_to_test)
    df_people_to_test.columns = ['p']
    mismatch_people = {}
    for _, row in tqdm(pd.merge(test_setup.pd_route, df_people_to_test).iterrows()):
        person_id = str(row["p"])

        active = int(test_setup.pd_people[test_setup.pd_people['p'] == int(person_id)]['active'].iloc[0])

        if not int(person_id) in test_setup.people_to_test:
            continue

        route = str(row["route"])
        route = route.replace("[", "").replace("]", "")
        route = route.split(",")
        route = route[:-1]  # delete last extra comma

        distance_sum_of_edges = 0
        for edge_id in route:
            try:
                distance_sum_of_edges += test_setup.pd_edges[test_setup.pd_edges['uniqueid'] == int(edge_id)]["length"].iloc[0]
            except:
                st()

        distance_people_info = test_setup.pd_people[test_setup.pd_people['p'] == int(person_id)]['distance'].iloc[0]

        if test_setup.verbose:
            print("Person_id: {} distance_sum_of_edges: {}, distance_people_info: {}, active: {}".format(person_id, distance_sum_of_edges, distance_people_info, active))
        
        if (active == 2 and abs(distance_people_info - distance_sum_of_edges) > test_setup.distance_epsilon):
            mismatch_people[person_id] = {"distance_people_info": distance_people_info,
                                        "distance_sum_of_edges": distance_sum_of_edges,
                                        "active": active}

    log("{} people have mismatched distances. The allowed number in this test is {}."\
        .format(len(mismatch_people), test_setup.maximum_percentage_distance_mismatch_allowed * len(test_setup.pd_people)))

    if (len(mismatch_people.keys()) > len(test_setup.pd_people) * test_setup.maximum_percentage_distance_mismatch_allowed):
        print("Too many people have mismatched distances")
        print(mismatch_people)
        st()

    log("Passed test 02")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test03_first_edge_of_trip_should_come_from_init_intersection(run_all_test_setups, test_setup):
    log("Test 03: The first edge of each trip should come from the init intersection for setup {}".format(test_setup.name))
    number_of_empty_routes = 0
    for person_id in test_setup.people_to_test:
        init_intersection_according_to_people_file = \
            int(test_setup.pd_people[test_setup.pd_people['p'] == person_id]['init_intersection'].iloc[0])
        route_list = route_csv_string_to_list(str(test_setup.pd_route[test_setup.pd_route['p'] == person_id]['route'].iloc[0]))
        if route_list == []:
            number_of_empty_routes += 1
            continue
        first_edge = int(route_list[0])
        init_intersection_according_to_route_and_edges_csv = \
            int(test_setup.pd_edges[test_setup.pd_edges['uniqueid'] == first_edge]['osmid_u'])

        try:
            assert init_intersection_according_to_people_file == init_intersection_according_to_route_and_edges_csv, \
                "Initial intersections do not match! According to the peoples' file, init intersection is " + \
                str(init_intersection_according_to_people_file) + \
                " but according to the routes' and edges files, the init intersection is " + \
                str(init_intersection_according_to_route_and_edges_csv) + \
                " for person " + str(person_id)
        except:
            st()

    log("Number of empty routes skipped: {}".format(number_of_empty_routes))

    log("Passed test 03")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test04_last_edge_of_trip_should_go_to_end_intersection(run_all_test_setups, test_setup):
    log("Test 04: The last edge of each trip should end at the end intersection for setup {}".format(test_setup.name))
    number_of_empty_routes = 0
    for person_id in test_setup.people_to_test:
        end_intersection_according_to_people_file = \
            int(test_setup.pd_people[test_setup.pd_people['p'] == person_id]['end_intersection'].iloc[0])
        route_list = route_csv_string_to_list(str(test_setup.pd_route[test_setup.pd_route['p'] == person_id]['route'].iloc[0]))
        if route_list == []:
            number_of_empty_routes += 1
            continue

        last_edge = int(route_list[-1])
        end_intersection_according_to_route_and_edges_csv = \
            int(test_setup.pd_edges[test_setup.pd_edges['uniqueid'] == last_edge]['osmid_v'])

        assert end_intersection_according_to_people_file == end_intersection_according_to_route_and_edges_csv, \
            "End intersections do not match! According to the peoples' file, end intersection is " + \
            str(end_intersection_according_to_people_file) + \
            " but according to the routes' and edges files, the end intersection is " + \
            str(end_intersection_according_to_route_and_edges_csv) + \
            " for person " + str(person_id)
    log("Number of empty routes skipped: {}".format(number_of_empty_routes))
    log("Passed test 04")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_05_each_edge_of_trip_should_connect_to_the_following_edge(run_all_test_setups, test_setup):
    log("Test 05: Each edge of each trip should connect to the following edge for setup {}".format(test_setup.name))
    number_of_empty_routes = 0
    for person_id in tqdm(test_setup.people_to_test):
        route_list = route_csv_string_to_list(str(test_setup.pd_route[test_setup.pd_route['p'] == person_id]['route'].iloc[0]))
        if route_list == []:
            number_of_empty_routes += 1
            continue

        for first_edge_index in range(len(route_list)-1):
            first_edge = int(route_list[first_edge_index])
            first_edge_end_intersection = int(test_setup.pd_edges.loc[test_setup.pd_edges['uniqueid'] == first_edge]['osmid_v'])

            second_edge = int(route_list[first_edge_index + 1])
            second_edge_init_intersection = int(test_setup.pd_edges.loc[test_setup.pd_edges['uniqueid'] == second_edge]['osmid_u'])
            
            error_message = "route for person {} does is not a valid path. Route is {} but edge {} and {} are not" + \
                "connected in the network. Edge {} finishes at intersection {} but edge {} starts at intersection {}" \
                .format(person_id, route_list, first_edge, second_edge, first_edge,
                first_edge_end_intersection, second_edge, second_edge_init_intersection)
            assert first_edge_end_intersection == second_edge_init_intersection, error_message
    log("Number of empty routes skipped: {}".format(number_of_empty_routes))
    log("Passed test 05")

# An empty route at the route file and a nan avg speed at the people file both mean that the trip has not finished
@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_06_empty_routes_should_match_nan_avg_speed_and_same_origin_and_destination(run_all_test_setups, test_setup):
    log("Test 06: People with empty routes should have nan avg speed for test setup {}".format(test_setup))
    people_with_empty_routes = test_setup.pd_route[test_setup.pd_route['route'] == '[]']['p']
    log("Found {} people with empty routes".format(len(people_with_empty_routes)))

    people_with_nan_avg_speed = test_setup.pd_people[np.isnan(test_setup.pd_people['avg_v(mph)'])]['p']
    log("Found {} people with NaN avg speed".format(len(people_with_nan_avg_speed)))

    try:
        assert people_with_empty_routes.equals(people_with_nan_avg_speed), 'The people with empty routes are not ' + \
                'the same people that have NaN avg_v(mph). They should be the same.'
    except:
        st()
    log("They match.")

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def deprecated_test_07_first_edge_of_trip_should_be_the_correct_initial_edge(run_all_test_setups, test_setup):
    log("Testing if the first edge of each trip is the correct initial edge...")
    number_of_empty_routes = 0
    for person_id in test_setup.people_to_test:
        if person_id == 0:
            continue
        route_list = route_csv_string_to_list(str(test_setup.pd_route.iloc[person_id][1]))
        if route_list == []:
            number_of_empty_routes += 1
            continue
        
        init_edge_according_to_route_file = int(route_list[0])

        init_edge_according_to_person_to_init_edge = int(test_setup.pd_indexPathInit[test_setup.pd_indexPathInit['p'] == person_id]['indexPathInit'].iloc[0])

        assert init_edge_according_to_route_file == init_edge_according_to_person_to_init_edge, \
            "Initial edges do not match! For person " + str(person_id) + ", according to the routes' file, init edge is " + \
            str(init_edge_according_to_route_file) + \
            " but according to the person to init file, the init edge is " + \
            str(init_edge_according_to_person_to_init_edge)
    log("Number of empty routes skipped: {}".format(number_of_empty_routes))

    log("Passed")

@pytest.mark.parametrize("test_setup", pytest.test_setups)
def test_08_indexPathInit_should_have_all_different_values(run_all_test_setups, test_setup):
    log("Testing if indexPathInit has all different values...")

    assert not test_setup.pd_indexPathInit['indexPathInit'].duplicated().any(), \
        "All indexPathInit values should be different"

    log("Passed")

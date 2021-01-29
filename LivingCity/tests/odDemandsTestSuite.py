import math
from utils import *



""" We create a random od_demand file where each trip has 1 edge that uniquely connects two nodes.
    Then corroborate that both the people and route files follow the generated demands.
    The filename has a 3 digit random number appended to it, to check if the parametrization works properly.
"""
def randomly_obtain_edges_that_uniquely_connect_two_nodes(df_edges, number_of_edges):
    result_edges = []
    iteration_order = list(range(len(df_edges)))
    random.shuffle(iteration_order)
    for i in iteration_order:
        an_edge = df_edges.iloc[i]

        if an_edge['uniqueid'] in [edge['uniqueid'] for edge in result_edges]:
            continue

        if len(df_edges[(df_edges['osmid_u'] == an_edge['osmid_u']) & (df_edges['osmid_v'] == an_edge['osmid_v'])]) == 1:
            result_edges.append(an_edge)

        if len(result_edges) == number_of_edges:
            break
    return result_edges

def test_01_randomly_generated_demands_of_one_edge_trips_corresponds_to_people_and_route_files():
    min_dep_time = 16000
    max_dep_time = 20000
    number_of_trips = 20

    network_path = "berkeley_2018/new_full_network/"
    log("Running system test setup on output files for network {}...".format(network_path))
    output_files = ["route", "people", "indexPathVec"]

    log("Removing previous run files...")
    for output_file_name in output_files:
        for run in ["", "demands_test"]:
            filename = "./0_{}5to12{}.csv".format(output_file_name, run)
            if os.path.exists(filename):
                log("Removing old {}...".format(filename))
                os.remove(filename)

    # Create the od_demand file
    od_demand_test_list = ['SAMPN,PERNO,origin_osmid,destination_osmid,dep_time,origin,destination']
    df_edges = pd.read_csv(os.path.join(network_path, 'edges.csv'))
    log("Randomly generating demands...")
    edges_for_trips = randomly_obtain_edges_that_uniquely_connect_two_nodes(df_edges, number_of_trips)
    for SAMPN, current_trip_edge in enumerate(edges_for_trips):
        PERNO = 1
        origin_osmid = current_trip_edge['osmid_u']
        destination_osmid = current_trip_edge['osmid_v']

        dep_time = random.randint(min_dep_time, max_dep_time)

        origin = current_trip_edge['u']
        destination = current_trip_edge['v']

        current_trip = [SAMPN, PERNO, origin_osmid, destination_osmid, dep_time, origin, destination]
        od_demand_test_list.append(",".join(map(str,map(int,current_trip))))

    od_demand_test = '\n'.join(od_demand_test_list)

    od_demand_test_filename = 'od_demand_test_{}.csv'.format(random.randint(100, 999))
    with open(os.path.join(network_path, od_demand_test_filename), 'w') as od_demand_file:
        od_demand_file.write(od_demand_test)

    log("Generated the following random demands:")
    log(od_demand_test, color='green')

    log("Running simulation with od_demand_test as the od_demand_file...")
    write_options_file({"NETWORK_PATH": "{}".format(network_path), \
                        "OD_DEMAND_FILENAME": od_demand_test_filename, \
                        "REROUTE_INCREMENT": 0})
    subprocess.run(["./LivingCity", "&"], check=True)

    output_files = ["route", "people", "indexPathVec"]
    for output_file_name in output_files:
        os.rename("./0_{}5to12.csv".format(output_file_name), "./0_{}5to12_od_demand_test.csv".format(output_file_name))
    log("Finished setup.")

    log("Checking that the people and route files follow the generated demands...")
    df_people = pd.read_csv("0_people5to12_od_demand_test.csv")
    df_route = pd.read_csv("0_route5to12_od_demand_test.csv", sep = ':')

    df_edges_for_trips = pd.DataFrame(edges_for_trips)
    for a_person_index, a_person in df_people.iterrows():
        edge = df_edges_for_trips[df_edges_for_trips['osmid_u'] == int(a_person['init_intersection'])]

        assert all(edge['osmid_v'] == int(a_person['end_intersection'])), \
            'Destination for person {} does not match the inputted od_demand file. '.format(a_person_index) + \
            'od_demand file saved at {}.'.format(os.path.join(network_path, od_demand_test_filename))

        a_person_id = int(a_person['p'])
        route_according_to_route_file = df_route.loc[df_route['p'] == a_person_id]['route'].iloc[0]
        route_according_to_input_od_demands = '[{},]'.format(int(edge['uniqueid'].iloc[0]))
        assert route_according_to_route_file == route_according_to_input_od_demands, \
            'For person {} input od_demands specify a trip of {} while the route has a trip of {}.'.format( \
            a_person_id, route_according_to_input_od_demands, route_according_to_route_file)

        assert math.isclose(a_person['distance'],edge['length']), \
            'Distance for person {} does not match the inputted od_demand file. '.format(a_person_index) + \
            'od_demand file saved at {}.'.format(os.path.join(network_path, od_demand_test_filename))

    log("Removing {}...".format(od_demand_test_filename))
    os.remove(os.path.join(network_path, od_demand_test_filename))
    log("Passed.")

    return network_path
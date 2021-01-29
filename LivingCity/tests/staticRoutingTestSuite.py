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
pytest.number_of_people = 3441952
pytest.distance_margin_between_route_and_people_file = 100
pytest.edges_path = "berkeley_2018/new_full_network/edges.csv"

"""
default_network_setup runs automatically before the tests that take default_network_setup as a parameter
For each network specified it does the following:
    * Runs LivingCity with prev_paths = false and then with prev_paths = true
    * Produces 0_people5to12_first_run.csv, 0_indexPathVec5to12_first_run.csv, 0_route5to12_first_run.csv on the first run
    * Produces 0_people5to12_second_run.csv, 0_indexPathVec5to12_second_run.csv, 0_route5to12_second_run.csv on the second run
"""
@pytest.fixture(params=["berkeley_2018/new_full_network/"], ids=["new_full_network"])
def default_network_setup(request):
    log("Running default network setup...")
    network_path = pytest.network_path
    if pytest.default_network_setup_has_run:
        return network_path
    pytest.default_network_setup_has_run = True

    log("Running system test setup on output files for network {}".format(network_path))
    output_files = ["route", "people", "indexPathVec"]

    log("Removing previous run files...")
    for output_file_name in output_files:
        for run in ["", "_first_run", "_second_run"]:
            filename = "./0_{}5to12{}.csv".format(output_file_name, run)
            if os.path.exists(filename):
                log("Removing old {}...".format(filename))
                os.remove(filename)

    log("Running first simulation with use_prev_paths=false...")
    write_options_file({"NETWORK_PATH": "{}".format(
        network_path), "USE_PREV_PATHS": "false", "REROUTE_INCREMENT": 0})
    subprocess.run(["./LivingCity", "&"], check=True)

    for output_file_name in output_files:
        os.rename("./0_{}5to12.csv".format(output_file_name), "./0_{}5to12_first_run.csv".format(output_file_name))

    log("Running second simulation with use_prev_paths=true...")
    write_options_file({"NETWORK_PATH": "{}".format(
        network_path), "USE_PREV_PATHS": "true", "REROUTE_INCREMENT": 0})
    subprocess.run(["./LivingCity", "&"], check=True)
    for output_file_name in output_files:
        os.rename("./0_{}5to12.csv".format(output_file_name), "./0_{}5to12_second_run.csv".format(output_file_name))

    log("Finished setup.")

    return network_path


def test01_all_output_files_should_have_the_same_length(default_network_setup):
    lengths = []
    log("Checking that the length of people and route files are the same...")
    for output_file in ["route", "people"]:
        delimiter = ":" if output_file == "route" else ","
        lengths.append(length_of_csv(
            '0_{}5to12_first_run.csv'.format(output_file), delimiter))
    assert all([lengths[0] == elem for elem in lengths])
    pytest.number_of_people = lengths[0]
    log("Passed")


@pytest.mark.parametrize("file_to_preserve", ["route", "indexPathVec"])
def test02_prev_paths_should_preserve_route_and_indexPathVec_files(default_network_setup, file_to_preserve):
    log("Comparing {} files between the two runs...".format(file_to_preserve))
    assert filecmp.cmp("0_{0}5to12_first_run.csv".format(file_to_preserve), "0_{0}5to12_second_run.csv".format(
        file_to_preserve), shallow=False), "{} file are not equal between runs".format(file_to_preserve)
    log("Passed")


def test03_prev_paths_should_have_consistent_people_files(default_network_setup):
    log("Comparing people files between the two runs...")

    df_first_run = pd.read_csv("0_people5to12_first_run.csv")
    df_second_run = pd.read_csv("0_people5to12_second_run.csv")
    for (id_row, row_first_run), (_, row_second_run) in tqdm(zip(df_first_run.iterrows(), df_second_run.iterrows()),
                                                             total=pytest.number_of_people):
        parameters_that_should_be_equal = [
            "p", "init_intersection", "end_intersection", "time_departure", "a", "b", "T"]
        for param in parameters_that_should_be_equal:
            assert row_first_run[param] == row_second_run[param], \
                "{} parameter is not equal in the people file between runs for row {}.".format(param, id_row) +\
                "First column has {} while second column has {}.".format(
                row_first_run[param], row_second_run[param])
    log("Passed")


def test_04_distance_in_people_file_should_match_the_sum_of_the_edges_in_the_route_file(default_network_setup):
    log("Comparing that the distance in the people file matches the sum of the edges in the route file with a margin of {}...".format(
        pytest.distance_margin_between_route_and_people_file))
    log("(Since testing for > 2 million people takes several hours, it is tested on 1.000 random people).")
    pd_people = pd.read_csv("0_people5to12_first_run.csv")
    pd_edges = pd.read_csv(pytest.edges_path)

    people_to_test = [random.randint(
        0, pytest.number_of_people) for i in range(1000)]

    for chunk_route in tqdm(pd.read_csv("0_route5to12_first_run.csv", sep=":", chunksize=1000),
                            total=pytest.number_of_people/1000):
        for _, row in chunk_route.iterrows():
            person_id = str(row["p"])

            if (not person_id in people_to_test):
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
            assert abs(distance_people_info -
                       distance_sum_of_edges) < pytest.distance_margin_between_route_and_people_file, error_message
    log("Passed")
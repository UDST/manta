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

@pytest.mark.parametrize("file_to_preserve", ["route", "indexPathVec"])
def test01_prev_paths_should_preserve_route_and_indexPathVec_files(default_network_setup, file_to_preserve):
    log("Comparing {} files between the two runs...".format(file_to_preserve))
    assert filecmp.cmp("0_{0}5to12_first_run.csv".format(file_to_preserve), "0_{0}5to12_second_run.csv".format(
        file_to_preserve), shallow=False), "{} file are not equal between runs".format(file_to_preserve)
    log("Passed")


def test02_prev_paths_should_have_consistent_people_files(default_network_setup):
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

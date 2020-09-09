import pandas as pd
from pdb import set_trace as st
import config

"""
  Given a csv with 3 columns (person_id,distance_sum_of_edges,distance_people_info),
  calculates the intersection between the two sets of distances distance_sum_of_edges and distance_people_info.
  This script was wrote in order to debug the discrepancy of distances between the distance according to the sum of
  the edges travelled by a person and the distance according to the output of LivingCity.

  Args:
    distances_file: Path to the csv file. Default given by config.py
      
    filter_out_distance_0: If true, filters out matches with distance 0.0. Useful since it's
      likely a match at a distance of exactly 0 is more of a coincidence than an index mismatch.
  
  Returns: 
    Nothing
"""

def intersect_distances_routes_vs_people(distances_file=config.distance_merge_file):
  print("Reading distance files...")
  distances = pd.read_csv(distances_file)

  pd_distances_sum_of_edges = distances[['person_id', 'distance_sum_of_edges']]
  len_distances_sum_of_edges_equal_to_0 = len(pd_distances_sum_of_edges[pd_distances_sum_of_edges.distance_sum_of_edges == 0.0])
  pd_distances_people_info = distances[['person_id', 'distance_people_info']]
  len_distances_people_info_equal_to_0 = len(pd_distances_people_info[pd_distances_people_info.distance_people_info == 0.0])

  print("Calculating inner join between distance_sum_of_edges and distance_people_info...")
  print("""Specifically:
        SELECT d_people_info.person_id as person_id, d_people_info.distance_people_info as distance_matched
        FROM distances as d_sum_of_edges
        INNER JOIN distances as d_people_info ON d_sum_of_edges.distance_sum_of_edges = d_people_info.distance_people_info;\n""")

  distances_inner_join = pd.merge(pd_distances_sum_of_edges,
                                  pd_distances_people_info,
                                  left_on='distance_sum_of_edges',
                                  right_on='distance_people_info')
  distances_inner_join = distances_inner_join.drop(["distance_people_info"],axis=1)

  distances_inner_join_people_info = group_by(distances_inner_join,"person_id_people_info")
  print("Distances inner joined grouped by the person id of the people info:")
  print(distances_inner_join_people_info)
  print("Conclusion: {} out of {} distances from the people information file match with a distance calculated by the sum of the edges."\
        .format(len(distances_inner_join_people_info), len(distances)))
  print("{} out of those matches have a distance of 0.".format(len_distances_sum_of_edges_equal_to_0))

  print("---")

  distances_inner_join_sum_of_edges = group_by(distances_inner_join,"person_id_sum_of_edges")
  print("Distances inner joined grouped by the person id of the sum of the edges:")
  print(distances_inner_join_sum_of_edges)
  print("Conclusion: {} out of {} distances calculated from the sum of the edges match with a distance obtained from the people information file."\
        .format(len(distances_inner_join_sum_of_edges), len(distances)))
  print("{} out of those matches have a distance of 0.".format(len_distances_people_info_equal_to_0))

def group_by(distances_inner_join, by):
  distances_inner_join.columns = ["person_id_sum_of_edges", "distance_matched", "person_id_people_info"]
  distances_inner_join_grouped = distances_inner_join.groupby([by],as_index = False).count()
  distances_inner_join_grouped = distances_inner_join_grouped.drop(["distance_matched"],axis=1)
  distances_inner_join_grouped.columns = [by,"number_of_matches"]

  return distances_inner_join_grouped


if __name__ == '__main__':
  intersect_distances_routes_vs_people()

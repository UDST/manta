import pandas as pd
from tqdm import tqdm
from pdb import set_trace as st
import os,sys
import config

"""
  Given a csv with 3 columns (person_id,distance_sum_of_edges,distance_people_info),
  calculates the intersection between the two sets of distances distance_sum_of_edges and distance_people_info.
  This script was wrote in order to debug the discrepancy of distances between the distance according to the sum of
  the edges travelled by a person and the distance according to the output of LivingCity.

  Args:
    distances_file: Path to the csv file.

    only_show_intersection_size: If true, output the intersection size to stdout.
        If false, output the size to stdout while outputting the intersection itself to a output_intersection_file
      
    output_intersection_file: (Optional) File to output the intersection of the distances.
  
  Returns: 
    Nothing
"""

def intersect_distances_routes_vs_people(distances_file=config.DEFAULT_DISTANCE_MERGE_OUTPUT_FILE_PATH, \
      only_show_intersection_size=False, output_intersection_file = config.DEFAULT_DISTANCE_INTERSECTION_FILE_PATH):

  print("Reading distance files...")
  distances = pd.read_csv(distances_file)

  pd_distances_sum_of_edges = distances[['person_id', 'distance_sum_of_edges']]
  pd_distance_people_info = distances[['person_id', 'distance_people_info']]

  print("Calculating inner join between distance_sum_of_edges and distance_people_info...")
  print("""Specifically:
        SELECT d_people_info.person_id as person_id, d_people_info.distance_people_info as distance_matched
        FROM distances as d_sum_of_edges
        INNER JOIN distances as d_people_info ON d_sum_of_edges.distance_sum_of_edges = d_people_info.distance_people_info;\n""")

  distances_inner_join = pd.merge(pd_distances_sum_of_edges,
                                  pd_distance_people_info,
                                  left_on='distance_sum_of_edges',
                                  right_on='distance_people_info')
  distances_inner_join = distances_inner_join.drop(["distance_people_info"],axis=1)
  
  distances_inner_join.columns = ["person_id_sum_of_edges", "distance_matched", "person_id_people_info"]
 
  print("Number of people ") = distances_inner_join.groupby(['distance_people_info']).count()
  

if __name__ == '__main__':
  intersect_distances_routes_vs_people()

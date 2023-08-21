import os
import shutil
import sys

# Add parent directory to sys.path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pandas as pd
import json
import tempfile
import argparse
import warnings
import src.predictor_utils as utils
from Bio.PDB import PDBIO, is_aa, PDBParser
from src.tmalign import Tmalign

# Ignore warnings
warnings.filterwarnings("ignore")

####################################################################################################################

if __name__ == '__main__':
    # Create an argument parser object
    parser = argparse.ArgumentParser(description='Generate updated collection for entries')

    # Add command line arguments
    parser.add_argument('--in_path', '-i', type=str,
                        help='Path to input JSON file')
    parser.add_argument('--out_dir', '-o', type=str,
                        help='Path to output dir')
    parser.add_argument('--target_db', '-t', type=str,
                        help='Path to target db')
    parser.add_argument('--query_db', '-q', type=str,
                        help='Path to query db')

    # Parse command line arguments
    args = parser.parse_args()

    # Access the variable values
    in_path = args.in_path
    out_dir = args.out_dir
    target_db = args.target_db
    query_db = args.query_db

    # Instantiate essential modules
    pdb_parser = PDBParser(QUIET=True)
    io_handler = PDBIO()
    tm_align = Tmalign()

    ####################################################################################################################

    # Read the input JSON file
    with open(in_path, 'r') as fp:
        input_dict = json.load(fp)

    # Load the essential data
    query_name = input_dict["query"]
    target_name = input_dict["target"]

    # Specify the path to the query and target structures
    query_structure_path = os.path.join(query_db, query_name[1:3], f"{query_name[:5]}.pdb")
    target_structure_path = os.path.join(target_db, target_name[1:3], f"{target_name}_tmax.pdb")

    ####################################################################################################################

    # Parse the structure from the PDB file
    qstructure = pdb_parser.get_structure(query_name, query_structure_path)
    # Designate the chain from the structure
    qchain = qstructure[0][query_name[4]]
    # Extract residues from the chain
    qchain_residues = [residue for residue in qchain.get_residues() if is_aa(residue.get_resname())]

    # Do the same for the target structure
    tstructure = pdb_parser.get_structure(target_name, target_structure_path)
    tchain = tstructure[0][target_name[4]]
    tchain_residues = [residue for residue in tchain.get_residues() if is_aa(residue.get_resname())]

    ####################################################################################################################

    # Specify the step of framing the fragments
    # (frame_step = 1 --> moving the frame 1 residue per time --> highest resolution)
    frame_step = 1
    fragment_frame_list = utils.get_res_frames(res_range=qchain_residues, length=len(tchain), step=frame_step)

    # Define the temporary dir to save the fragment structures
    fragment_dir = tempfile.mkdtemp()

    # Create a list to store the path to fragment structures
    fragment_path_list = []
    # Iterate the fragment frame list, extract and store each fragment structure
    for fragment in fragment_frame_list:
        fragment_start = fragment[0].id[1]
        fragment_end = fragment[-1].id[1]
        fragment_out_path = os.path.join(fragment_dir, f"{query_name}_{fragment_start}_{fragment_end}.pdb")
        utils.get_structure(res_range=fragment, res_chain=query_name[4], structure=qstructure,
                            out_path=fragment_out_path, io_handler=io_handler)
        fragment_path_list.append(fragment_out_path)

    ####################################################################################################################

    # Create a list to store tm-scores
    tmscore_results = []
    # Iterate over fragment structure paths
    for fragment_path in fragment_path_list:
        # Identify the first residue of the fragment structure (x)
        fragment_start = os.path.basename(fragment_path).split("_")[1]
        # Align the target and the query fragment structure
        command_output = tm_align(target_structure_path, fragment_path)
        # Parse and extract the normalized tmscore (y)
        lines = command_output.split("\n")
        tm_score = None
        for line in lines:
            if line.startswith('TM-score=') and 'if normalized by average length' in line:
                # Extract the TM-score value
                tm_score = float(line.split('=')[1].split()[0])
                # Stop looping after finding the desired TM-score line
                break

        # Add the results to the list as a pair
        tmscore_results.append([fragment_start, tm_score])

    ####################################################################################################################

    # Sort the result list based on the query fragment start residue (x)
    tm_score_results = [[int(i[0]), i[1]] for i in tmscore_results]
    tm_score_results.sort(key=lambda x: x[0])

    # Divide and store the sorted data
    x = [i[0] for i in tm_score_results]
    y = [i[1] for i in tm_score_results]

    # Save the results a csv file at the output directory
    graph_df = pd.DataFrame({"x": x, "y": y})

    out_subdir = os.path.join(out_dir, query_name[1:3])
    os.makedirs(out_subdir, exist_ok=True)
    out_path = os.path.join(out_subdir, f"{query_name}_vs_{target_name}_graph.csv")
    graph_df.to_csv(out_path)

    # Remove the temporary directory
    shutil.rmtree(fragment_dir)

    # Finish

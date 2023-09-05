import os
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import warnings
import src.utils as utils

warnings.filterwarnings("ignore")

# Load the RepeatsDB dataframe
repeatsdb_df = pd.read_csv("../data/repeatsdb_supplemented_updated.csv", dtype={"ct": str})
# Create a list of region_id's with more than one entry with same pdb and class.topology (ct)
multi_region_list = repeatsdb_df.groupby(['pdb', 'ct']).filter(lambda group: len(group) > 1)["region_id"].tolist()
# Specify the input and output directories
in_dir = "/home/soroushm/Documents/ML-TRP-Parser_workframe/dataframes/graph_tmscore_output"
out_dir = "/home/soroushm/Documents/ML-TRP-Parser_workframe/dataframes/graph_features"

# Iterate the files in the input directory and extract the file paths
file_path_list = []
for root, dirs, files in os.walk(in_dir):
    for filename in files:
        if filename.endswith("_graph.csv"):
            file_path = os.path.join(root, filename)
            file_path_list.append(file_path)

# Iterate over each file path and plot, annotated, save the graph
for file_path in file_path_list:
    try:
        # Identify the region_id and load the associated data
        region_id = "_".join(os.path.basename(file_path).split("_")[:3])
        graph_path = os.path.join(in_dir, f"{region_id[1:3]}", f"{region_id}_vs_{region_id}_graph.csv")
        graph_df = pd.read_csv(graph_path)

        # If the region_id in multi_region_list, limit the graph data to the region_id start and end positions
        if region_id in multi_region_list:
            region_start = int(region_id.split("_")[1])
            region_end = int(region_id.split("_")[2])
            graph_df = graph_df[region_start <= graph_df["x_start"]]
            graph_df = graph_df[region_end >= graph_df["x_start"]]

        # Extract x, y and length from the graph dataframe
        # Convert into numpy arrays
        x = np.array(graph_df["x_start"])
        y = np.array(graph_df["y"])

        # Smooth the graph by n residue window size
        # (small size --> noise reduction & not compromising resolution)
        window_len = graph_df.iloc[0, 1] - graph_df.iloc[0, 0] + 1  # Last residue is inclusive
        smooth_window_size = round(window_len * 0.2)
        if smooth_window_size % 2 == 0:
            smooth_window_size -= 1
        poly_order = 0
        y = savgol_filter(y, smooth_window_size, poly_order)

        x, y = utils.adjust_graph_ends(x, y)

        # Detect peaks and associated properties
        peaks, properties = find_peaks(y, height=0.3, distance=round(window_len * 0.6),
                                       threshold=0, prominence=0, width=0)

        # Mark the peaks on y and x-axis
        peaks_x = [x[peak_idx] for peak_idx in peaks]
        peaks_y = [y[peak_idx] for peak_idx in peaks]

        peaks_df = pd.DataFrame(properties)

        # Extract data of the reference units
        ref_units = repeatsdb_df[repeatsdb_df["region_id"] == region_id]["units"].values[0].split(",")
        ref_unit_start_list = [int(i.split("_")[0]) for i in ref_units]

        # Create a list of reference unit accepted margin range
        unit_start_range_list = []
        margin = round(window_len * 0.2)
        for ref_unit_start in ref_unit_start_list:
            range_start = ref_unit_start - margin
            range_end = ref_unit_start + margin
            unit_start_range = (range_start, range_end)
            unit_start_range_list.append(unit_start_range)

        # Check if the detected peaks fall into accepted margin range
        peak_labels = []
        for x_start in peaks_x:
            truth = 0
            for unit_start_range in unit_start_range_list:
                range_start, range_end = unit_start_range
                if range_start <= x_start <= range_end:
                    truth = 1
            peak_labels.append(truth)

        # Extract extra data regarding region_id
        ct = str(repeatsdb_df[repeatsdb_df["region_id"] == region_id]["ct"].values[0])
        window_avg = repeatsdb_df[repeatsdb_df["region_id"] == region_id]["units_avg"].values[0]

        # Add extra data to the peaks dataframe
        peaks_df["peak_x"] = peaks_x
        peaks_df["ct"] = [ct for i in range(len(peaks_df))]
        peaks_df["window_avg"] = [window_avg for i in range(len(peaks_df))]
        peaks_df["window_len"] = [window_len for i in range(len(peaks_df))]
        peaks_df["label"] = peak_labels

        # Save the results as CSV at the output directory
        out_subdir = os.path.join(out_dir, region_id[1:3])
        os.makedirs(out_subdir, exist_ok=True)
        out_path = os.path.join(out_subdir, f"{region_id}_features.csv")
        peaks_df.to_csv(out_path, index=False)

    # In case of error, print the error and the file path causing it
    except Exception as error:
        print(error)
        print(file_path)

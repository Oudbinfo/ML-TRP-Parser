import os
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from src.predictor_utils import adjust_graph_ends
import warnings

warnings.filterwarnings("ignore")

# Define the essential inputs
repeatsdb_df = pd.read_csv("../data/repeatsdb_unique_updated.csv", dtype={'ct': str})
in_dir = "/home/soroushm/Documents/ML-TRP-Parser_workframe/dataframes/graph_tmscore_output"
out_dir = "/home/soroushm/Documents/ML-TRP-Parser_workframe/dataframes/graph_features"

# Iterate the files in the input directory and extract the file paths
file_path_list = []
for root, dirs, files in os.walk(in_dir):
    for filename in files:
        file_path = os.path.join(root, filename)
        file_path_list.append(file_path)

# Iterate over each file path and plot, annotated, save the graph
for file_path in file_path_list:
    try:
        # Identify the region id and load the associated data
        region_id = "_".join(os.path.basename(file_path).split("_")[:3])
        graph_df_path = os.path.join(in_dir, f"{region_id[1:3]}", f"{region_id}_vs_{region_id}_graph.csv")
        graph_df = pd.read_csv(graph_df_path)

        # Extract x, y and length from the graph dataframe
        # Convert into numpy arrays
        x = np.array(graph_df["x_start"])
        y = np.array(graph_df["y"])
        length = graph_df.loc[0, "x_end"] - graph_df.loc[0, "x_start"] + 1  # Last residue is inclusive

        # Extract the data regarding the annotated units and their starting point
        units = repeatsdb_df[repeatsdb_df["region_id"] == region_id]["units"].values[0].split(",")
        unit_start_list = [int(i.split("_")[0]) for i in units]

        ct = str(repeatsdb_df[repeatsdb_df["region_id"] == region_id]["ct"].values[0])

        # Smooth the graph by n residue window size
        # (small window size --> reducing noise / not compromising resolution)
        window_size = round(length * 0.2)
        if window_size % 2 == 0:
            window_size -= 1
        poly_order = 0
        y = savgol_filter(y, window_size, poly_order)

        x, y = adjust_graph_ends(x, y)

        # Detect peaks and associated properties
        peaks, properties = find_peaks(y, height=0.2, distance=round(length * 0.5),
                                       threshold=0, prominence=0, width=0, plateau_size=0)

        # Mark the peaks on y and x-axis
        peaks_x = [x[peak_idx] for peak_idx in peaks]
        peaks_y = [y[peak_idx] for peak_idx in peaks]

        peaks_df = pd.DataFrame(properties)
        peaks_df["peak_x"] = peaks_x
        column_list = ["left_edges", "right_edges", "left_ips", "right_ips", "left_bases", "right_bases"]
        for column in column_list:
            peaks_df[column] = (peaks_df['peak_x'] - peaks_df[column]) / length

        # Create shifted dataframes
        peaks_df_previous1 = peaks_df.shift(1)
        peaks_df_next1 = peaks_df.shift(-1)

        # Rename columns
        peaks_df_previous1.columns = ['pre1_' + col for col in peaks_df.columns]
        peaks_df_next1.columns = ['post1_' + col for col in peaks_df.columns]

        # Concatenate
        peaks_df = pd.concat([peaks_df_previous1, peaks_df, peaks_df_next1], axis=1)
        peaks_df.fillna(0, inplace=True)

        unit_start_range_list = []
        margin = round(length * 0.2)
        for unit_start in unit_start_list:
            range_start = unit_start - margin
            range_end = unit_start + margin
            unit_start_range = (range_start, range_end)
            unit_start_range_list.append(unit_start_range)


        def is_in_range(number, range_tuple):
            start, end = range_tuple
            return start <= number <= end


        peak_labels = []
        for x_start in peaks_x:
            truth = 0
            for unit_start_range in unit_start_range_list:
                if is_in_range(x_start, unit_start_range):
                    truth = 1
            peak_labels.append(truth)

        peaks_df["peak_x"] = peaks_x
        peaks_df["ct"] = [ct for i in range(len(peaks_df))]
        peaks_df["window_len"] = [length for i in range(len(peaks_df))]
        peaks_df["label"] = peak_labels

        peaks_df.drop(["pre1_peak_x", "peak_x", "post1_peak_x"], axis=1, inplace=True)

        out_subdir = os.path.join(out_dir, region_id[1:3])
        os.makedirs(out_subdir, exist_ok=True)
        out_path = os.path.join(out_subdir, f"{region_id}_features.csv")
        peaks_df.to_csv(out_path, index=False)

    # In case of error, print the error and the file path causing it
    except Exception as error:
        print(error)
        print(file_path)
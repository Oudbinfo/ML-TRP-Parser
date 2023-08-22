import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from src.predictor_utils import adjust_graph_ends
import warnings

warnings.filterwarnings("ignore")

# Define the essential inputs
repeatsdb_df = pd.read_csv("../data/repeatsdb_unique_updated.csv", dtype={'ct': str})
in_dir = "/home/soroushm/Documents/ML-TRP-Parser_workframe/dataframes/graph_tmscore_output"
out_dir = "/home/soroushm/Documents/ML-TRP-Parser_workframe/graphs"

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
        unit_starts = [int(i.split("_")[0]) for i in units]

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

        # Plot the results
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(12, 7))
        # Graph of all the residue vs tm-score points
        ax.plot(x, y)
        # Scatter of all the residue vs tm-score points
        ax.scatter(x, y, marker='o', color='blue')
        # Scatter of the detected peaks
        ax.scatter(peaks_x, peaks_y, marker="*", color="red", s=100)
        for unit_start in unit_starts:
            ax.axvline(x=unit_start, color='r', linestyle='--')

        # Highlight a margin of 7 residues per each unit start annotation
        lim_1 = ax.get_ylim()[1]
        lim_2 = ax.get_ylim()[0]
        highlight_range = round(length * 0.2)
        for point in peaks_x:
            start_x = max(point - highlight_range, min(x))
            end_x = min(point + highlight_range, max(x))
            ax.fill_between([start_x, end_x], lim_1, lim_2, color='green', alpha=0.3)

        # Save the graph at the output directory
        out_subdir = os.path.join(out_dir, region_id[1:3])
        os.makedirs(out_subdir, exist_ok=True)
        out_path = os.path.join(out_subdir, f"{region_id}.png")
        plt.savefig(out_path, format="png")
        plt.close()

    # In case of error, print the error and the file path causing it
    except Exception as error:
        print(error)
        print(file_path)

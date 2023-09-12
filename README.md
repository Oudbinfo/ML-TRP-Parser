# Machine Learning for Tandem Repeat Proteins Mapping
This repository is focused on automating the identification and mapping of tandem repeat proteins using machine learning techniques.

# Modules Overview:

## graph_tmscore.py & graph_tmscore.sh:
Aligns target structures with query structures in sliding window fashion.
Generates graphs of structural similarities at each query position.
Bash script provides compatibility with SGE systems.

## extract_graph_features.ipynb & extract_graph_features.py:
Detects and analyses peaks on the generated graphs.
Extracts features from peaks and labels them based on RepeatsDB reference data.

## evaluate_graph_features.ipynb:
Assesses the quality of the extracted features.
Uses Logistic Regression and Random Forest Classifier to develop a predictive model for peak detection.

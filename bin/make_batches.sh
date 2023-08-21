#!/bin/bash


# Define input path
in_path=$1
# Define output path
out_path=$2
# Define file basename_regex
basename_regex=$3
# Define batch size
batch_size=$4
# Define whole files list
input_list=$(find "$(realpath "${in_path}")" -name "${basename_regex}")
# Define output directory
mkdir -p "${out_path}"
# Make batches of initial list
echo "${input_list}" | awk -v subs=${batch_size} -v out="${out_path}" '{if (((FNR-1) % subs)==0) {fn=fn+1}; print $0 >> ((out "/list." fn ".dat"))}'
# Define list of lists
find "$(realpath "${out_path}")" -name "list.*.dat" >"${out_path}/lists.dat"

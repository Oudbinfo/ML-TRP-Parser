#!/bin/bash
#$ -o /home/$USER/logs/
#$ -e /home/$USER/logs/
#$ -cwd
#$ -q high@bronte,high@cerbero,high@chimera,high@gerione,high@nemeo,high@ortro,high@sfinge,high@tartaro

# Load modules
module load python3/3.8.5
module load tmalign/20190822

# Define path to input lists directory
in_dir=$1
# Define path to output directory
out_dir=$2
# Define path to target structure database
t_db=$3
# Define path to query structure database
q_db=$4
# Define path to list of lists
in_lists="${in_dir}/lists.dat"
# Define path to input list
curr_list=$(sed "${SGE_TASK_ID}q;d" "${in_lists}")
# Loop through each input file
while read in_path; do
  # Debug
  echo "Current input file is ${in_path}"
  # Run Python script (substitute input with output)
  python3 ./bin/graph_tmscore.py -i "${in_path}" -o "${out_dir}" -t "${t_db}" -q "${q_db}"
# Read current line
done < "${curr_list}"
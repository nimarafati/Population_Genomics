#!/usr/bash
input=$1
output_muscle="$1.muscle"
output_tree="$1.tree"

echo "module load muscle
module load FastTree
muscle -in $input -out $output_muscle
FastTree < $output_muscle >$output_tree" >run_Muscle.sh
bash run_Muscle.sh

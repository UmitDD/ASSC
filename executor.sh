#!/usr/bin/env bash
set -euo pipefail

echo "==========Now You will provide information that generate INCAR files for Spin-Spiral Calculations======="
read -p "Enter Number of magnetic basis interaction in unit-cell: " basis
read -p "Enter file name include propagation vectors: " pro_vectors
read -p "Enter lattice constant of unit-cell: " lat_cont
read -p "Enter name of Bravais lattice: " system
read -p "Enter magnetic moments list (e.g [[interaction 1->[m1x,m1y,m1z],[m2x,m2y,m2z]], interaction 2->[[m'1x,m'1y,m'1z],[m'2x,m'2y,m'2z]]]): " mag_moment
read -p "Enter constrained list (e.g magnetic moment list format): " const_moment
read -p "Enter fractional atomic positions (e.g [[p1x,p1y,p1z]]): " positions
read -p "Enter atomic types (e.g atom1->0, atom2->1): " types
read -p "Enter unfolding mesh of IBZ points (e.g [10,10,10]): " mesh
read -p "Enter Scaling factor of energy, mev or mRy: " scaling_factor
read -p "Enter file name of output INCAR (e.g incar_ss): " file_name
read -p "Enter number repeat for Job submission:" repeat

source /home/users/daglum1/p_python/a_library/myenv/bin/activate
python3 Input_generator.py "$basis" "$pro_vectors" "$lat_cont" "$system" "$mag_moment" "$const_moment" "$file_name" || { echo "Sorry Input generation failed" >&2; exit 1; }

job_ids=()

# ==== All jobs go to the queue in parallel ==== #
for ((index=1; index<=basis; index++)); do
  for ((r=0; r<repeat; r++)); do
    jid=$(sbatch --parsable "job_SS${index}.sh") \
      || { echo "Failed to submit job_SS${index}" >&2; exit 1; }
    job_ids+=("$jid")
    echo "Submitted job_SS${index} (repeat $((r+1))/$repeat) -> jobid $jid"
  done
done

# ==== Block until all submitted jobs leave the queue in HPC ==== #
if (( ${#job_ids[@]} > 0 )); then
  id_list=$(IFS=,; echo "${job_ids[*]}")
  echo "Waiting for ${#job_ids[@]} jobs to finish..."
  while squeue -h -j "$id_list" -o "%i" 2>/dev/null | grep -q .; do
    sleep 30
  done
  echo "All jobs finished."
fi

# ==== Merge All Spin Spiral Energies ==== #
for ((index=1; index<=basis; index++)); do
  chmod +x "cat${index}.sh" || { echo "cat${index}.sh is not exist" >&2; exit 1; }
  bash "cat${index}.sh" || { echo "cat${index}.sh failed" >&2; exit 1; }
done

source /home/users/daglum1/p_python/a_library/myenv/bin/activate
python3 Magnetic_script.py "$basis" "$lat_cont" "$system" "$positions" "$types" "$mesh" "$scaling_factor" || { echo "Sorry Magnetic_script failed" >&2;}
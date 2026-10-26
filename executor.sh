#!/usr/bin/env bash
set -euo pipefail

echo "==========Now you will provide information that generates INCAR files for Spin-Spiral Calculations======="
read -p "Enter the number of magnetic basis interactions in the unit cell: " basis
read -p "Enter position index of magnetic atom in unit-cell for each basis interaction: " mag_basis
read -p "Enter file name including propagation vectors: " pro_vectors
read -p "Enter lattice constant of unit-cell: " lat_cont
read -p "Enter name of Bravais lattice: " system
read -p "Enter atomic types: " types
read -p "Enter unfolding mesh of IBZ points: " mesh
read -p "Enter scaling factor of energy, meV or mRy: " scaling_factor
read -p "Enter file name of output INCAR (e.g incar_ss): " file_name
echo "==========Enter these two attributes if there is more than one magnetic basis in the system. Otherwise, you might give 0======="
read -p "Enter SS energy as a function of 90 degrees magnetic moment phase between magnetic basis: " enSS90
read -p "Enter positional differences of two magnetic basis: " delta
read -p "Enter number of repeats for Job submission: " repeat

source /home/users/daglum1/p_python/a_library/myenv/bin/activate
python3 Input_generator.py "$basis" "$pro_vectors" "$lat_cont" "$system" "$file_name" || { echo "Sorry Input generation failed" >&2; exit 1; }

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
python3 Exchange_calculator.py "$basis" "$system" "$delta" "$types" "$mesh" "$scaling_factor" "$mag_basis" "$enSS90" || { echo "Sorry Magnetic_script failed" >&2;}
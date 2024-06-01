#!/bin/bash

# Define the directory path
temperature_dir="../_INPUT/temp_values.in"
execution_dir="../_INPUT/input_orig.dat"

# Clear the content of temperature_dir and create an empty file
> "$temperature_dir"

# clean previous sim data and reexecute
"cd" "../_SOURCE"
"make" "remove"
"make" "clean"
"make"
"cd" "../_INPUT"


# Define the cycling variable
start=3.0
end=0.2
increment=-0.02 # Adjust increment value as needed

# Generate cycling variable values
while (( $(bc <<< "$end <= $start") )); do
   echo "$(printf "%.3f" "$start")" >> "$temperature_dir"
   start=$(bc <<< "$start + $increment")
done

# Check if the temperature values file exists
if [ ! -f "$temperature_dir" ]; then
   echo "File $temperature_dir not found."
   exit 1
fi

# Loop over the two different methods, Gibbs and Metropolis
for method in 2 3; do

	# Modify the value of H in the execution file
	awk -v m="$method" 'NR==1{$2=m}1' "$execution_dir" > "${execution_dir}.tmp" && mv "${execution_dir}.tmp" "$execution_dir" || { echo "Failed to modify file."; exit 1; }

	# Loop over the two different values of H
	for H_value in 0 0.02; do

		# Counter variable to track the number of iterations
		iteration=0

		# Modify the value of H in the execution file
		awk -v H="$H_value" 'NR==1{$4=H}1' "$execution_dir" > "${execution_dir}.tmp" && mv "${execution_dir}.tmp" "$execution_dir" || { echo "Failed to modify file."; exit 1; }

		# Loop over each line in the temperature file
		while IFS= read -r execution_temperature || [ -n "$execution_temperature" ]; do

			# Modify the equilibration file directly
			awk -v temp="$execution_temperature" 'NR==2{$2=temp}1' "$execution_dir" > "${execution_dir}.tmp" && mv "${execution_dir}.tmp" "$execution_dir" || { echo "Failed to modify file."; exit 1; }

			# Check if it's the first or second iteration
			if [ "$iteration" -eq 0 ]; then
				awk -v restart="0" 'NR==3{$2=restart}1' "$execution_dir" > "${execution_dir}.tmp" && mv "${execution_dir}.tmp" "$execution_dir" || { echo "Failed to modify file."; exit 1; }
			fi
			if [ "$iteration" -eq 1 ]; then
				awk -v restart="1" 'NR==3{$2=restart}1' "$execution_dir" > "${execution_dir}.tmp" && mv "${execution_dir}.tmp" "$execution_dir" || { echo "Failed to modify file."; exit 1; }
			fi

			# Execute simulator.exe with the modified file
			"../_SOURCE/simulator.exe" "$execution_dir" || { echo "Failed to execute simulator."; exit 1; }

			# Increment the iteration counter
			((iteration++))

		done < "$temperature_dir"
		
	done

done
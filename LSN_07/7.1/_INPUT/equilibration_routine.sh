#!/bin/bash

# Check if the user provided an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <phase> = {Solid,Liquid,Gas}"
    exit 1
fi

PHASE=$1  # Assign the first argument

# Execute the Python script to generate temp_values.in
python3 ./generate_temp_values.py "$PHASE"

# Check if the Python script executed successfully
if [ $? -ne 0 ]; then
    echo "Error in Python script. Please ensure you provided a valid phase name."
    exit 1
fi

# clear previous equilibration
"cd" "../_SOURCE"
#"make" "$PHASE"
"make" "clean"
"make"
"cd" "../_INPUT"

# Construct the file paths
dir="../$PHASE/INPUT"
equilib_file="$dir/input_equilibration.in"
temp_values_file="$dir/temp_values.in"

# Check if the equilibration file exists
if [ ! -f "$equilib_file" ]; then
    echo "File $equilib_file not found."
    exit 1
fi

# Check if the temperature values file exists
if [ ! -f "$temp_values_file" ]; then
    echo "File $temp_values_file not found."
    exit 1
fi

# Loop over each line in the temperature values file
while IFS= read -r equilib_temperature || [ -n "$equilib_temperature" ]; do
    # Modify the equilibration file directly using sed
    awk -v temp="$equilib_temperature" 'NR==2{$2="\t\t\t\t\t\t\t"temp}1' "$equilib_file" > "${equilib_file}.tmp" && mv "${equilib_file}.tmp" "$equilib_file" || { echo "Failed to modify file."; exit 1; }
    # Execute simulator.exe with the modified file
    "../_SOURCE/simulator.exe" "$equilib_file" || { echo "Failed to execute simulator."; exit 1; }
done < "$temp_values_file"

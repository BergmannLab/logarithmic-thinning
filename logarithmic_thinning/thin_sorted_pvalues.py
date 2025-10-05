#!/bin/python3
#SBATCH --account=<account>
#SBATCH --job-name=thin
#SBATCH --error=slurm-%j.err
#SBATCH --output=slurm-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 100G
#SBATCH --time 1:00:00
#SBATCH --partition urblauna

import time
import subprocess
import argparse
import os
import math
import tempfile

# Function to determine the threshold for logarithmic thinning
def determine_threshold(factor):
    return round(factor / (factor - 1.0))

# Logarithmic thinning function to determine row indices to keep
def logarithmic_thinning(N, factor):
    # Print initial number of p-values
    print("Number of p-values before thinning:", N)

    # Initialize the list of indices
    indices = []

    # Start from the initial value N
    index = N
    indices.append(index)

    # Calculate threshold once
    threshold = determine_threshold(factor)
    print("Number of unthinned top p-values:", threshold)

    # Perform logarithmic thinning until the threshold
    while index > threshold:
        index = math.floor(index / factor)
        indices.append(index)

    # Add the remaining indices down to 1
    indices.extend(range(math.floor(threshold) - 1, 0, -1))

    print("Total number of p-values after thinning:", len(indices))

    return sorted(indices)

# Python version to filter specific rows based on logarithmic thinning indices
def filter_specific_rows_python(input_file, output_file, rows_to_keep):
    rows_set = set(rows_to_keep)  # Convert to a set for faster lookups
    with open(input_file, 'r', buffering=16 * 1024 * 1024) as infile, \
         open(output_file, 'w', buffering=16 * 1024 * 1024) as outfile:

        for i, line in enumerate(infile, start=1):
            if i in rows_set:
                line = line.rstrip() + f',{i}\n'
                outfile.write(line)

# Optimized AWK version to filter specific rows using temporary file for large numbers of rows
def filter_specific_rows_awk(input_file, output_file, rows_to_keep):
    # Write row numbers to a temporary file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        for row in rows_to_keep:
            temp_file.write(f'{row}\n')
        temp_filename = temp_file.name

    # Use AWK with a temporary file to match row numbers for performance
    awk_command = f"awk 'FNR==NR {{ lines[$1]; next }} (FNR in lines)' {temp_filename} {input_file} > {output_file}"
    subprocess.run(awk_command, shell=True, check=True)

    # Clean up the temporary file
    os.remove(temp_filename)

# Function to time the execution of a given function
def time_execution(function, *args):
    start_time = time.time()
    function(*args)
    end_time = time.time()
    elapsed_time = end_time - start_time
    return elapsed_time

# Main function to handle command-line arguments and execute the selected method
def main():
    parser = argparse.ArgumentParser(description="Thinning specific rows from a large file.")
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('thinning_factor', type=float, help='Thinning factor.')
    parser.add_argument('method', type=str, choices=['python', 'awk'], help='Method to use: python for Python, awk for AWK. (In this implementation, Python seems to be significantly faster)')

    args = parser.parse_args()

    input_file = args.input_file
    thinning_factor = args.thinning_factor
    method = args.method

    # Generate output filename based on method in the same path as input file
    output_file = os.path.join(os.path.dirname(input_file), f"thinned_python_{os.path.basename(input_file)}" if method == 'python' else f"thinned_awk_{os.path.basename(input_file)}")
    
    # Get the number of rows in the input file
    #with open(input_file, 'r') as f:
    #    N = sum(1 for _ in f)
    N = 15974975488

    # Calculate which rows to keep using logarithmic thinning
    rows_to_keep = logarithmic_thinning(N, thinning_factor)

    # Run the selected method and time it
    if method == 'python':
        # Run the Python version
        print("Running Python version...")
        elapsed_time = time_execution(filter_specific_rows_python, input_file, output_file, rows_to_keep)
        print(f"Python version took {elapsed_time:.2f} seconds.")

    elif method == 'awk':
        # Run the AWK version
        print("Running AWK version...")
        elapsed_time = time_execution(filter_specific_rows_awk, input_file, output_file, rows_to_keep)
        print(f"AWK version took {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()


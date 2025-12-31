#!/bin/bash

# Check for correct usage
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 -i <slha_file> -o <output_file> -n <number_of_events> -s <random_seed>"
    exit 1
fi

# Parse command line arguments
while getopts ":i:o:n:s:" opt; do
  case $opt in
    i)
      slha_file="$OPTARG"
      ;;
    o)
      output_file="$OPTARG"
      ;;
    n)
      num_events="$OPTARG"
      ;;
    s)
      seed="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Ensure the generator binary exists
if [ ! -f "./generator" ]; then
    echo "Error: generator binary not found. Please compile it first."
    exit 1
fi

# Set the dynamic linker path
export DYLD_LIBRARY_PATH=/Users/buding/Workshop/HEPTools/pythia8/lib:/Users/buding/Workshop/HEPTools/HepMC3/build/lib:$DYLD_LIBRARY_PATH

# Run the generator with the new parameters
echo "Input SLHA file   -> $slha_file"
echo "Write event file  -> $output_file"
echo "Number of events  -> $num_events"
echo "Random seed       -> $seed"
echo "Running the generator..."
./generator -i "$slha_file" -o "$output_file" -n "$num_events" -s "$seed"

# Check if the generator ran successfully
if [ $? -eq 0 ]; then
    echo "Generation completed successfully."
else
    echo "Generation failed"
    exit 1
fi
#!/bin/bash
set -e  # exit immediately if any command fails, so a bad mv doesn't silently corrupt the structure

# -- Pass 1 ---------------
# Move every flat SRR*.bw file into a first-level subfolder named SRR + first 3 digits.
# Example: SRR9295900_pshifted_forward.bw  →  SRR929/SRR9295900_pshifted_forward.bw
for file in ./SRR*; do
    [ -f "$file" ] || continue                            # skip if it's already a directory (from a previous run)
    base=$(basename "$file")                              # strip the ./ prefix  →  SRR9295900_pshifted_forward.bw
    num=$(echo "$base" | grep -oP 'SRR\K\d{3}')          # \K resets the match start, so we capture only the 3 digits: 929
    mkdir -p "./SRR$num"                                  # create SRR929/ if it doesn't exist yet (-p silences the error if it does)
    mv "$file" "./SRR$num/"                               # move the file one level down
done

# -- Pass 2 ---------------
# Descend into the SRR929/ folders and push files one level deeper using the first 5 digits.
# Example: SRR929/SRR9295900_…  →  SRR929/SRR92959/SRR9295900_…
for file in ./SRR*/SRR*; do
    [ -f "$file" ] || continue                            # skip subdirectories themselves
    base=$(basename "$file")                              # filename only, no path
    num=$(echo "$base" | grep -oP 'SRR\K\d{5}')          # first 5 digits after SRR: 92959
    parent=$(dirname "$file")                             # the SRR929/ folder this file currently lives in
    mkdir -p "$parent/SRR$num"                            # create SRR929/SRR92959/ if needed
    mv "$file" "$parent/SRR$num/"                         # move the file into the deeper folder
done

echo "Done!" 

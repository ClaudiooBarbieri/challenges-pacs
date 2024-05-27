#!/bin/bash

# directory of executable
DIRECTORY="../src"

# compile and link parallel version
make parallel

# Define the values for X
values=(1 2 4)

# Loop over each value and run mpiexec
for X in "${values[@]}"
do
  echo "Running mpiexec with $X"
  mpiexec -np $X "$DIRECTORY/main_parallel"
done

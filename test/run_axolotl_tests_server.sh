#!/bin/bash

# Change to the correct directory
cd /data/axela/BuildChaste/

# Make the test
make -j$1 $2

start_run=$3;
end_run=$4;

DRAG=1.0;
ADHESION=10.0;
MIGRATION=10.0;

echo "Beginning runs from ${start_run} to ${end_run}."
for ((i=start_run; i<=end_run; i++))
do
# "&" on the end lets the script carry on and not wait until this has finished.
projects/$3/test/$2 -seed ${i} &
done

echo "Jobs submitted"
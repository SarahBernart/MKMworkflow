#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

python extracted.py
echo "Finished running extracted.py."

python ped.py
echo "Finished running ped.py."

rm -r temp

python temp.py
echo "Finished running temp.py."

python plot.py

#python plot_temp.py

python input.py
echo "Finished running input.py."

python input_rates.py
echo "Finished running input_rates.py."

python input_combined.py
echo "Finished running input_combined.py."

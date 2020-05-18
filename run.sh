#!/bin/bash
set -e

# Build binaries
bash build.sh

# prep parallel command
parallel --citation

# Run Computational Experiments
bash RunFigure3Simulations.sh
bash RunFigure4Simulations.sh
bash RunFigure5Simulations.sh
bash RunFigure6Simulations.sh
bash RunTable1Simulations.sh

# Run Visualization Scripts
bash Figure3.py
bash Figure4.py
bash Figure5.py
bash Figure6.py
bash Table1.py

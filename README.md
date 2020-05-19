# Single-node second-order boundary schemes for the lattice Boltzmann method

This is a Reproduction Package as described in the manuscript "Three Empirical
Principles for Computational Reproducibility and their Implementation:
*The Reproduction Package*” by M. S. Krafczyk, A. Shi, A. Bhaskar,
D. Marinov, & V. Stodden

["Single-node second-order boundary schemes for the lattice Boltzmann method"](https://dx.doi.org/10.1016/j.jcp.2016.10.049)
describes a scheme for obtaining second order
accuracy on boundaries in Lattice Boltzmann simulations. This method especially
treats curved boundaries better than other methods. The authors demonstrate
this performance through a series of three toy problems for which the analytic
solutions is known. These examples are the Poiseuille flow, the Taylor-Green
vortex flow, and the Taylor-Couette flow. Additionally they demonstrate
stability of their technique to converge in a reasonable amount of time for a
wide variety of input parameters with the Poiseuille flow.

## Build Instructions

This is a simple package to build and run containing only simple source code files.

### Requirements

The following are required to build and run this software:

* A c++ compiler such as the gnu compiler `g++`.
* Gnu make
* Gnu parallel

### Building with docker

    docker build -t ${DOCKER_IMAGE_NAME} .

### Running with docker

    docker run -it --rm -v $(pwd):/Scratch ${DOCKER_IMAGE_NAME}

#### Run everything script

We provide the script `run.sh` which will build the software, run the computational experiments and run the visualization output assuming that all necessary dependencies are available or installed.

Before executing this script, please read the (Computational Effort)[computational_effort.md] document to be sure that you have the necessary resources available to execute this script.

More details about what is run in this script follow.

#### Building with linux

The build scripts assume you are using the gnu compiler suite and have `g++` available. Changing to a different compiler is as simple as editing the `Makefile` to change instances of `g++` to your favorite compiler.

Simply execute the `build.sh` script to perform this build process if you don't want to think about it.

#### Scripts for Computational Experiments.

The computational experiments backing each of the Figures and Tables each have a catch-all script which executes the necessary experiments for them. These are:

* `RunFigure3Simulations.sh`
* `RunFigure4Simulations.sh`
* `RunFigure5Simulations.sh`
* `RunFigure6Simulations.sh`
* `RunTable1Simulations.sh`

Each of these scripts utilizes gnu parallel run experiments simultaneously on multiple cores. Please adjust the number of cores available within the scripts to improve performance on your system.

Each of these scripts produces a set of outputs. Generally speaking each figure and table draws from multiple output files. These output files have the following form:

    FigureX_Y.out
    Table1_Y.out

Output produced by our machine can be found in the directory `expected_output`

#### Scripts for Visualizations

Once the computational experiments have been performed, the produced data can be analyzed to produce the Figures and plots. The following scripts can be run to produce these:

* `Figure3.py`
* `Figure4.py`
* `Figure5.py`
* `Figure6.py`
* `Table1.py`

The output of these scripts can be found in the main directory with endings like `.png` in the case of Figures and printed to the terminal in the case of the Table.

## Reproduction Notes

We kept track of our progress and issues inside `notes.txt`. We also have a jupyter notebook showing this progress over time `ReproducibilityPlot.ipynb`.

## Acknowledgement of the Authors

We want to acknowledge the authors Weifeng Zhao and Wen-An Yong for their fine work on this computational experiment. We succeeded to reproduce most aspects of their paper which would not have been possible if their work was not already well done.

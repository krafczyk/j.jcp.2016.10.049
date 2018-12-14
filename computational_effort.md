# Computational Effort

To produce the five assets from this article, we provide 5 scripts to automate
the generation of computational results and and 4 scripts to produce plots as
visualizations. We provide here, estimates of how long these experiments should
take to complete with our hardware, An Intel(R) Core(TM) i7-6900K CPU @ 3.20GHz.

Each experiment requires a small amount of RAM < 200MB.

## Computational Experiments

We provide driver scripts to execute each experiment for each Figure with
parameters for each of the experiments stored in a file.

* Figure 3, convex and nonconvex: <24 hours each @ 8 cores
* Figure 4, convex and nonconvex: <24 hours each @ 8 cores
* Figure 5, convex and nonconvex: <24 hours each @ 8 cores
* Figure 6, convex: < 24 hours @ 8 cores

For the single table, we provide a driver script for an individual component of
the table as well as another driver script to coordinate the experiments for
all portions of the table.

* Table 1: ~48 hours @ 8 cores

## Visualizations

All Visualizations will take < 5 minutes each @ 1 core. A visualization script
written in python using the matplotlib library is written for each asset.

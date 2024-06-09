# L-SRTDE algorithm for GNBG 2024 benchmark

C++ implementation of the L-SRTDE algrithm (Linear population size reduction Success RaTe-based adaptive Differential Evolution)

This version is for the GECCO 2024 Numerical Global Optimization Competition on GNBG-generated Test Suite (https://competition-hub.github.io/GNBG-Competition/) and in particular implemented in C++ (https://github.com/VladimirStanovov/GNBG-Instance-C)

The algorithm code is in "L-SRTDE.cpp" file.

# Compilation and usage

In order to run the experiments, first run convert.py to convert GNBG instance from .mat to .txt with specific format, available to read by C++ code.

Compilation is simple using gcc/g++:

g++ -std=c++11 -O3 L-SRTDE.cpp -o L-SRTDE.exe

or depending on hardware

g++ -std=c++11 -O3 -march=corei7-avx L-SRTDE.cpp -o L-SRTDE.exe

Please note that the compilation requires support of C++11 standard.

This will create L-SRTDE executable, available for running.

Data will be written to "L-NTADE_GNBG_F#_D#.txt", where F and DIM are the function number and problem dimention.

process_results.py will summarize all results and give averages and standard deviation.

make_graphs.py will generate convergence graphs for all 24 functions and 31 runs.

Results folder contains the results recevived by running this code.

# GNBG 2024 competition info

Competition page: https://competition-hub.github.io/GNBG-Competition/

Reference:

D. Yazdani, M. N. Omidvar, D. Yazdani, K. Deb, and A. H. Gandomi, "GNBG: A Generalized
  and Configurable Benchmark Generator for Continuous Numerical Optimization," arXiv prepring	arXiv:2312.07083, 2023.
  
A. H. Gandomi, D. Yazdani, M. N. Omidvar, and K. Deb, "GNBG-Generated Test Suite for Box-Constrained Numerical Global
  Optimization," arXiv preprint arXiv:2312.07034, 2023.

# Example code for Projected Rayleigh Quotient Iteration

This is the code that was used to produce the results and figures for the _Projected Rayleigh Quotient Iteration_ (PRQI) paper. All examples were only tested on Linux.

## Instructions for Example 1
This example requires a C++ 11 compiler and Python 3. The C++ program generates the results, the Python script then visualises them. 

To compile the C++ executable the [Blaze library](https://bitbucket.org/blaze-lib/blaze) is required which itself relies on LAPACK. See the Blaze documentation for details on how to install both. To compile the program (e.g., using the Clang compiler) inside the `example1` directory run
```bash
$ clang++ -O3 -ogen_results -llapack main.cc
```
and then execute
```bash
$ ./gen_results 2000
```
to generate the results using 2000 runs (the more runs, the more accurate the final figures). To generate the figures, run (this assumes that both `numpy` and `matplotlib` are available)
```
$ mkdir figures
$ python main.py
```
This generates the figures in a subdirectory `figures`. Note that the generated figures differ slightly from those in the paper. The published figures were post-processed to increase readability when printed in grayscale. 

## Instructions for Example 2
This example only requires MATLAB. It was tested in MATLAB R2023b. To generate the results just run the script `generate_results.m` in MATLAB.

To generate the figure directly from the command line, run
```bash
$ matlab -nodisplay -nosplash -nodesktop -r "run('generate_results.m');exit;"
```

## Instructions for Example 3
This example only requires MATLAB. It was tested in MATLAB R2023b. To generate the results just run the script `generate_results.m` in MATLAB. Note that this only creates Figure 5 of the paper and the results that were used to create Table 1. The other figures are simple plots of MATLAB vectors and are therefore not included here. 

To generate the figure directly from the command line, run
```bash
$ matlab -nodisplay -nosplash -nodesktop -r "run('generate_results.m');exit;"
```

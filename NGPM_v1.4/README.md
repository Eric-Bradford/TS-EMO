# NGPM -- A NSGA-II Program in Matlab v1.4 by Song Lin

This program is an implementation of nondominated sorting genetic algorithm II (NSGA-II) proposed by K. Deb. Capabilities:
1. R-NSGA-II: Reference-point-based NSGA-II. 
2. Coding: real, integer. 
3. GA operator: Intermediate crossover, Gaussian mutation. 
4. Constraint handling. 
5. Parallel computation of objective function evaluation. 
6. Population plot in window. 
I write this program because Aravind Seshadri’ program (File ID#10429) could not satisfy my request. I need constraint handling, integer coding to solve a finite element optimization problem. The finite element solution is very time-expensive, thus the parallel computation of objective evaluation is implemented in the code.

See https://se.mathworks.com/matlabcentral/fileexchange/31166-ngpm-a-nsga-ii-program-in-matlab-v1-4 for downloading the original unmodified NSGA-II algorithm and http://www.codelooker.com/dfilec/6987NGPMv11.4/NGPMmanualv1.4.pdf for more information

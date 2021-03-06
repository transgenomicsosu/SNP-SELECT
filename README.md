# SNP-SELECT

Shuzhen Sun, Zhuqi Miao, Blaise Ratcliffe, Polly Campbell, Bret Pasch, Yousry A. El-Kassaby, Balabhaskar Balasundaram, Charles Chen. 2018. "Graph domination for SNP variable selection".

SNP variable selection tool

This software is developed to select significan SNP variables without phenotype data.

This program is implemented for numeric genotype data that is coded as 0, 1, and 2. 

"kdom.cpp" is the source file;

"data_test.csv" is the data file used to test the program, and the results are included in the output branch;

"simulation.zip" are the data files corresponding to the simulation study;

"mice.zip" is the data file corresponding to the clustering analysis;

"Douglas_fir.zip" is the data file corresponding to the pedigree recovery analysis;

Users can adjust the paramaters, including input file (INPUT), threshold value (THRE) and k value (K), based on their need.

Input and output pathway need to be adjusted as well.

To successfully run this C++ program, it requires Gurobi Optimizer being installed. Please refer to http://www.gurobi.com/documentation/ for the details about installing Gurobi Optimizer.

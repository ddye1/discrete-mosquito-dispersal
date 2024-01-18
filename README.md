# discrete-mosquito-dispersal
Code for "Efficacy of Wolbachia-based mosquito control: Predictions of a spatially discrete mathematical model" (Dye and Cain, 2024, PLOS ONE). DOI: 10.1371/journal.pone.0297964

This directory contains code for evaluating and plotting models presented in the referenced paper. The C Code directory contains 
three programs written in C++ and some auxiliary files that can be used (optionally) to plot output using the 
freely available gnuplot software. The MATLAB Code directory contains two MATLAB scripts that solve the two model equations and
produce additional figures from the referenced paper. Both the C++ code and the MATLAB code produce identical results.

## C Code

1) hu-model.C:  Solves the model equations of Hu et al., Bulletin of Mathematical Biology 
(2021) 83:58, adapted to include two diffusively-coupled habitats (A and B).  This program 
calls the gnuplot script huplot, stored in this directory, to generate plots of the plain text 
output files.  If you prefer to use something other than gnuplot to generate figures, simply 
comment out the line at the end of hu-model.C in which gnuplot is called.  To compile this 
program using the g++ compiler, use the command

               g++ -o hu-model hu-model.C
   
2) qu-model.C:  Needs to be stored in the same directory as parameters.h.  Solves the model 
equations of Qu et al. (2018) Siam J Appl Math, pp. 826--852 that accounts for sexes and age 
class of mosquitoes, adapted to include two diffusively coupled habitats (A and B).  This 
program calls the gnuplot script quplot, stored in this directory, to generate plots of the 
plain text output files.  If you prefer to use something other than gnuplot to generate 
figures, simply comment out the line at the end of qu-model.C in which gnuplot is called.  To 
compile this program using the g++ compiler, use the command

               g++ -o qu-model qu-model.C

3) bifurcation-diagram.C:  Needs to be stored in the same directory as parameters.h.  Solves 
the model equations of Qu et al. (2018) Siam J Appl Math, pp. 826--852 that accounts for sexes 
and age class of mosquitoes, adapted to include two diffusively coupled habitats (A and B).  
This program computes the steady-state proportion of Wolbachia-infected mosquitoes in each 
habitat for various choices of the migration parameter m.  The output can then be used to 
generate a bifurcation diagram showing infected populations versus m. This program calls the 
gnuplot script bifplot, stored in this directory, to generate plots of the plain text output 
files.  If you prefer to use something other than gnuplot to generate figures, simply comment 
out the line at the end of bifurcation-diagram.C in which gnuplot is called.  To compile this 
program using the g++ compiler, use the command

               g++ -o bifurcation-diagram bifurcation-diagram.C

## MATLAB Code

1) model_1.m: Solves the model equations of Hu et al., Bulletin of Mathematical Biology 
(2021) 83:58, adapted to include either two diffusively-coupled habitats (A and B) or N
diffusively-coupled habitats. A heatmap representation of lognormal and exponential dispersal
kernels is also given for two habitats with equal initial populations of mosquitoes. A function
is provided to create a video of Wolbachia-infected mosquitoes invading an N-habitat model.

2) model_2.m: Solves the model equations of Qu et al. (2018) Siam J Appl Math, pp. 826--852 that
accounts for sexes and age class of mosquitoes, adapted to include two diffusively coupled habitats
(A and B). Functions for creating the plots presented in the paper are provided.

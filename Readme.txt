Title: Non-rigid Puzzles
Authors: Or Litany, Emanuele Rodola, Alex Bronstein, Michael Bronstein, Daniel Cremers
OS: Windows 10

Instructions:
- The implementation was done with Matlab 2015b 64bit
- To reproduce the results simply use the script in the main directory: reproduce_paper_figures.m
- The only non self-contained dependency is the manopt toolbox. It should be downloaded from: http://www.manopt.org/
- After you download manopt, make sure to change the path to where you saved it. 

Remark:
In case you get an error related to one of the files: mumford_shah_wrapper.mexw64 or norm_21_gradient.mexw64 please recompile by typing in Matlab:
  mex mumford_shah_wrapper.cpp COMPFLAGS="/openmp $COMPFLAGS"
  mex norm_21_gradient.cpp COMPFLAGS="/openmp $COMPFLAGS"

License: 
You are free to use and modify this code for academic purposes. If you do, please cite the paper using:

@article{litany16,
  author = {O. Litany and E. Rodola and A. M. Bronstein and M. M. Bronstein and D. Cremers},
  title = {Non-Rigid Puzzles},
  journal = {Computer Graphics Forum},
  volume = {35},
  number = {5},
  year = {2016},
  publisher = {Wiley},
  topic = {Shape Analysis, Segmentation},
  titleurl = {litany16.pdf},
}
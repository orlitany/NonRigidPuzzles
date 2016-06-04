% This script reproduces the experiment results shown in the paper "Non-rigid Puzzles"
% The required functions are availiable at https://github.com/orlitany/NonRigidPuzzles
% The only non self-contained dependency is the manopt toolbox. It should be downloaded from: http://www.manopt.org/
% After you download manopt, make sure to change the path to where you saved it. 
% To run an experiment, simply uncomment the relevant subsection. 
% Results will be saved automatically under a new folder. 
% If you encounter problems, feel free to contact us at: orlitany <at> post <dot> tau <dot> ac <dot> il

clear all;close all;
%% dependencies
addpath(genpath('./'));
addpath(genpath('./../../../manopt/')) %Replace the the correct path to manopt folder

%% Figure 1 
load('./Data\Figure1_data.mat')
main_params.output_folder = ['./Results\Figure1\' datestr(now, 'ddmmHHMM') '\']; 

%% Figure 6 (dog parts) [~45 mins]
% render_figure = 'Figure6'
% load(['./Data\' render_figure '_data.mat'])
% main_params.output_folder = ['./Results\' render_figure '\']; 

%% Figure 7 (real scan) [~8.5 mins]
% render_figure = 'Figure7';
% load(['./Data\' render_figure '_data.mat'])
% main_params.output_folder = ['./Results\' render_figure '\']; 
%% Figure 8 [~ 1 hour]
% load('./Data\Figure8_data.mat')
% main_params.output_folder = ['./Results\Figure8\' datestr(now, 'ddmmHHMM') '\']; 

%% Figure 9 (overlapping parts) [~ 35 min]
% render_figure = 'Figure9'
% load(['./Data\' render_figure '_data.mat'])
% main_params.output_folder = ['./Results\' render_figure '\']; 

%% Figure 10 (wolf with extra part) [~48 mins]
% render_figure = 'Figure10'
% load(['./Data\' render_figure '_data.mat'])
% main_params.output_folder = ['./Results\' render_figure '\']; 

%% run 
mkdir(main_params.output_folder);
tic
[C,u,v] = solveNonRigidPuzzle(model,parts,corr_functions,main_params,algo_params);
toc

%% render
switch(render_figure)
    case 'Figure6'        
        render_Figure6
    case 'Figure10'        
        render_Figure10
    case 'Figure9'        
        render_Figure9
end


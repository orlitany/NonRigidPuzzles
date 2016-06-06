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
% load('./Data\Figure1a_data.mat')
% main_params.output_folder = ['./Results\Figure1\']; 
% main_params.is_clutter = 1;
% main_params.is_missing_part = 0;
% algo_params.num_iter.C_step = 10
% algo_params.lambda_area_parts = 1000
% algo_params.parts_area_thresh = 0.85;
% [C,u,v] = solveNonRigidPuzzle(model,parts,corr_functions,main_params,algo_params);
% 
% irrelevant_parts_ind = [];
% irrelevant_parts = []; jj=1;
% for prt=1:main_params.num_parts    
%     [~,~,~,arr] = extract_eigen_functions_new(parts{prt}.shape,1);
%     arr = full(diag(arr));
%     if sum(arr.*v{prt})/sum(arr) < algo_params.parts_area_thresh,
%         irrelevant_parts_ind = [irrelevant_parts_ind;prt];
%         irrelevant_parts{jj} = parts{prt};
%         jj = jj+1;
%     end
% end
% load('./Data\Figure1a_data.mat')
% main_params.output_folder = ['./Results\Figure1\']; 
% parts(irrelevant_parts_ind) = [];
% corr_functions(irrelevant_parts_ind) = [];
% main_params.num_parts = numel(parts);
%% Figure 6 (dog parts) [~45 mins]
% render_figure = 'Figure6'
% load(['./Data\' render_figure '_data.mat'])
% main_params.output_folder = ['./Results\' render_figure '\']; 

%% Figure 7 (real scan) [~8.5 mins]
render_figure = 'Figure7';
load(['./Data\' render_figure '_data.mat'])
main_params.output_folder = ['./Results\' render_figure '\']; 
%% Figure 8 [~ 1 hour]
% render_figure = 'Figure8';
% load(['./Data\' render_figure '_data.mat'])
% main_params.output_folder = ['./Results\' render_figure '\']; 
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

%% render (Figure 7 is generated during the code - thus no need for a separate rendering)
switch(render_figure)
    case 'Figure1'        
        render_Figure1
    case 'Figure6'        
        render_Figure6
    case 'Figure8'        
        render_Figure8
    case 'Figure9'        
        render_Figure9
    case 'Figure10'        
        render_Figure10    
end
%% Demo of the non-rigid puzzle algorithm.
clear all;
%% Parameters
main_params.input_folder = './Data\';
main_params.input_file = 'Figure7_data';
main_params.output_folder = ['./Results\Figure7\' datestr(now, 'ddmmHHMM') '\']; 
mkdir(main_params.output_folder);
main_params.num_eigen = 50;
main_params.is_clutter = 0;
main_params.is_missing_part = 0;
%%
%-
algo_params.num_iter.C_step = 1;
algo_params.num_iter.ransac = 2;
algo_params.num_iter.manopt_maxiter_C = 5e3;
algo_params.num_iter.manopt_maxiter_u = 1500;
algo_params.num_iter.icp = 0;
algo_params.lambda_area_model = 1;
algo_params.lambda_ones = 1e7;
algo_params.lambda_reg_model = 1e2;
algo_params.mu_1 = 1;
algo_params.mu_2 = 1e3;
algo_params.num_ransac_corr = 50;

%% Input shapes
%- load model, parts and correspondence 
load([main_params.input_folder main_params.input_file]);
main_params.num_subspace = numel(model.shape.X);
main_params.num_parts = numel(parts);

%% Solve
close all;
save([main_params.output_folder 'parameters']);                             
[C,u,v] = solveNonRigidPuzzle(model,parts,corr_functions,main_params,algo_params);

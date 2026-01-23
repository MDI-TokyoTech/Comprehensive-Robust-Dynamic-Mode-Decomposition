% main_cylinder.m
% Script to demonstrate CR-DMD preprocessing on Cylinder data

%% 1. Setup Environment
clc; clear; close all;

% Get the directory of this script to ensure relative paths work
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);

% Add necessary paths
addpath(genpath(fullfile(rootDir, 'src')));
addpath(genpath(fullfile(rootDir, 'utils')));

%% 2. Load Data
dataName = 'cylinder';
[X, h, v] = load_cylinder_data(rootDir, dataName);
[N, M] = size(X);

%% 2.5 Add Synthetic Noise
sigma = 0.05;
outlier_ratio = 0.05;
[X_noisy, X_clean_gt] = add_synthetic_noise(X, sigma, outlier_ratio);

%% 3. Configure CR-DMD
alpha = 0.95;
params.pre.w = 0.9;
params.pre.epsilon = alpha * sigma * sqrt((1-outlier_ratio) * N * M); % Adjust based on noise level
params.pre.eta = 0.5 * alpha * outlier_ratio * N * M;     % Adjust based on outlier sparseness
params.pre.max_iter = 5000; % Set lower for quick test, higher for proper convergence
params.pre.criteria = 1e-4;
params.me.r = 21;            % Number of DMD modes to extract

% Initialize CR-DMD object
crdmd = CR_DMD(params);

% Adjust display settings
crdmd.display_freq = 5; 

%% 4. Run Preprocessing
tic;
[X_clean, S_outliers, Phi, Lambda, b] = crdmd.preprocess(X_noisy, h, v);
elapsedTime = toc;
fprintf('Preprocessing finished in %.2f seconds.\n', elapsedTime);

%% 5. Visualization Preprocessing
visualize_preprocessing(crdmd, X_clean_gt, X_noisy, X_clean, S_outliers, Lambda, h, v, sigma, outlier_ratio);

%% 6. Dimensional Reduction
% Configure and run Dimensional Reduction
crdmd.params.dr.w = 0.9;
crdmd.params.dr.epsilon = alpha * sigma * sqrt((1-outlier_ratio) * N * M); % Adjust based on noise level
crdmd.params.dr.eta = 0.5 * alpha * outlier_ratio * N * M;     % Adjust based on outlier sparseness
crdmd.params.dr.max_iter = 5000; % Set lower for quick test, higher for proper convergence
crdmd.params.dr.criteria = 1e-4;
crdmd.params.dr.mu = 1;

fprintf('Running Dimensional Reduction...\n');
tic;
% Note: dimensional_reduction now accepts h and v
[xi, S_outliers_dr] = crdmd.dimensional_reduction(X_noisy, Phi, Lambda, b, h, v);
elapsedTime = toc;
fprintf('Dimensional Reduction finished in %.2f seconds.\n', elapsedTime);

% Visualize Dimensional Reduction
visualize_dimensional_reduction(X_clean_gt, xi, S_outliers_dr, Phi, Lambda, h, v);


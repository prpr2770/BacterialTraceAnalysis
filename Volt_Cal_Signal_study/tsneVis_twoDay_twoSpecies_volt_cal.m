%{
Script to extract model parameters from volt-ca signals of Ecoli-Styph
dataset, on two different days, and then observe the
clustering/visualization obtained. 

%}
clear all; close all; clc;
model_fit_threshold = 51;
% --------------------------------------------------------
% Extract datasets:

inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\ecoli\transferFunctionStudy\inputVOLT';
file_nm = 'volt_input_model_parameters.mat';
load([inDIR filesep file_nm]);

ecoli_volt_params = model_est_pvecs(2:end-1,:);
ecoli_members = ones(1,size(ecoli_volt_params,2));
ecoli_fit = model_est_fit;

%

numEcoli_total_1 = size(model_est_pvecs,2);
good_ecoli_indx_1 = find(ecoli_fit > model_fit_threshold);
numEcoli_good_1 = length(good_ecoli_indx_1);
GOOD_ECOLI_1 = ecoli_volt_params(:,good_ecoli_indx_1)';
% --------------------------------------------------------
inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\ecoli\transferFunctionStudy\inputVOLT';
file_nm = 'volt_input_model_parameters.mat';
load([inDIR filesep file_nm]);

ecoli_volt_params = model_est_pvecs(2:end-1,:);
ecoli_members = ones(1,size(ecoli_volt_params,2));
ecoli_fit = model_est_fit;

%
numEcoli_total_2 = size(model_est_pvecs,2);
good_ecoli_indx_2 = find(ecoli_fit > model_fit_threshold);
numEcoli_good_2 = length(good_ecoli_indx_2);
GOOD_ECOLI_2 = ecoli_volt_params(:,good_ecoli_indx_2)';

% --------------------------------------------------------
% --------------------------------------------------------
inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\styph\transferFunctionStudy\inputVOLT';
file_nm = 'volt_input_model_parameters.mat';
load([inDIR filesep file_nm]);

styph_volt_params = model_est_pvecs(2:end-1,:);
styph_members = 2*ones(1,size(styph_volt_params,2));
styph_fit = model_est_fit;

numStyph_total_1 = size(model_est_pvecs,2);
good_styph_indx_1 = find(styph_fit > model_fit_threshold);
numStyph_good_1 = length(good_styph_indx_1);
GOOD_STYPH_1 = styph_volt_params(:,good_styph_indx_1)';
% --------------------------------------------------------
inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\styph\transferFunctionStudy\inputVOLT';
file_nm = 'volt_input_model_parameters.mat';
load([inDIR filesep file_nm]);

styph_volt_params = model_est_pvecs(2:end-1,:);
styph_members = 2*ones(1,size(styph_volt_params,2));
styph_fit = model_est_fit;

numStyph_total_2 = size(model_est_pvecs,2);
good_styph_indx_2 = find(styph_fit > model_fit_threshold);
numStyph_good_2 = length(good_styph_indx_2);
GOOD_STYPH_2 = styph_volt_params(:,good_styph_indx_2)';

% --------------------------------------------------------

saveDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\';


% %{

%% tSNE - visualization:

% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 8;

figure;
num_ecoli_good = numEcoli_good_1 + numEcoli_good_2;
num_styph_good = numStyph_good_1 + numStyph_good_2;
cellFamily_membership = [1*ones(num_ecoli_good,1); 2*ones(num_styph_good,1)];
GOOD_DATA = [GOOD_ECOLI_1; GOOD_ECOLI_2; GOOD_STYPH_1; GOOD_STYPH_2]; % row-vectors

DATA = tsne(GOOD_DATA, cellFamily_membership, out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params');

% scatter plot for all the data
fig = figure;
colormap jet
h1 = scatter3(DATA(1:numEcoli_good_1,1),DATA(1:numEcoli_good_1,2),DATA(1:numEcoli_good_1,3),'filled','red');
hold on
h2 = scatter3(DATA(1+numEcoli_good_1:num_ecoli_good,1),DATA(1+numEcoli_good_1:num_ecoli_good,2),DATA(1+numEcoli_good_1:num_ecoli_good,3),'filled','magenta');
hold on
h3 = scatter3(DATA(1+num_ecoli_good:numStyph_good_1+num_ecoli_good,1),DATA(1+num_ecoli_good:numStyph_good_1+num_ecoli_good,2),DATA(1+num_ecoli_good:numStyph_good_1+num_ecoli_good,3),'filled','blue');
hold on
h4 = scatter3(DATA(1+numStyph_good_1+num_ecoli_good:end,1),DATA(1+numStyph_good_1+num_ecoli_good:end,2),DATA(1+numStyph_good_1+num_ecoli_good:end,3),'filled','cyan');
hold off
title(PLOT_TITLE)


figName = [saveDIR filesep 'twoday_model_est_params.fig'];
savefig(figName);

plt_nm = sprintf('twoday_model_est_params.png');
saveas(fig, fullfile(saveDIR, plt_nm), 'png');

% %}

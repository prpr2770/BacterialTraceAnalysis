clc; close all; clear all;

% %{

%%

%{
inputType = 'volt'; %'volt' 'ca'

cellFamily = 'styph'; %'styph'
count = 0;

while count < 2
    switch cellFamily
        case 'ecoli'
            clearvars -except cellFamily count inputType
            sprintf('Processing>> inputType: %s, cellFamily: %s', inputType, cellFamily)
                   
            inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\ecoli';
            file_nm = 'ecoli_volt_ca_sigs.mat';
            load([inDIR filesep file_nm]);
            
            transferFunctionStudy
            count = count + 1;
                        
        case 'styph'
            clearvars -except cellFamily count inputType
            sprintf('Processing>> inputType: %s, cellFamily: %s', inputType, cellFamily)
            
            inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\styph';
            file_nm = 'styph_volt_ca_sigs.mat';
            load([inDIR filesep file_nm]);
            
            transferFunctionStudy

            
            cellFamily = 'ecoli';
            count = count + 1;
    end % switch
    
end % while

clear all; close all; 

% %}
%%

inputType = 'ca'; %'volt' 'ca'

cellFamily = 'styph'; %'styph'
count = 0;

homeDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun';

while count < 2
    switch cellFamily
        case 'ecoli'
            clearvars -except cellFamily count inputType homeDIR
            sprintf('Processing>> inputType: %s, cellFamily: %s', inputType, cellFamily)
            
            inDIR = [homeDIR filesep cellFamily];
            file_nm = [cellFamily '_volt_ca_sigs.mat'];
            load([inDIR filesep file_nm]);
            
            transferFunctionStudy
            count = count + 1;
            
        case 'styph'
            clearvars -except cellFamily count inputType homeDIR
            sprintf('Processing>> inputType: %s, cellFamily: %s', inputType, cellFamily)
            
            inDIR = [homeDIR filesep cellFamily];
            file_nm = [cellFamily '_volt_ca_sigs.mat'];
            load([inDIR filesep file_nm]);
            
            transferFunctionStudy
            
            
            cellFamily = 'ecoli';
            count = count + 1;
    end % switch
    
end % while

% close all; clear all; 

%}

% %{

%% tSNE Plots of the transfer function parameters

model_fit_threshold = 50;

inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\ecoli\transferFunctionStudy\inputVOLT';
file_nm = 'volt_input_model_parameters.mat';
load([inDIR filesep file_nm]);

numEcoli_total = size(model_est_pvecs,2);
ecoli_volt_params = model_est_pvecs(2:end-1,:);
ecoli_members = ones(1,size(ecoli_volt_params,2));


ecoli_fit = model_est_fit;
good_ecoli_indx = find(ecoli_fit > model_fit_threshold);

numEcoli_good = length(good_ecoli_indx);
GOOD_ECOLI = ecoli_volt_params(:,good_ecoli_indx)';


% --------------------------------------------------------
inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\styph\transferFunctionStudy\inputVOLT';
file_nm = 'volt_input_model_parameters.mat';
load([inDIR filesep file_nm]);

numStyph_total = size(model_est_pvecs,2);

styph_volt_params = model_est_pvecs(2:end-1,:);
styph_members = 2*ones(1,size(styph_volt_params,2));
styph_fit = model_est_fit;
good_styph_indx = find(styph_fit > model_fit_threshold);

numStyph_good = length(good_styph_indx);
GOOD_STYPH = styph_volt_params(:,good_styph_indx)';


saveDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\';
%% tSNE - visualization

% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 8;

figure;
GOOD_DATA = [ecoli_volt_params(:,good_ecoli_indx) styph_volt_params(:,good_styph_indx)]'; % row-vectors
DATA = tsne(GOOD_DATA, [ecoli_members(:,good_ecoli_indx) styph_members(:,good_styph_indx)]', out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params');

% scatter plot for all the data
fig = figure;
colormap jet
h1 = scatter3(DATA(1:length(good_ecoli_indx),1),DATA(1:length(good_ecoli_indx),2),DATA(1:length(good_ecoli_indx),3),'filled','red');
hold on
h2 = scatter3(DATA(1+length(good_ecoli_indx):end,1),DATA(1+length(good_ecoli_indx):end,2),DATA(1+length(good_ecoli_indx):end,3),'filled','blue');
hold off
title(PLOT_TITLE)


figName = [saveDIR filesep 'model_est_params.fig'];
savefig(figName);

plt_nm = sprintf('model_est_params.png');
saveas(fig, fullfile(saveDIR, plt_nm), 'png');

% %}



%% Exploring for geographical clustering

%{
1. kmeans: identify 4 clusters in the figure.
2. extract IDs for the cells. 
3. Obtain their original volt-calcium dynamics, and then extract tfest()
for the collection of data. 
    + Compare the parameters extracted, to the average model obtained via
    kmeans clustering. 

%}

% kmeans of DATA
[km_cluster_idx, km_cluster_centers] = kmeans(DATA,4);

ecoli_good_clusterID = zeros(1,length(GOOD_ECOLI));
styph_good_clusterID = zeros(1,length(GOOD_STYPH));

% determine clusterIDs for the good_ecoli and good_styph data.
countEcoli = 0;
countStyph = 0;
for dataID = 1:length(km_cluster_idx)
   if dataID <= length(GOOD_ECOLI)
       countEcoli = countEcoli + 1;
       ecoli_good_clusterID(countEcoli) = km_cluster_idx(dataID);
   else
       countStyph = countStyph + 1;
       styph_good_clusterID(countStyph) = km_cluster_idx(dataID);
   end
end

% Ecoli Data

ecoli_good_clusterID  % clusterID for each of the ecoli cells. 
GOOD_ECOLI % model parameters of the good ecoli.
good_ecoli_indx % index of each Good_Ecoli in original dataset. 

% Identify original cellID for each cell in cluster
members_cluster_1 = [];
members_cluster_2 = [];
members_cluster_3 = [];
members_cluster_4 = [];

for idx = 1:length(ecoli_good_clusterID)
   switch ecoli_good_clusterID(idx)
       case 1
           members_cluster_1 = [members_cluster_1 good_ecoli_indx(idx)];
       case 2
           members_cluster_2 = [members_cluster_2 good_ecoli_indx(idx)];
       case 3
           members_cluster_3 = [members_cluster_3 good_ecoli_indx(idx)];
       case 4
           members_cluster_4 = [members_cluster_4 good_ecoli_indx(idx)];
   end
end

%%
% Extract the volt-cal signals for cells in a cluster
inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\ecoli';
file_nm = 'ecoli_volt_ca_sigs.mat';


[pvec_1, dpvec_1] = getTransferFunctionParams_clusterMembers(members_cluster_1,inDIR, file_nm)
[pvec_2, dpvec_2] = getTransferFunctionParams_clusterMembers(members_cluster_2,inDIR, file_nm)
[pvec_3, dpvec_3] = getTransferFunctionParams_clusterMembers(members_cluster_3,inDIR, file_nm)
[pvec_4, dpvec_4] = getTransferFunctionParams_clusterMembers(members_cluster_4,inDIR, file_nm)


%%
% Generate the tSNE plots


% %}
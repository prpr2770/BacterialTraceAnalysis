function [pvec, dpvec] = getTransferFunctionParams_clusterMembers(members_cluster,inDIR, file_nm)
%{
Script to extract transfer function model parameters for a given group of cells.

inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun\ecoli';
file_nm = 'ecoli_volt_ca_sigs.mat';

%}

load([inDIR filesep file_nm]);
Ts = 5;

% Cluster 1 : members_cluster_1
volt_cluster = zeros(length(members_cluster),1200);
ca_cluster = zeros(length(members_cluster),1200);

for idx = 1:length(members_cluster)
    cellID = members_cluster(idx);
    volt_cluster(idx,:) = volt_sigs_noisy(cellID).sig_dn_wav;
    ca_cluster(idx,:) = ca_sigs_noisy(cellID).sig_dn_wav;
end

A = volt_cluster';
u = A - repmat(mean(A),size(A,1),1);


A = ca_cluster';
y = A - repmat(mean(A),size(A,1),1);

% y = ca_cluster' - mean(ca_cluster');

num_poles = 5;
num_zeros = 3;
iodelay = NaN; % unknown transport delay

data = iddata(y,u,Ts);

% assign the data-subsets for tasks.
estimation = data;
validation = data;  % determine validation dataset.

% estimated - model parameters.
size(u)
size(y)

model_est = tfest(estimation,num_poles, num_zeros, iodelay,'Ts',data.Ts);

[~,fit_goodness,~] = compare(estimation, model_est);
[pvec, dpvec] =  getpvec(model_est);

end
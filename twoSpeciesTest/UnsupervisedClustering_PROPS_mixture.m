% Unsupervised Clustering: PROPS_Mixed
% Compute: CodeWordHistogram_TSNE 

%%
close all; clear all;

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella\NoisyResults_18-Mar-2016\SetA\ChangeScores_20-Mar-2016\spectG_21-Mar-2016\04_122826_PROPS_200ms_900fr_0V_561_Mix_mix_3_Traces';
cd(dAnz)

flist = dir('NORM_SPECTG_*.mat');
nfiles = length(flist);

fs = 5; % Sampling frequency. 5Hz

% % Directory to save processed data
% saveDir = [dAnz filesep 'codeWords_' num2str(date)];
% if ~exist(saveDir)
%     mkdir(saveDir)
% end

% create file to store CodeWords-Histograms of each Trace
cwHist_fileName = [dAnz filesep 'cwHIST_SPECTG.mat'];
cwHist_mFile = matfile(cwHist_fileName,'Writable',true);

%% Compile all NORM_SPECTG 
ALL_NORM_SPECTG = []; % Col-Vectors
for fid = 1:nfiles
   fdata = load(flist(fid).name);
   ALL_NORM_SPECTG = [ALL_NORM_SPECTG fdata(1).NORM_SPECTG];
end

%% Use K-Means clustering to detect Code-Words.

numCodeWords = 256; 
iter = 20;                      % for knn-ma-toolbox.

CODEWORDS = ma_kmeans(ALL_NORM_SPECTG', iter, numCodeWords);
CODEWORDS = CODEWORDS'; % Convert into Col-vectors.

cwHist_mFile.CODEWORDS = CODEWORDS;
%% Compute CODEWORD_histogram for each trace.

tau = 5;
arr_ALL_CW_HIST = [];

for fid = 1:nfiles
   fdata = load(flist(fid).name);
   DATA = fdata(1).NORM_SPECTG;     % col-vectors
   numFrames = size(DATA,2);
   
   nbrs_of_Frames = knnsearch(CODEWORDS',DATA','k',tau,'distance','euclidean');
   nbrs_of_Frames = reshape(nbrs_of_Frames,[],1); % col-vector 1D
   
   % extract histogram of occurence and re-structure into COL_VEC
   hist_vec = histc(nbrs_of_Frames,1:numCodeWords);
   hist_vec = reshape(hist_vec,[],1);
   hist_vec = (1/numFrames) * (1/tau)* hist_vec;
   
   % archive the code-word histogram for each trace
   ALL_CW_HIST(fid).CW_HIST = hist_vec;
   arr_ALL_CW_HIST = [arr_ALL_CW_HIST hist_vec];
end

cwHist_mFile.ALL_CW_HIST = ALL_CW_HIST;

%% tSNE - visualization

% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 30;

figure;
DATA = tsne(arr_ALL_CW_HIST', [], out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of CodeWordsHistogram-Spectogram');

% -------------------------------------------------------------------------
% scatter plot for all the data
fig = figure()
colormap jet
h = scatter3(DATA(:,1),DATA(:,2),DATA(:,3),'filled');
title(PLOT_TITLE)

figName = [dAnz filesep 'cwHist_tsne_data.fig'];
savefig(figName);

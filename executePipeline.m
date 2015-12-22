% =========================================================================
%{
pipeline: codeword-histogram + tsne

itml :> Requires Supervision!

%}
close all;clear all;


tracesDirName = 'H:\KraljLab\';

% create file to store CodeWords-Histograms of each Trace
cwHist_fileName = strcat(tracesDirName,'cwHIST_SPECTG.mat');
cwHist_mFile = matfile(cwHist_fileName);

CWHIST = cwHist_mFile.CODEWORD_HIST;        % col-vectors.

% -------------------------------------------------------------------------
% Execute tSNE

% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 30;

DATA = tsne(CWHIST, [], out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of CodeWordsHistogram-Spectogram');

% -------------------------------------------------------------------------
% scatter plot for all the data
fig = figure()
colormap jet
h = scatter3(DATA(:,1),DATA(:,2),DATA(:,3),'filled');
title(PLOT_TITLE)

figName = strcat(tracesDirName,'cwHist_tsne_data.fig');
savefig(figName);









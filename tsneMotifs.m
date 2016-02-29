% motif-mining with tSNE
% for every trace extract the W-sized windows, that should contain the motifs. 
%     Create a bag of words for one/all traces. Then run tSNE.
%     Observe output

close all; 
clear all;

try
tracesDirName = 'H:\KraljLab\';
tracesFileName = strcat(tracesDirName,'Data_20-Jan-2016.mat');
load(tracesFileName); %Data_20-Jan-2016.mat;
Fs = 5; % samplingFreq 5Hz

ORG_DATA = data;
% extract the values
tt=[]; for i = 1:size(ORG_DATA,2) tt = [tt; ORG_DATA(i).TimeTrace]; end
cntrs=[]; for i = 1:size(ORG_DATA,2) cntrs = [cntrs; ORG_DATA(i).Centroid]; end
catch
 warning('Error Reading Input Data');
end

indx = ceil(rand()*size(tt,1));
sig = tt(indx,:);
wnd = 16;

sig = sig- mean(sig);

words = [];
for id = 1:length(sig)-wnd+1
    words = [words; sig(1,id:id+wnd-1)];
end

% specWords = fft(words');
%==========================================================
% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 8;

% DATA = tsne(CWHIST, [], out_dims, initial_dims, perplexity);
figure
DATA = tsne(words, [], out_dims, initial_dims, perplexity);

size(words)
size(DATA)

scatter3(DATA(:,1),DATA(:,2),DATA(:,3));


% study intrinsic dimension of the dataset.

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


% cols = unique(ceil(rand(1,120)*120));
cols = [100 200 300 400 600 800 900];
dimEstimates = [];

for i=1:length(cols)
    col = cols(i);
X = tt(:,1:col);
dims = GetDim(X);

dimEstimates = [dimEstimates; dims];
end

plot(dimEstimates)


X=tt';
N = size(X,2);
% N = 5;
synchIndex = zeros(N);
for i = 1:N
    i
    tic
    for j= i+1:N
sig1 = X(:,i);
sig2 = X(:,j);
cpr = phasesynchro(sig1,sig2);
synchIndex(i,j) = cpr;
synchIndex(j,i) = cpr;
    end
    toc
end


fName = strcat(tracesDirName,'SynchMatrix.mat');
save(fName,'synchIndex');
figure;
imagesc(synchIndex>0.8);
colorbar

% =========================================================================
% determining threshold for detecting clusters.
threshold = 0.75 ; % this generates 2 distinct communities. 

% derive adjacency matrix for the edges satisfying threshold.
A1 = synchIndex>threshold ;
% generate plot of the graph from adjacency matrix. 
G1 = graph(A1) ;

% what does this do? Implements DFS in the Graph generated. 
L = [] ;
for k=1:N
    if ~sum(L==k) % if node k hasn't been visited then do following
        L = [L unique(dfsearch(G1,k))'] ;
    end
end

sortedSynchMatrix = synchIndex(L,L) ;

A2 = A1(L,L);
G2 = graph(A2)

% ------------- Plots
% visualize graphs generated. 
% fig4 = figure(4); 
figure;
subplot(1,2,1)
plot(G1) ;
subplot(1,2,2)
plot(G2)

% fig5 = figure(5)
figure;
colormap gray
imagesc(sortedSynchMatrix )
colorbar
% % 
% % fName = strcat(tracesDirName,'SynchMatrix_NegCorr.mat');
% % save(fName,'sortedSynchMatrix','synchIndex');


% =================================================================
weak_concomps = conncomp(G1);
totalNumberOfComponents = max(unique(weak_concomps))

memberCounts = zeros(1,totalNumberOfComponents);
for k = 1:totalNumberOfComponents 
   memberCounts(k) = sum(weak_concomps == k );
end

figure;
plot(memberCounts)

[maxMemb, maxClusterID] = max(memberCounts);
% finding componentID for each node. 
% memberCounts(1:50) % 11, 18 : major ComponentIDs

group2 = find(weak_concomps ==  maxClusterID);

posSynch = group2;
negSynch = group1;
fName = strcat(tracesDirName,'Synch_ClusterMembers.mat');
save(fName,'posSynch','negSynch');

% =================================================================
% plotting signals from Clusters

figure;
plot(X(:,posSynch),'Color',[0.8, 0.8, 0.8])

asynchGroup = setdiff((1:N),posSynch);
figure;
plot(X(:,asynchGroup ),'Color',[0.8, 0.8, 0.8])

% get 10 of the synch memebers, and plot their changescores
indx = unique(ceil(rand(1,10)*length(posSynch)));
Y = X(:,indx);
ChgScores_Y = [];

alpha = .0;
n = 50; % 2*n+k-2 is the size of the "buffer zone".
k = 10; % window width

for m = 1:length(indx)
    sig = Y(:,m)';
    
    score1 = change_detection(sig,n,k,alpha);
    score2 = change_detection(sig(:,end:-1:1),n,k,alpha);
    score2 = score2(end:-1:1); % why do this?
    cs_y = [zeros(1,n-1+k/2), (score1 + score2), zeros(1,n-1+k/2)];
    
    ChgScores_Y  = [ChgScores_Y;  cs_y];
end

figure;
plot(ChgScores_Y','Color',[0.8, 0.8, 0.8])

    
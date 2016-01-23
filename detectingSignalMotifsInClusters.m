% Script 2: Determining signal motifs in the clusters.

% Detecting whether there exists synchronization between the cell pulses
% -------------------------------------------------------------------------

close all;clear all;


tracesDirName = 'H:\KraljLab\';
tracesFileName = strcat(tracesDirName,'PROPS_data.mat');
load(tracesFileName);
fs = 5; % samplingFreq 5Hz

% datafile
A = intens;

% -----------remove the mean
numCols = size(A,2);
meanA = sum(A,2)/numCols;
meanA_mat = repmat(meanA,1,numCols);
data = A - meanA_mat;

% normalizing the traces
A = data; 
normA = sum(A.*A,2).^(0.5);
normA_mat = repmat(normA,1,numCols);
Data = A./normA_mat;
% 
% % concise way to compute norm
% norm_Data = diag(sqrt(sum(data.^2')).^(-1))*data ;
% diff_norm_data = Data - norm_Data; 

% -----------

[numTraces,numSamplePoints] = size(data);

% The below corrMatrices haven't been used to compute the Graph/etc
CorrelationMatrix = zeros(numTraces,numTraces);
CorrelationMatrix = Data*Data';


% ========================================================================
% Plots :  Finding Correlation Matrix. 

fig1 = figure(1)
subplot(2,1,1)
colormap gray
% colormap parula
imagesc(CorrelationMatrix)
colorbar
subplot(2,1,2)
colormap gray
% colormap parula
imagesc(CorrelationMatrix)
colorbar

fig2 = figure(2)
subplot(2,1,1)
colormap gray
perm = symrcm(CorrelationMatrix);
sortedCorrMatrix = CorrelationMatrix(perm,perm);
imagesc(sortedCorrMatrix)
colorbar
subplot(2,1,2)
colormap gray
perm = symrcm(CorrelationMatrix);
sortedCorrMatrix = CorrelationMatrix(perm,perm);
imagesc(sortedCorrMatrix)
colorbar

% ========================================================================
% Determine Cluster existence. 

% assigning spatial coordinates to the points
xCoords = 1:numTraces;
yCoords = ceil(rand(size(xCoords))*numTraces);
pointCoords = [xCoords' yCoords'];
size(pointCoords)


% determining threshold for detecting clusters.
threshold = 0.69 ; % this generates 2 distinct communities. 

% derive adjacency matrix for the edges satisfying threshold.
A1 = CorrelationMatrix>threshold ;
% generate plot of the graph from adjacency matrix. 
G1 = graph(A1) ;

% what does this do? Implements DFS in the Graph generated. 
L = [] ;
for k=1:numTraces
    if ~sum(L==k) % if node k hasn't been visited then do following
        L = [L unique(dfsearch(G1,k))'] ;
    end
end

sortedCorrMatrix = CorrelationMatrix(L,L) ;

A2 = A1(L,L);
G2 = graph(A2)

% ------------- Plots
% visualize graphs generated. 
fig4 = figure(4); 
subplot(1,2,1)
plot(G1) ;
subplot(1,2,2)
plot(G2)

fig5 = figure(5)
colormap gray
imagesc(sortedCorrMatrix )
colorbar

% Generate spatial plots of graph to determine clusters
fig6 = figure(6)
subplot(1,2,1)
gplot(A1,pointCoords)
subplot(1,2,2)
gplot(A2,pointCoords)


% ========================================================================
% Determine Cluster existence and cluster nodes

weak_concomps = conncomp(G2);
totalNumberOfComponents = max(unique(weak_concomps))

memberCounts = zeros(1,totalNumberOfComponents);
for k = 1:totalNumberOfComponents 
   memberCounts(k) = sum(weak_concomps == k );
end

figure;
plot(memberCounts)

% finding componentID for each node. 
memberCounts(1:20) % 11, 18 : major ComponentIDs

group1 = find(weak_concomps == 10);
group2 = find(weak_concomps == 17);

% ========================================================================
avg_grp1 = sum(data(group1,:),1)/length(group1);
avg_grp2 = sum(data(group2,:),1)/length(group2);

figure;
subplot(2,1,1)
plot(avg_grp1)
subplot(2,1,2)
plot(avg_grp2)


figure;
subplot(2,1,1)
plot(data(group1,:))
subplot(2,1,2)
plot(data(group2,:))
% ========================================================================



% Detecting whether there exists synchronization between the cell pulses
% -------------------------------------------------------------------------

close all;clear all;


tracesDirName = 'H:\KraljLab\';
% tracesFileName = strcat(tracesDirName,'PROPS_data.mat');
tracesFileName = strcat(tracesDirName,'Data_20-Jan-2016.mat');
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
% -----------

[numTraces,numSamplePoints] = size(data);

CorrelationMatrix = zeros(numTraces,numTraces);
for row = 1:numTraces
    for col = row: numTraces
        x1 = data(row,:);
        x2 = data(col,:);
        innerProd = sum(x1.*x2)/(sum(x1.*x1)*sum(x2.*x2))^0.5; 
        CorrelationMatrix(row,col) = innerProd;
        CorrelationMatrix(col,row) = innerProd;
    end
end

% The below corrMatrices haven't been used to compute the Graph/etc
fastCorrMat = zeros(numTraces,numTraces);
fastCorrMat = Data*Data';

diffMat = fastCorrMat - CorrelationMatrix;
find(diffMat)


% concise way to compute norm
norm_Data = diag(sqrt(sum(data.^2')).^(-1))*data ;
diff_norm_data = Data - norm_Data; 


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
imagesc(fastCorrMat)
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
perm = symrcm(fastCorrMat);
sortedCorrMatrix = fastCorrMat(perm,perm);
imagesc(sortedCorrMatrix)
colorbar

% Plot difference of computing the Correlation Matrix using Matrix
% Multiplication and Loop based. 
fig3 = figure(3)
subplot(2,1,1)
colormap gray
% colormap parula
imagesc(diffMat)
colorbar
subplot(2,1,2)
colormap gray
% colormap parula
imagesc(diff_norm_data)
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

group1 = find(weak_concomps ==  maxClusterID);
group2 = find(weak_concomps == 17);

% ========================================================================
avg_group1 = sum(Data(group1,:),1)/length(group1);
avg_group2 = sum(Data(group2,:),1)/length(group2);

figure;
subplot(2,1,1)
plot(Data(group1,:)')
title('traces of members: Cmp1')
subplot(2,1,2)
plot(Data(group2,:)')
title('traces of members: Cmp2')

figure;
subplot(2,1,1)
plot(avg_group1)
title('average of traces: Cmp1')
subplot(2,1,2)
plot(avg_group2)
title('average of traces: Cmp2')

% ========================================================================
diff_orig_frm_avg_grp1 = Data(group1,:) - repmat(avg_group1,length(group1),1);
diff_orig_frm_avg_grp2 = Data(group2,:) - repmat(avg_group2,length(group2),1);

noise_grp1 = reshape(diff_orig_frm_avg_grp1,1,[]);
noise_grp2 = reshape(diff_orig_frm_avg_grp2,1,[]);

figure;
subplot(1,2,1)
histogram(noise_grp1,'Normalization','probability')
title('pdf of noise in signals from avg-signal: Cmp1')
subplot(1,2,2)
histogram(noise_grp2,'Normalization','probability')
title('pdf of noise in signals from avg-signal: Cmp2')

% alternative approach to determining histogram and pdf.
[ftshist1, binpos1] = hist(noise_grp1, 512);
[ftshist2, binpos2] = hist(noise_grp2, 512);

% ========================================================================
% determine a line-fit for the data. 
% 
% ft = fittype({'x','1'})
% data1 = [1:length(avg_group1); avg_group1]';
% fitobject = fit(data1,ft)

% best line-fit for group1
pf1 = polyfit(1:length(avg_group1),avg_group1,1)
bstFit1 = polyval(pf1,1:length(avg_group1));
figure;
plot(1:length(avg_group1),avg_group1,'r+')
hold on
plot(1:length(avg_group1),bstFit1,'g');
hold off
title('Determine best-fit curve to approximate avg-trace of Cmp1')

% best line-fit for group1
pf2 = polyfit(1:length(avg_group2),avg_group2,1)
bstFit2 = polyval(pf2,1:length(avg_group2));
figure;
plot(1:length(avg_group2),avg_group2,'r+')
hold on
plot(1:length(avg_group2),bstFit2,'g');
hold off
title('Determine best-fit curve to approximate avg-trace of Cmp2')

% ========================================================================
% find noise-distribution from best-fit line.
% Shouldn't this be GAUSSIAN, because that's how the POLYFIT works!!

diff_orig_frm_bstFt1 = Data(group1,:) - repmat(bstFit1,length(group1),1);
diff_orig_frm_bstFt2 = Data(group2,:) - repmat(bstFit2,length(group2),1);

frmBstFit_noise_grp1 = reshape(diff_orig_frm_bstFt1,1,[]);
frmBstFit_noise_grp2 = reshape(diff_orig_frm_bstFt2,1,[]);

figure;
subplot(1,2,1)
histogram(frmBstFit_noise_grp1,'Normalization','probability')
title('pdf of noise in signals from best-fit: Cmp1')
subplot(1,2,2)
histogram(frmBstFit_noise_grp2,'Normalization','probability')
title('pdf of noise in signals from best-fit: Cmp2')

% alternative approach to determining histogram and pdf.
[ftshist1, binpos1] = hist(noise_grp1, 512);
[ftshist2, binpos2] = hist(noise_grp2, 512);

% ========================================================================

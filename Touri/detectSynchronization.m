% Detecting whether there exists synchronization between the cell pulses
% -------------------------------------------------------------------------

%close all;clear all;


%tracesDirName = 'H:\KraljLab\';
%tracesFileName = strcat(tracesDirName,'PROPS_data.mat');
%load(tracesFileName);
fs = 5; % samplingFreq 5Hz

% datafile
A = intens;

% -----------remove the mean
numCols = size(A,2);
meanA = sum(A,2)/numCols;
meanA_mat = repmat(meanA,1,numCols);
data = A - meanA_mat;

% -----------

[numTraces,numSamplePoints] = size(data);
data_norm = diag(sqrt(sum(data.^2')).^(-1))*data ;
CorrelationMatrix  = data_norm*data_norm' ;
% CorrelationMatrix = zeros(numTraces,numTraces);
% for row = 1:numTraces
%     for col = row: numTraces
%         x1 = data(row,:);
%         x2 = data(col,:);
%         innerProd = (x1*x2')/(norm(x1)*norm(x2)); 
%         CorrelationMatrix(row,col) = innerProd;
%         CorrelationMatrix(col,row) = innerProd;
%     end
% end

fig1 = figure(1)
colormap gray
% colormap parula
imagesc(CorrelationMatrix)
colorbar


% fig2 = figure(2)
% colormap gray
% perm = symrcm(CorrelationMatrix);
% sortedCorrMatrix = CorrelationMatrix(perm,perm);
% imagesc(sortedCorrMatrix)
% colorbar


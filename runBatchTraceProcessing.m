% Main script calling all future scripts to be run on batch-data
close all;clear all;


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
%--------------------------------------------------------------------------
% detect MaxEventChange Distribution 

counts = 200;
ids = unique(ceil(rand(1,counts)*size(tt,1)));
traces = tt;

k = 10; %window size 10


tic
% [maxEventChangeTime_A ,allPeakIntervals_A] = getDistributionMaxEventChangeTime(traces,Fs, 40);
getDistributionMaxEventChangeTime
toc





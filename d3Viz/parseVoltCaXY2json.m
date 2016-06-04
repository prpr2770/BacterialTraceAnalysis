% Reading volt_ca_xy_data.mat and then generating json file.
clear all; close all; clc;
load('C:\Users\Prahladan\Documents\MATLAB\KraljLab\d3Viz\volt_ca_xy_data.mat');

% data : [3x 3x 95] = [drug x videoFrame x numCells]

[totalDrugs, totalFrames, totalCells] = size(data);

drugID = 2;
frameID = 1;

all_volt_sigs = [];
all_ca_sigs = [];

% struct array, with each struct representing a cell.
cellData = [];
cellID = 0;
for idx = 1:totalCells
    
    if(size(data(drugID, frameID, idx).centroid,2)>1)
        cellID = cellID + 1;
        % update cell-struct
        clear aCell;
        aCell.voltSig = data(drugID, frameID, cellID).volt;
        aCell.caSig = data(drugID, frameID, cellID).ca;
        
        aCell.xPos = data(drugID, frameID, cellID).centroid(1);
        aCell.yPos = data(drugID, frameID, cellID).centroid(2);
        aCell.cellID = cellID;
        % update cell-struct into struct-array
        cellData = [cellData aCell];
        
        % compile matrix of all signals.
        all_volt_sigs = [all_volt_sigs aCell.voltSig'];
        all_ca_sigs = [all_ca_sigs aCell.caSig'];
        
    end
end

% =========================================================================
% compute voltLinkData
% runCausality(all_volt_sigs);
cov_volt_sigs = corr(all_volt_sigs);
cov_ca_sigs = corr(all_ca_sigs);


voltLinkEdges = [];
caLinkEdges = [];

for cellID = 1:length(cellData)
    
    for nbr = cellID: length(cellData)
        
        if (cov_volt_sigs(cellID,nbr) >0.4)
            
            clear voltEdge;
            voltEdge.src   = cellID;
            voltEdge.dest = nbr;
            voltEdge.value = cov_volt_sigs(cellID,nbr);
            voltLinkEdges = [voltLinkEdges voltEdge];
        end
        
        if (cov_volt_sigs(cellID,nbr) >0.4)
            clear caEdge;
            caEdge.src   = cellID;
            caEdge.dest = nbr;
            caEdge.value = cov_ca_sigs(cellID,nbr);
            caLinkEdges = [caLinkEdges caEdge];
            
        end
    end
    
end


% saving json file
json1.cellData = cellData;
json1.voltLinkData = voltLinkEdges;
json1.caLinkData = caLinkEdges;
filename = 'C:\Users\Prahladan\Desktop\GoogleInterviewPre\Data_20-Jan-2016_frame3.json';
savejson('',json1,filename);


% compute caLinkData
% =========================================================================
%%

%{
% Randomly generating network edges!

% determine the number of cells being handled.
numCells = length(cellData);

% randomly generate edges for the network.
numEdges = ceil(rand()*numCells);

all_edges = [];
for edge_indx = 1:numEdges
    clear edge;
    edge.src = ceil(rand()*numCells);
    edge.dest = ceil(rand()*numCells);
    edge.value = rand();
    all_edges = [all_edges edge];
end

linkData_1 = all_edges;

% --------------------------------------------------------
% --------------------------------------------------------
numEdges = ceil(rand()*numCells);

all_edges = [];
for edge_indx = 1:numEdges
    clear edge;
    edge.src = ceil(rand()*numCells);
    edge.dest = ceil(rand()*numCells);
    edge.value = rand();
    all_edges = [all_edges edge];
end
linkData_2 = all_edges;

% -------------------------------------------------------------------------

% saving json file
json1.cellData = cellData;
json1.voltLinkData = linkData_1;
json1.caLinkData = linkData_2;
filename = 'C:\Users\Prahladan\Desktop\GoogleInterviewPre\Data_20-Jan-2016_frame2.json';
savejson('',json1,filename);
%}

% =========================================================================
%%




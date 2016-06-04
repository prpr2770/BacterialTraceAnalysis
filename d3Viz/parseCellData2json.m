%{
Script: Read the spatial time-series data and then export into JSON format,
to visualize inter-cell relationships.

Intended JSON formats:
var cellData = [
{cellID: 1, voltSig : [1, 2, 3, 4, 5,1, 2, 3, 4, 5], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 1, yPos: 1 },
{cellID: 2, voltSig : [1, 2, 3, 4, 5,5, 4, 3, 2, 1], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 2, yPos: 2 },
{cellID: 3, voltSig : [1, 2, 3, 4, 5,6, 7, 8, 7, 6], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 3, yPos: 3 },
{cellID: 4, voltSig : [1, 2, 3, 4, 5,5, 6, 7, 8, 9], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 4, yPos: 4 },
{cellID: 5, voltSig : [1, 2, 3, 4, 5,4, 3, 2, 1, 0], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 5, yPos: 5 }];

var linkData = [{src:1, dest:2, value:0.8},{src:3, dest:4, value:0.2},{src:1, dest:5, value:0.1}];

% In MATLAB, the above data-structure is called a CELL!

%}

% =========================================================================
% read jsondata
jsondata = loadjson('C:\Users\Prahladan\Desktop\GoogleInterviewPre\cellData.json');

% write jsondata
filename = 'C:\Users\Prahladan\Desktop\GoogleInterviewPre\cellData_frm_mat.json';
savejson('',jsondata,filename)

% =========================================================================
% reformat mat-file to json data
load('H:\KraljLab\Data_20-Jan-2016.mat');

% writing mat-data into json and reading it.
filename = 'C:\Users\Prahladan\Desktop\GoogleInterviewPre\Data_20-Jan-2016.json';
savejson('',data,filename)
voltData = loadjson('C:\Users\Prahladan\Desktop\GoogleInterviewPre\Data_20-Jan-2016.json');


% reading and parsing mat file into required json format.
cellData_1 = [];
cellData_2 = [];
cellData_3 = [];

for i=1:length(data)
    clear aCell;
    aCell.voltSig = data(i).TimeTrace;
    aCell.caSig = 16000*rand(size(aCell.voltSig));
    aCell.xPos = data(i).Centroid(1);
    aCell.yPos = data(i).Centroid(2);
    aCell.cellID = i;
    
    switch data(i).Frame
        case 1
            cellData_1 = [cellData_1 aCell];
        case 2
            cellData_2 = [cellData_2 aCell];
        case 3
            cellData_3 = [cellData_3 aCell];
        otherwise
            % fprintf('Error! in data.Frame\n');
            % warning('Error! in data.Frame\n');
            error('Error! in data.Frame\n');
    end
    
end


%%{
% randomly create edges between cells for sake of visualization
linkData_1 = [];
linkData_2 = [];
linkData_3 = [];

for frame = 1:3
    
    % eval(sprintf('cellData = cellData_%d',frame);
    clear cellData;
    switch frame
        case 1
            cellData = cellData_1;
        case 2
            cellData = cellData_2;
        case 3
            cellData = cellData_3;
    end
    % --------------------------------------------------------
    % --------------------------------------------------------
    % determine the number of cells being handled.
    numCells = length(cellData);
    
    
    %{
    all_cell_sigs = [];
    for cellID = 1:numCells
        sig = cellData(cellID).voltSig;
        all_cell_sigs = [all_cell_sigs sig'];
    end
    
    CovarMatrix = covar(all_cell_sigs);
    %}
    
%     %{
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
%     %}
    

    switch frame
        case 1
            linkData_1 = all_edges;
        case 2
            linkData_2 = all_edges;
        case 3
            linkData_3 = all_edges;
    end
    
    % --------------------------------------------------------
    % --------------------------------------------------------
    %%{
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
    
end
%%}
% -------------------------------------------------------------------------

% Extracting link-edge information:

% representt het eimeseries as columns in a matrix. After that compute the
% covariance:






% saving json file 
json1.cellData = cellData_1;
json1.voltLinkData = linkData_1;
json1.caLinkData = linkData_2;
filename = 'C:\Users\Prahladan\Desktop\GoogleInterviewPre\Data_20-Jan-2016_frame1.json';
savejson('',json1,filename);




% queryBoxes of tSNE projections
% -------------------------------------------------------------------------
% Find the x-y-z coordinates of the points. For all points that lie inside
% a cube of 10-units side, see if there's any pattern. 


% Determine center of Query
queryCenter = [10 20 30];

% xyzDims =  [-60 60];        % Expresses the range for the plotted DATA
boxSide = 10;
boxDims = (boxSide/2)*[-1 -1 -1; 1 1 1];

% Coordinates of the 
queryBox = repmat(queryCenter,2,1) + boxDims

% Determine Indices of Points lying inside QueryBox
queryPoints = [];
for i = 1:3
    validPoints = find((DATA(:,i)>= queryBox(1,i)) .*  (DATA(:,i) <= queryBox(2,i) ));
    queryPoints = [queryPoints; validPoints];
end

% Determine Unique Indices
queryPoints = unique(queryPoints);

% load traces
tracesDirName = 'H:\KraljLab\';
tracesFileName = strcat(tracesDirName,'PROPS_data.mat');
load(tracesFileName);


% plot the traces
fig2 = figure(2)
plot(intens(queryPoints,:)')
axis([0 920 0 25000])

fig3 = figure(3)
plot(intens')

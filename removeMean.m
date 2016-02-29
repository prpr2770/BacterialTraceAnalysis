function sigs = removeMean(data)
% remove dc-component from the row-vectors in data-array.

meanData = mean(data,2);
numCols = size(data,2);
meanMat = repmat(meanData,1,numCols);
sigs = data - meanMat;
end
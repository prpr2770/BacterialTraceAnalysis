function DistanceMatrix = getDistanceMatrixFromMetric(X,distance_metric);

W = distance_metric;
numPoints = length(X);
D = zeros(numPoints,numPoints);
for i = 1:numPoints
    for j = i+1:numPoints
    x_point = X(i,:);
    y_point = X(j,:);
    d_ij = (x_point - y_point)*W*(x_point - y_point)';
    if (d_ij >= 0)
        D(i,j) = d_ij.^(1/2);
        D(j,i) = D(i,j);
    else
        warning('distance are complex')
    end
    
    end
end

DistanceMatrix = D;

end
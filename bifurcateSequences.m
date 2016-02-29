function [seqA_group1, seqB_group1, seqA_group2, seqB_group2] = bifurcateSequences(seqA, seqB, bandWidth)
% group1 : elements of seqB that correspond to events of groupA, lying
% within given bandwidth.

% group2 : elements of seqB that DO NOT correspond to events of groupA, lying
% within given bandwidth.

bounds_seqA = repmat(seqA,2,1);
bounds_seqA(1,:) = bounds_seqA(1,:) - bandWidth*ones(size(seqA));
bounds_seqA(2,:) = bounds_seqA(2,:) + bandWidth*ones(size(seqA));

seqA_group1 = [];
seqA_group2 = [];

seqB_group1 = [];
seqB_group2 = [];

for k=1:length(seqB)
    
    lowerBound = [seqB(k) bounds_seqA(1,:)];
    upperBound = [seqB(k) bounds_seqA(2,:)];
    
    lowIDX = find(sort(lowerBound)==seqB(k));
    highIDX = find(sort(upperBound)==seqB(k));
    
    if lowIDX > highIDX
        seqA_group1 = [seqA_group1  seqA(highIDX)];
        seqB_group1 = [seqB_group1  seqB(k)];
    end
    
end

seqA_group2 = setdiff(seqA, seqA_group1);
seqB_group2 = setdiff(seqB, seqA_group1);

end
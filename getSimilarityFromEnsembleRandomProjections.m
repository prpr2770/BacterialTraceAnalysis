function SimilarityMatrix = getSimilarityFromEnsembleRandomProjections(rndProjections , X, numEnsembles)
% % computes similarity matrix, from the binary-hash projections of points X,
% % onto K- random planes. 

% % Input:
% % + X: WinLen x N column-vectors
% % + rndProjections: K x N col-vectors
% % 
% % Output:
% % + SimilarityMatrix

lambda = 5;
N = size(X,2);
SimilarityMatrix = zeros(N);
for i = 1: numEnsembles 

    rndProjections;
    X;
    
    % sort rndProjections lexicographically.
    decProj = bi2de(rndProjections');
    [Y, sortIDx] = sort(decProj);
    % use BANDWIDTH-Lambda to compute nearest neighbor distances. 
     
end


end
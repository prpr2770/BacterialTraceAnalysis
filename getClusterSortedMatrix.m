function SortedMatrix = getClusterSortedMatrix(Mtrx, threshold)

numTraces = length(Mtrx);

% derive adjacency matrix for the edges satisfying threshold.
A1 = Mtrx>threshold ;

% generate plot of the graph from adjacency matrix. 
G1 = graph(A1) ;

% what does this do? Implements DFS in the Graph generated. 
L = [] ;
for k=1:numTraces
    if ~sum(L==k) % if node k hasn't been visited then do following
        L = [L unique(dfsearch(G1,k))'] ;
    end
end

SortedMatrix = Mtrx(L,L) ;

end
threshold = 0.7 ;
A = CorrelationMatrix>threshold ;
G = graph(A) ;
figure; 
plot(G) ;
%[B,I] = sort(sum(A),'descend') ;
%C = CorrelationMatrix(I,I);
L = [] ;
for k=1:col
    if ~sum(L==k)
        L = [L unique(dfsearch(G,k))'] ;
    end
end
C = CorrelationMatrix(L,L) ;

%fig2 = figure(2)
figure 
colormap gray
%perm = symrcm(CorrelationMatrix);
%sortedCorrMatrix = CorrelationMatrix(perm,perm);
imagesc(C)
colorbar
function [Lnorm, SimilarityMatrix, KernelMatrix, referenceWords] = getGraphLaplacian(features_win_sigs, K, maxN)
% Compute the GraphLaplacian for the dataset.


% parameters:
sigma = 2; % Gaussian RBF Kernel
method = 'kNN'; % LSH: Locality Sensitive Hashing

% --------------------------------------------------------
% extract 2-D matrix bag-of-words representation.

% collect all features into columns of matrix
X = reshape(features_win_sigs,size(features_win_sigs,1),[],1);

N = min(maxN, size(X,2));     % size(Bag-of_windowed-code-words)

% --------------------------------------------------------
% extract : CODE-WORDS 
% selecting first N- randomly chosen elements of the bag-of-words as
% code-words.

% % p = randperm(size(X,2));
% % size(X)
% % size(p)
% % X = X(:,p(1:N)); % col-vectors
% % referenceWords = X'; % row-vecs

[IDX, referenceWords] = kmeans(X', N);

% can be done by k-means !!!

% --------------------------------------------------------
SimilarityMatrix = zeros(N,N);
KernelMatrix = zeros(N,N);
posChecked = zeros(N,N);

switch method
    
    case 'kNN'
        Y = referenceWords; % row-vectors
        
        size(Y)
        size(referenceWords)
        IDX = knnsearch(Y,referenceWords,'K',K); % size NxK
        
        %iterate over each node
        for i = 1:size(IDX,1)
            
            try
            % extract coordinates of root-node.
            m = IDX(i,1);
            
            x1 = referenceWords(m,:);
            catch
               warning('value IDX not found: %d', m) 
            end
            
            % iterate over each nbr of root.
            for j = 1:size(IDX,2)
                n = IDX(i,j);
                
                if posChecked(m,n) == 0 && posChecked(n,m) == 0
                    % compute the following, only if the pair wasn't
                    % visited earlier. 
                    
                    x2 = Y(n,:);
                    
                    dist = norm(x1-x2,2);
                    
                    % Sigmoid function
                    similarity = 1; %1/(1 + dist^2);
                    
                    % RBF-Kernel
                    kernelSim = exp(-dist^2/(2*sigma^2));
                    
                    if m ~= n
                        % insert symmetrically
                        SimilarityMatrix(m,n) = similarity;
                        SimilarityMatrix(n,m) = similarity;
                        
                        KernelMatrix(m,n) = kernelSim;
                        KernelMatrix(n,m) = kernelSim;
                    
                        posChecked(m,n) =1;
                        posChecked(n,m) =1;
                    else
                        % insert along diagonal.
                        SimilarityMatrix(m,n) = similarity;
                        KernelMatrix(m,n) = kernelSim;
                        posChecked(m,n) =1;
                    end
                    
                    
                end
                
            end %j
        end % i: iterate over each Node
        
        
        
    case 'kMeans'
        X = X'; %ensure row-vector configuration.
        [IDX, C] = kmeans(X, K);
        for idx = 1:K % iterate over clusters
            ids = find(IDX == idx);
            
            for m = 1:length(ids)
                x1 = X(ids(m),:);
                for n = m:length(ids)
                    x2 = X(ids(n),:);
                    dist = norm(x1-x2,2);
                    
                    % Sigmoid function
                    similarity = 1/(1 + dist^2);
                    
                    % RBF-Kernel
                    
                    kernelSim = exp(-dist^2/(2*sigma^2));
                    
                    if m ~= n
                        SimilarityMatrix(ids(m),ids(n)) = similarity;
                        SimilarityMatrix(ids(n),ids(m)) = similarity;
                        
                        KernelMatrix(ids(n),ids(m)) = kernelSim;
                        KernelMatrix(ids(m),ids(n)) = kernelSim;
                    else
                        SimilarityMatrix(ids(m),ids(n)) = similarity;
                        
                        KernelMatrix(ids(m),ids(n)) = kernelSim;
                    end
                    
                end
            end % iterate inside a cluster
            
            
        end % iterate over clusters
        
    case 'LSH'
        
        % %         % Create Hash-function Half-Planes
        % %         T = randn(K,winLen);        % KxW matrix
        % %         rndProjections = (T*X>=0);   % KxN matrix
        % %
        % %         % sort set of binary-arrays.
        % %         SimilarityMatrix = getSimilarityFromEnsembleRandomProjections(rndProjections , X);
        
        
        
end %switch

referenceWords = referenceWords'; % col-vectors

D = diag(sum(SimilarityMatrix,2));
L = D - SimilarityMatrix;
Lnorm = eye(size(L)) - D^(-0.5)*SimilarityMatrix*D^(-0.5);
KernelMatrix;


end
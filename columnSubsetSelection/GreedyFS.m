function S = GreedyFS(data, k, c)
% Code for Greedy Unsupervised Feature Selection 
% Author: Ahmed K. Farahat (afarahat@uwaterloo.ca)
% Date: November 4th, 2013

% If you use this code, please cite any of the following papers:
% 1.Ahmed K. Farahat, Ali Ghodsi, Mohamed S. Kamel, “An Efficient Greedy 
%   Method for Unsupervised Feature Selection”, In Proceedings of the 
%   Eleventh IEEE International Conference on Data Mining, pp.161-170, 2011
%   http://dx.doi.org/10.1109/ICDM.2011.22
% 2.Ahmed K. Farahat, Ali Ghodsi, Mohamed S. Kamel, "Efficient Greedy Feature 
%   Selection for Unsupervised Learning. Knowledge and Information Systems. 
%   35(2): 285-310, 2013
%   http://dx.doi.org/10.1007/s10115-012-0538-1

% Inputs:
%   data: An n-by-m data matrix. The rows represent features and the 
%         columns represent data points.
%   k: The number of features to be selected.
%   c (OPTIONAL): The number of feature partitions used by the 
%      partition-based method. Default: no partitions (i.e., c=n). 
%      Decreasing the number of partitions decreases the accuracy of 
%      the method but improves the run time. 

% Outputs:
%   S: The indexes of the selected features (rows)

n = size(data, 1);

% Calculate the weighted centroids of random partitions
if ~exist('c', 'var')
    c = n;
    part = data;        % No partitions
else
    % Divide features into c partitions
    M = ceil(rand(n, 1)*c);
    R = zeros(c, n);
    R([1:c:numel(R)]' + M - 1) = 1;
    clear M;
       
    % Exclude empty partitions
    vSize = sum(R, 2);
    I = vSize>0;
    R = R(I, :);
    vSize = vSize(I, :);
    c = length(vSize);
    clear I vSize;
    
    part = R*data;
    clear R;
end

% Calculate inner-products between all features and centroids of partitions
H = data*part';   

% Score variables used to select features
f = sum(H.^2, 2);
g = sum(data.^2, 2);

W = zeros(n, k);
V = zeros(c, k);
S = zeros(k, 1);
for t=1:k
    % Calculate score function
    score = f./g;
    
    score(S(1:t-1)) = 0;        % To exclude previously selected features
    score(isnan(score)) = 0;    % To exclude 0/0 entries
    score(isinf(score)) = 0;    % To exclude small-value/0 entries
    
    [~, l] = max(score);
    
    if nnz(score)==0
        W = W(:, 1:t-1);
        S = S(1:t-1, :);
        break;
    end
    if t == 1
        delta = data*data(l, :)';
        gamma = H(l, :)';
    else
        delta = data*data(l, :)' - W * W(l, :)';
        gamma = H(l, :)' - V * W(l, :)';
    end
    
    alpha_sqrt = sqrt(delta(l));
    w = delta./alpha_sqrt;
    v = gamma./alpha_sqrt;
    
    if t == 1
        r2 = 0;
    else
        r2 = W * (v'*V)';
    end
    
    if sum(delta.^2) < 1e-40 || alpha_sqrt < 1e-40      % t > rank(K)
        break;
    end
    
    r1 = H*v;
    r3 = w.*w;
    
    f = f - 2*(w .* (r1 - r2)) + (norm(v)).^2 * r3;
    g = g - r3;
    
    f(l) = 0; f(f<1e-10) = 0;
    g(l) = 0; g(g<1e-10) = 0;
    
    W(:, t) = w;
    V(:, t) = v;
    S(t) = l;
end

S = S(S>0);
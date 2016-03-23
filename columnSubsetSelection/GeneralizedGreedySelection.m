function [S W]= GeneralizedGreedySelection(target, source, k)
% Code for Greedy Generalized Column Subset Selection
% Select rows of source that best represent rows of target

% Author: Ahmed K. Farahat (afarahat@uwaterloo.ca)
% Date: October 30, 2014

% If you use this code, please cite any of the following papers:
% 1.Ahmed K. Farahat, Ali Ghodsi, Mohamed S. Kamel, "Greedy Column Subset 
%   Selection for Large-scale Data Sets."  Knowledge and Information Systems. 
%   To appear, 2014
%   http://arxiv.org/abs/1312.6838
% 2.Ahmed K. Farahat, Ali Ghodsi, Mohamed S. Kamel, "A Fast Greedy Algorithm
%   for Generalized Column Subset Selection." Workshop on Greedy Algorithms, 
%   Frank-Wolfe and Friends at the Neural Information Processing Systems, 2013 
%   http://www.afarahat.com/NIPSW13_Farahat_Greedy.pdf

% Inputs:
    % target: The pxm target matrix
    % source: The nxm source matrix
    % k: The number of rows to select
% Outputs:
    % S: The indexes of selected rows
    % W: The representation of all rows in the subspace of seleced ones

H = source*target';   

[n, p] = size(H);

f = sum(H.^2, 2);
g = sum(source.^2, 2);

W = zeros(n, k);
V = zeros(p, k);
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
        delta = source*source(l, :)';
        gamma = H(l, :)';
    else
        delta = source*source(l, :)' - W * W(l, :)';
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
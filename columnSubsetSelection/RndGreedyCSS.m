function S = RndGreedyCSS(data, k, r)
% Code for "Greedy Column Subset Selection for Large-scale Data Sets"
% Implementation of RndGreedyCSS method in [1]

% Author: Ahmed K. Farahat (afarahat@uwaterloo.ca)
% Date: June 10, 2014

% If you use this code, please cite any of the following papers:
%
%[1] Farahat, A.K., Elgohary, A., Ghodsi, A., and Kamel, M.S. (2014) 
%    Greedy Column Subset Selection for Large-scale Data Sets.  
%    Knowledge and Information Systems. 33 pages. Under Review. Accepted
%    http://arxiv.org/abs/1312.6838

% Inputs:
%   data: An m-by-n data matrix. The code selects a subset of columns of 'data'.
%   k: The number of columns to be selected.
%   r (OPTIONAL): The number of random vectors used for random projection.
%      Default value is 100. If r=0, this code runs the exact GreedyCSS algorithm.

% Outputs:
%   S: The indexes of the selected columns.

if ~exist('r', 'var')
    r = 100;
end

% Prepare source and target matrices of the generalized CSS problem

% Source: Transpose the data matrix as the following code selects a subset of rows of 'source'
source = data';   

% Target: Project the columns of data into into a random subspace
n = size(data, 2);
if r == 0       % No random projection
    target = source;
    r = n;
else
    R = normrnd(0,1/sqrt(r), r, n);
    target = R*data';
end

% Calculate inner-product between source and target
H = source*target';   

% Calculate selection score variables
f = sum(H.^2, 2);
g = sum(source.^2, 2);

W = zeros(n, k);
V = zeros(r, k);
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
    
    % Update delta and gamma
    if t == 1
        delta = source*source(l, :)';
        gamma = H(l, :)';
    else
        delta = source*source(l, :)' - W * W(l, :)';
        gamma = H(l, :)' - V * W(l, :)';
    end
    
    % Calculate w and v
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
    
    % Update f and g
    f = f - 2*(w .* (r1 - r2)) + (norm(v)).^2 * r3;
    g = g - r3;
    
    f(l) = 0; f(f<1e-10) = 0;
    g(l) = 0; g(g<1e-10) = 0;
    
    W(:, t) = w;
    V(:, t) = v;
    S(t) = l;
end

% Return the selected indexes
S = S(S>0);


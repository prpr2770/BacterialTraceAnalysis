clear;

dim=2;
samples=10;

c=3; s=1;

%% 
% Generate Data-set using GaussianMixtureModels

% Mean-centers of the clusters. % Component means
myu = [0 c zeros(1,dim-2);-c -c zeros(1,dim-2);c -c zeros(1,dim-2)];   

% Concatatenate Identity Matrices
Smat = cat(3,eye(dim),eye(dim),eye(dim)); % component covariance

% What's this about? 
p = [0.4,0.3,0.3];  % component mixing proportion

% Create the GMM Model - Object.
GMMobj = gmdistribution(myu,Smat,p);

% Extract random samples from the distribution. 
[X,idx]=random(GMMobj,samples); 

X=X'; % Col-vector representation for each data-point.

%% Perform the LSLDG Clustering.

% op.bnum=min(op.samples,100); % number of basis functions
% op.sigma_list=logspace(-1,1,10); % candidates for the Gaussian width parameter
% op.lambda_list=logspace(-2,1,10); % candidates for the regularization parameter
Y=LSLDGClust(X);

%% Plot the dataset. 
if dim==2; plotData(X,Y); end;
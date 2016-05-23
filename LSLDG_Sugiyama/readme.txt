1. demo.m
X : col-vector representation of all datapoints. 
Y=LSLDGClust(X);
plotData.m

2. LSLDGClust.m
% Clustering via mode seeking based on LSLDG

C: 
op: op.tol = 1e-5; op.maxiter = 200; [op.dim, op.samples] = size(X)

a) [~,theta,C,sigma]=LSLDG(X,op,C);
b) Choose a vector, Z_old = X(:,ii);
c) Determine probabilities for each x_i centered at Z_old, under Gaussian Kernel. 
d) Determine Z_new, using theta, C and above probabilities. 
e) repeat until Z stabilies, within tolerance values. 
(Is Z the fixed point for the given X(:,ii)? It seems like it!!)


3. LSLDG.m
[g,theta,C,sigma,lambda,cind]=LSLDG(X,op,C)
% Estimating Log-Density Gradients

C >> Center points of kernels extracted randomly from the dataset. 
op.bnum = min(op.samples,100) % Number of Basis Functions.
op.sigma_list = logspace(-1,1,10)
op.lambda_list = logspace(-2,1,10)
op.cvfold = 5

XC_dist : distance of all data-vectors from centers. 
MedDim : Compute Median of Data for each dimension.  % why is this necessary? used to determine sigma for the gaussian kernel. 

cv_fold_iteration:
	psi_train
	phi_train
	K_train
	h_train

	psi_test
	phi_test
	K_test
	h_test

	for lambda_iteration:
		theta_train = linsolve(K_train + Lambda*I, h_train)
		term1 = theta_train' * K_test* theta_train
		term2 = theta_train'* h_test
		score = term1 - 2*term2

determine: 
	sigma_best
	lambda_best


% Compute these values from Data in the Reduced Dimensional Space!! 
theta = computeTheta(X,C,sigma,lambda,op);
g = grad(theta,X,C,sigma,op);

===========================================================================================
Strategy:
1. For given X_HD \in R^D, compute sigma, lambda by cross-validation. 
Iteratively do the following: 
2. For selected reduced_dim points X_LD \in R^d, 
compute Theta : of the cost-function representation via Kernel methods. 
compute Gradients : of the cost function
3. Update the coordinates: How is that done?
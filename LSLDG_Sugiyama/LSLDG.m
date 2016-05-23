function [g,theta,C,sigma,lambda,cind]=LSLDG(X,op,C)
%
% Estimating Log-Density Gradients
%
% X: (dim) by (sample) matrix
% op: options
% C: center points of kernels
%

[op.dim,op.samples]=size(X);

if ~isfield(op,'bnum') % number of basis functions
    op.bnum=min(op.samples,100);
end

if ~isfield(op,'sigma_list') % candidates for the Gaussian width parameter
    op.sigma_list=logspace(-1,1,10);
end

if ~isfield(op,'lambda_list') % candidates for the regularization parameter
    op.lambda_list=logspace(-2,1,10);
end

if ~isfield(op,'cvfold') % 5 fold cross-validation by a default setting
    op.cvfold=5;
end

%%
% Extract C: Center points of the Kernels. 
if ~exist('C','var') || isempty(C)
    cind=randperm(op.samples,op.bnum);
    C=X(:,cind);
else
    cind=[];
end

%%

% distance of all data-vectors from the centers!
XC_diff=repmat(permute(C,[2,3,1]),[1,op.samples,1])...
    -repmat(permute(X,[3,2,1]),[op.bnum,1,1]);
XC_dist=sum(XC_diff.^2,3);

% parameters for cross-validation
cv_fold=(1:op.cvfold);
cv_split=floor((0:op.samples-1)*op.cvfold./op.samples)+1;
cv_index=cv_split(randperm(op.samples));

% Compute Median of the Data - for each dimension
MedDim=MedianDiffDim(X);

sigma=zeros(op.dim,1); lambda=zeros(op.dim,1);


for dd=1:op.dim
    %%  Iterate over the dimensions of the data-set
    score_cv=zeros(length(op.sigma_list),length(op.lambda_list),length(cv_fold));
    
    
    for sigma_index=1:length(op.sigma_list)
        %% Iterate Over Sigma
        % Choose a sigma from list, and multiply with Median for given
        % dimension.
        sigma_tmp=MedDim(dd)*op.sigma_list(sigma_index);
        
        GauKer=exp(-XC_dist/(2*sigma_tmp^2));
        
        
        for kk=cv_fold
            %% Iterate over the number of times of cross-validation folds
            psi_train=XC_diff(:,cv_index~=kk,dd).*GauKer(:,cv_index~=kk);
            phi_train=(1-XC_diff(:,cv_index~=kk,dd).^2/sigma_tmp^2).*GauKer(:,cv_index~=kk);
            
            K_train=psi_train*psi_train'/size(psi_train,2);
            h_train=mean(phi_train,2);
            
            psi_test=XC_diff(:,cv_index==kk,dd).*GauKer(:,cv_index==kk);
            phi_test=(1-XC_diff(:,cv_index==kk,dd).^2/sigma_tmp^2).*GauKer(:,cv_index==kk);
            
            K_test=psi_test*psi_test'/size(psi_test,2);
            h_test=mean(phi_test,2);
            

            for lambda_index=1:length(op.lambda_list)
            %% Iterate Over Lambda
                lambda_tmp=op.lambda_list(lambda_index);
                
                thetah=linsolve(K_train+lambda_tmp*eye(size(K_train)),h_train);
                
                term1=thetah'*K_test*thetah;
                term2=thetah'*h_test;
                
                score_cv(sigma_index,lambda_index,kk)=term1-2*term2;
            end %lambda
        end % kk
    end % sigma
    
    [score_cv_tmp,lambda_index]=min(mean(score_cv,3),[],2);
    [~,sigma_index]=min(score_cv_tmp);
    
    % best values of lambda for given dimension, after cross-validation.
    lambda(dd)=op.lambda_list(lambda_index(sigma_index));
    sigma(dd)=MedDim(dd)*op.sigma_list(sigma_index);
    fprintf('dd=%g, sigma=%g, lambda=%g\n',dd,sigma(dd),lambda(dd));
end % dd

%% Critical parts
theta=computeTheta(X,C,sigma,lambda,op);
g=grad(theta,X,C,sigma,op);
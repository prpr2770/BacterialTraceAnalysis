function Y=LSLDGClust(X,op,C)
%
% Clustering via mode seeking based on LSLDG.
%

narginchk(1,3)

[op.dim,op.samples]=size(X);

if nargin<2 || ~isfield(op,'tol')
    op.tol=1e-5;
end

if nargin<2 || ~isfield(op,'maxiter')
    op.maxiter=200;
end

if nargin==3; 
    [~,theta,C,sigma]=LSLDG(X,op,C);
else
    [~,theta,C,sigma]=LSLDG(X,op);
end

%%
% What's happening here?
% Determine Modes
Y=zeros(size(X));

for ii=1:op.samples    
    Z=X(:,ii);
    iter=0; criteria=1; Ncall=0;
    while criteria
        iter=iter+1;
        dist = sum(bsxfun(@minus,Z,C).^2,1);
        
        % Find the Probabilities for each x_i centered at Z
        Gau=exp(-bsxfun(@times,dist,1./(sigma.^2))  /2.0);

        Zold=Z;
        
        % What does this do?
        % Is this the gradient descent? <theta,C> is obtained from LSLDG.
        Z=sum(theta.*C.*Gau,2)./sum(theta.*Gau,2);
        
        % For irregular samples
        if sum(isnan(sum(Z,1)))~=0 || sum(isinf(sum(Z,1)))~=0; Ncall=1; break; end;
        
        % convergence criteria
        d=sqrt(sum((Z-Zold).^2));                
        criteria = (d > op.tol) & (iter < op.maxiter);
    end
        
    % For irregular samples
    if Ncall; Y(:,ii)=NaN; else Y(:,ii)=Z; end;
end        


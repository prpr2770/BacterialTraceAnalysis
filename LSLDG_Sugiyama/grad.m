function g=grad(theta,X,C,sigma,op)
%
% Computing gradients 
%

narginchk(5,5);

n=size(X,2);

XC_diff=repmat(permute(C,[2,3,1]),[1,n,1])...
    -repmat(permute(X,[3,2,1]),[op.bnum,1,1]);
XC_dist=sum(XC_diff.^2,3);

GauKer3D=exp(-bsxfun(@times,XC_dist,1./(2*permute(sigma.^2,[2,3,1]))));
psi=XC_diff.*GauKer3D;
    
g=permute(sum(bsxfun(@times,theta',permute(psi,[1,3,2]))),[2,3,1]);

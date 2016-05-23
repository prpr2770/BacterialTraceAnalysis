function theta=computeTheta(X,C,sigma,lambda,op)
%
% Computing theta
%

XC_diff=repmat(permute(C,[2,3,1]),[1,op.samples,1])...
    -repmat(permute(X,[3,2,1]),[op.bnum,1,1]);
XC_dist=sum(XC_diff.^2,3);

theta=zeros(op.bnum,op.dim);
for dd=1:op.dim
    GauKer=exp(-XC_dist/(2*sigma(dd)^2));

    psi=XC_diff(:,:,dd).*GauKer;
    phi=(1-XC_diff(:,:,dd).^2/sigma(dd)^2).*GauKer;

    K=psi*psi'/size(psi,2);
    h=mean(phi,2);
        
    theta(:,dd)=linsolve(K+lambda(dd)*eye(size(K)),h);
end

theta=theta';

function r = get_r(alpha,m,W,nu, beta,X)
% This function calculates a matrix r which are parameters of the
% distribution over assignments for each row of the data matrix X.  The
% upated r is based on the values alpha, m, W, nu, and beta.
%
%@param alpha       : K x 1 matrix of positive dirichlet parameters
%@param m           : D x K matrix of means
%@param W           : K long cell array of D x D covariance matrics 
%@param nu          : K x 1 matrix of degrees of freedom for W matrices
%@param beta        : K x 1 matrix of scaling factors for NIW distributions
%@param X           : N x D data matrix
%
%@return r          : N x K matrix for distribution of z_i

K=size(alpha,1);
N=size(X,1);
D=size(X,2);

r=zeros(N,K);

psisum=psi(sum(alpha));

for n=1:N
    for k=1:K
        epi=exp(psi(alpha(k))-psisum);
        sumpsi=sum(psi(.5*(nu(k)+1-(1:D))));
        pre=exp(sumpsi + D*log(2) + log(det(W{k})) );
        
        power=D/(-2*beta(k))-nu(k)/2*((X(n,:)-m(:,k)')*W{k}*(X(n,:)-m(:,k)')');
        r(n,k)=epi*sqrt(pre)*exp(power);
    end
    
    %normalization
    sumri=sum(r(n,:));
    r(n,:)=r(n,:)/sumri;
end
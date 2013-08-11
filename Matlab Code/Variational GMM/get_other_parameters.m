function [alpha,m,W,nu,beta] = get_other_parameters(r, X)
% This function calculated the updated parameter values for alpha, m, W, 
% nu, and beta given the value of r and the data matrix X.
%
%@param r           : n x K matrix for distribution of z_i
%@param X           : n x D data matrix
%
%@return alpha       : K x 1 matrix of positive dirichlet parameters
%@return m           : D x K matrix of means
%@return W           : K long cell array of D x D covariance matrics 
%@return nu          : K x 1 matrix of degrees of freedom for W matrices
%@return beta        : K x 1 matrix of scaling factors for NIW distributions

global m_0 b_0 a_0 W_0 nu_0

N=size(X,1);
D=size(X,2);
K=size(r,2);

alpha=zeros(K,1);
m=zeros(D,K);
W = cell(1,K);
nu=zeros(K,1);
beta=zeros(K,1);

Nk=sum(r,1);
Xk=zeros(D,K);
Sk=cell(1,K);

for k=1:K
    sumrx=0;
    sumrxm=0;
    for n=1:N
        sumrx=sumrx+r(n,k)*X(n,:);
    end
    Xk(:,k)=sumrx/Nk(k);
    for n=1:N
        sumrxm=sumrxm+r(n,k)*(X(n,:)-Xk(:,k)')'*(X(n,:)-Xk(:,k)');
    end
    Sk{k}=sumrxm/Nk(k);
end

winv=inv(W_0);

for k=1:K
    alpha(k)=a_0+Nk(k);
    beta(k)=b_0+Nk(k);
    m(:,k)=(b_0*m_0+Nk(k)*Xk(:,k))/beta(k);
    nu(k)=nu_0+Nk(k);
    
    temp=winv+Nk(k)*Sk{k}+(b_0*Nk(k)/(b_0+Nk(k)))*(Xk(:,k)-m_0)*(Xk(:,k)-m_0)';
    
    W{k}=inv(temp);
    
end




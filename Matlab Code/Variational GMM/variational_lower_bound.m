function lb = variational_lower_bound(r,alpha,m,W,nu,beta,X)
% This function calculates the variational lower bound.  This value should
% go up as the algorithm progresses.  The algorithm only stops when the
% variational lower bound has stopped increasing.
%
%@param r           : N x K matrix for distribution of z_i
%@param alpha       : K x 1 matrix of positive dirichlet parameters
%@param m           : D x K matrix of means
%@param W           : K long cell array of D x D covariance matrics 
%@param nu          : K x 1 matrix of degrees of freedom for W matrices
%@param beta        : K x 1 matrix of scaling factors for NIW distributions
%@param X           : N x D data matrix
%
%@return lb         : calculated scalar lower bound

global m_0 b_0 a_0 W_0 nu_0

N=size(X,1);
K=size(r,2);
D=size(X,2);

alpha0=ones(K,1)*a_0;

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

lnpre=zeros(K,1);
sumalpha=sum(alpha);
lnpi=psi(alpha)-psi(sumalpha);

for k=1:K
    sumpsi=0;
    for h=1:D
        sumpsi=sumpsi+psi((nu(k)+1-h)/2);
    end
    lnpre(k)=sumpsi + D*log(2) + log(det(W{k})) ;
end


E1=0;
for k=1:K
    E1=E1+Nk(k)*(lnpre(k)-D/beta(k)-nu(k)*trace(Sk{k}*W{k})- nu(k)*(Xk(:,k)-m(:,k))'*W{k}*(Xk(:,k)-m(:,k))-D*log(2*pi));
end
E1=E1/2;

temp=0;
for n=1:N
    for k=1:K
        temp=temp+r(n,k)*lnpi(k);
    end
end

E2=temp;


dirconst=  log(gamma(sum(alpha0))/prod(gamma(alpha0)));
E3=dirconst+(a_0(1)-1)*sum(lnpi);
p1=0;
for k=1:K
    p1=p1+0.5*(D*log(b_0/(2*pi))+lnpre(k)-D*b_0/beta(k)-b_0*nu(k)*(m(:,k)-m_0)'*W{k}*(m(:,k)-m_0))... 
        +.5*(nu_0-D-1)*lnpre(k)...
        -.5*nu(k)*trace(inv(W_0)*W{k});
end


E4=p1+K*log(B(W_0,nu_0));

lnr=log(r);
E5=0;
for n=1:N
    for k=1:K
        E5=E5+r(n,k)*lnr(n,k);
    end
end

dirconst=  log(gamma(sum(alpha))/prod(gamma(alpha)));
E6=sum((alpha-1).*lnpi)+dirconst;

E7=0;
pre=exp(lnpre);
for k=1:K;
    
    E7=E7+lnpre(k)/2+D/2*log(beta(k)/(2*pi))-D/2-H(W{k},nu(k));
end

lb=E1+E2+E3+E4-E5-E6-E7;


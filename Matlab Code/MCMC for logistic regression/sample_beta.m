function beta = sample_beta(beta, sigma, X, Y, a, b, mu)
%@param beta  : k x 1 matrix column matrix of coefficient parameters of linear model
%@param sigma : scalar, variance of beta vector in prior
%@param X     : n x k design matrix
%@param Y     : n x 1 column vector of binary ({0,1}) outcomes
%@param a     : scalar prior parameter in gamma prior on precision
%@param b     : scalar prior parameter in gamma prior on precision
%@param mu    : k x 1 vector which is mean of prior on beta
%
%@param beta  : sampled values of beta

%using metroplis 
p=unifrnd(0,1);
beta_new=mvnrnd(beta,eye(8))';
P2=mvnpdf(beta,mu,sigma*eye(8))*prod((Y.*exp(X*beta)+(1-Y))./(1+exp(X*beta)));
P1=mvnpdf(beta_new,mu,sigma*eye(8))*prod((Y.*exp(X*beta_new)+(1-Y))./(1+exp(X*beta_new)));
if (p<=min(1,P1/P2))
    beta=beta_new;
else
    beta=beta;
end
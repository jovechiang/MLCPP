function log_likelihood = logistic_log_likelihood(beta, X, Y)
%Returns the loglikelihood term corresponding to the data part of the
%model.  This added to the log likelihood from the prior gives the joint
%log likelihood of the data
%
%@param beta            :   k x 1 vector of regression coefficients
%@param X               :   n x k matrix of covariates
%@param Y               :   n x 1 binary response vector
%
%@return log_likelihood :   log likelihood of response given beta and X

log_likelihood=0;
n=size(X,1);
for i=1:n
    power = X(i,:)*beta;
    prob=exp(power)/(1+exp(power));
    if(Y(i)==1)
        log_likelihood=log_likelihood+log(prob);

    else
        log_likelihood=log_likelihood+log(1-prob);

    end
end
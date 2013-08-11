function sigma = sample_sigma(beta, a, b)
%@param beta    :   k x 1 vector of regression parameters
%@param a       :   alpha parameter in gamma prior on inv sigma
%@param b       :   beta parameter in gamma prior on inv sigma
%
%@return sigma  :   sigma sampled from the condititional distribtion


%integration  
k=prod(normpdf(beta,0,sqrt(beta'*beta/8)));
while true
   sigma = 1./gamrnd(a,b);
   x= unifrnd(0,1);
   if (x<=prod(normpdf(beta,0,sigma))/k)
       break;
   end
end
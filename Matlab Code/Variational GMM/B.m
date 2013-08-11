function ret=B(W,nu)

d=size(W,1);
ret=power(det(W),-nu/2)*((2^(nu*d/2)*pi^(d*(d-1)*.25))*prod(gamma(0.5*(nu+1-(1:d)))))^-1;


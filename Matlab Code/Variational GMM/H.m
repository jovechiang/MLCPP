function ret=H(W,nu)

d=size(W,1);

ret=-log(B(W,nu))-.5*(nu-d-1)*(sum(psi(.5*(nu+1-(1:d))))+d*log(2)+log(det(W)))+.5*nu*d;
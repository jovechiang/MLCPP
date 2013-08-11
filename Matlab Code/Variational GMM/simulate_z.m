function z_sample = simulate_z(r)
% This function simulates Z values from the optimized Q distribution.
% Since the Z_i are independent, all we need is the r matrix and we can
% sample each Z_i independently from the multinomial distribution in the
% i'th row of r.
%
%@param r   : n x k matrix with rows = distribution for z_i
%
%@return z  : n x 1 sample of cluster allocations for each data point

n=size(r,1);
k=size(r,2);
z_sample=zeros(n,1);

for i=1:n
    pdf=round(r(i,:)*10);
    scope=sum(pdf);
    randnum=unidrnd(scope);
    cum=0;
    for j=1:k
        cum=cum+pdf(j);
        if(cum>=randnum)
            z_sample(i)=j;
            break;
        end
    end
end

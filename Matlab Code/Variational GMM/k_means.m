function assignments = k_means(k,data)
% This function returns one of k assignment variables for each data point
% in a given data matrix using the k means algorithm.
%
% @param k    : scalar value indicating the number of groups
% @param data : n x d matrix of data
%
% @return assignments : n x 1 matrix of cluster assignments in 1:k

n=size(data,1);
d=size(data,2);

assignments=zeros(n,1);

mu=zeros(k,d);
perm=randperm(n);
for i=1:k
    mu(i,:)=data(perm(i),:);
end
rounds=20;

while(true)
    [newmu assignments]=update_cluster_centers(mu,data);
    if(newmu==mu)
        break;
    end
    mu=newmu;
end


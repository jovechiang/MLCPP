function er = expected_rand_index(data, r, labels)
% This function calculates the expected rand index given the distributions
% over assignments represented by r.  The rand index is an indication of
% how often two data points which should be in the same group are in
% different groups and two data points which should be in seperate groups
% are actually found to be in the same group.
%
%@param data    : n x d data matrix
%@param r       : n x k matrix of distributions over z
%@param labels  : n x 1 matrix of labels

n=size(data,1);
d=size(data,2);
k=size(r,2);

indicator=0;
for i=1:n
    for j=i+1:n

        sameprob=sum(r(i,:).*r(j,:));
        
        if(labels(i)==labels(j))
            indicator=indicator+1-sameprob;
        else
            indicator=indicator+sameprob;
        end
    end
end

er=indicator/(n*(n-1)/2);

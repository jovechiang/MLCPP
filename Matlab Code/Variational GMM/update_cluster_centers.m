function [new_centers assignments] = update_cluster_centers(cluster_centers,data)
% This function takes a data set, and a set of cluster centers and returns
% both the updated set of cluster centers and a center assignment variable 
% for each data point.  The data set is assumed to be arranged such that
% each row is one data d dimensional data point.  In the documentation we
% wil use k to denote the number of clusters for the k-means analysis.
%
% @param cluster_centers : k x d matrix of cluster centers
% @param data            : n x d matrix of data
%
% @return new_centers    : k x d matrix of updated cluster centers
% @return assignments    : n x 1 matrix of cluster assignments in 1:k
%


k=size(cluster_centers,1);
n=size(data,1);
d=size(data,2);
r=zeros(n,k);
assignments=zeros(n,1);
new_centers=zeros(k,d);

for i=1:n
    x=data(i,:);
    mindis=inf;
    pos=0;
    for j=1:k
        t=x-cluster_centers(j,:);
        dis=sum(t.*t);
        if (dis<mindis)
            mindis=dis;
            pos=j;
        end
    end
    r(i,pos)=1;
    assignments(i)=pos;
end

sumr=sum(r,1);
for i=1:k
    sumdata=zeros(1,d );
    for j=1:n
        sumdata=sumdata+r(j,i)*data(j,:);
    end
    new_centers(i,:)=sumdata/sumr(i);
end
assign=sum(r,1);

for i=1:k
    if(assign(i)==0)
        disp('no assignment');
        randomnum=unidrnd(n);
        assign(i)=assign(i)+1;
        for j=1:k
            if(r(randomnum,j)==1)
                break;
            end
        end
        assign(j)=assign(j)-1;
        new_centers(i,:)=data(randomnum,:);
    end
end
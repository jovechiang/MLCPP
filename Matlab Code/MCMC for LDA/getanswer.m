
N=size(topic_assignment,2);
D=size(doc_counts,1);
K=size(topic_counts,1);
U=size(topic_counts,2);

top10=zeros(K,10);%store top10 for every topic
TOP10WO=cell(K,10);
for k=1:K
    [sortvals, sortidx] = sort(topic_counts(k,:),'descend');
    top10(k,:)=sortidx(1:10);
    TOP10WO(k,:)=WO(top10(k,:));
end

td=doc_counts;
TOP10T=cell(1,10);
for i=1:D
    td(i,:)=td(i,:)/sum(td(i,:));
end
similar=zeros(D,1);

for i=2:D
    similar(i)=dot(td(1,:),td(i,:));
end

[dsortvals,dsortdix]=sort(similar,'descend');
tensimilar=dsortdix(1:10);
tenscore=dsortvals(1:10);
TOP10T=titles(tensimilar);




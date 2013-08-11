function ll = joint_log_lik(doc_counts,topic_counts, alpha, gamma)
%Gets the joint log likelihood of the model
%
%@param doc_counts       : n_docs x n_topics vector of counts per document
%of unique topics
%@param topic_counts     : n_topics x as vector of counts per topic of unique
%words
%@param alpha            : prior dirichlet parameter on document specific
%distributions over topics
%@param gamma            : prior dirichlet parameter on topic specific
%distribuitons over words.
%
%@return ll              : joint log likelihood of model
 
%using the equation in lda lecture.pdf
D=size(doc_counts,1);
K=size(doc_counts,2);
M=size(topic_counts,2);
ll = 0;
N_d = sum(doc_counts,2);
for d = 1:D
    ll = ll + gammaln(K*alpha)- gammaln(alpha) * K;
    for k = 1:K
       ll = ll + gammaln(alpha+doc_counts(d,k));
    end
    ll = ll - gammaln(K*alpha + N_d(d));
end
M_k = sum(topic_counts,2);
for k = 1:K
    ll = ll + gammaln(M*gamma) - gammaln(gamma) * M;
    for m = 1:M
        ll = ll + gammaln(gamma+topic_counts(k,m));
        
    end
    ll = ll - gammaln(M*gamma + M_k(k));
end

function [mu,sigma,pi] = m_step_gaussian_mixture(data,gamma)
% Performs the M-step of the EM algorithm for gaussain mixture model.
%
% @param data   : n x d matrix with rows as d dimensional data points
% @param gamma  : n x k matrix of resposibilities
%
% @return mu    : d x k matrix of maximized cluster centers
% @return sigma : cell array of maximized 
%

k = size(gamma, 2);
d = size(data, 2);
n = size(gamma,1);
mu = zeros(d, k);
x = data;


nk = sum(gamma,1);

pi = zeros(1, k);
for i = 1 : k
    pi(i) = nk(i) / n;
end



for i = 1 : k
    sumn = 0;
    for j = 1 : n
        sumn = sumn + gamma(j, i) * x(j, :);
    end
    sumn = sumn / nk(i);
    mu(:,i) = sumn(:);
end



sigma = cell(1,k);
for i = 1 : k
    sumn = 0;
    for j = 1 : n
        t = x(j, :) - mu(:, i)';
        sumn = sumn + gamma(j, i) * (t' * t);
    end
    sumn = sumn / nk(i);
    %add eps by eps
    p = all(eig(sumn) > 0);
    while(p == 0)
        % matlab function eig from matlab document
        [V, D] = eig(sumn);
        %turn to be positive indefinite 
        D(D < 0) = eps;
        sumn = V * D * V';
        p = all(eig(sumn) > 0);
    end
    sigma{i} = sumn;
end



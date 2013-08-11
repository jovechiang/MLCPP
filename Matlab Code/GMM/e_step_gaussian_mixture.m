function gamma = e_step_gaussian_mixture(data,pi,mu,sigma)
% Returns a matrix of responsibilities.
%
% @param    data : data matrix n x d with rows as elements of data
% @param    pi   : column vector of probabilities for each class
% @param    mu   : d x k matrix of class centers listed as columns
% @param    sigma: cell array of class covariance matrices (d x d)
%
% @return   gamma: n x k matrix of responsibilities


k = size(mu, 2);
n = size(data, 1);
gamma = zeros(n, k);
x = data;


for i = 1 : n
    for j = 1 : k
        denom=0;
        for l = 1 : k
            denom = denom + pi(l) * mvnpdf(x(i, :), mu(:, l)', sigma{l});
        end
        gamma(i, j) = pi(j) * mvnpdf(x(i, :), mu(:, j)', sigma{j}) / denom;
    end
end
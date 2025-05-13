function R = computeR(phi, X1, X2)
% This function computes the correlation matrix btween X1 and X2
% phi - lengthscale parameters
% X1, X2 - n by d, n sample size, d dim
% correlation is calculated using squared exponential formula


[n1,d] = size(X1);
n2 = size(X2,1);

d2 = 0;
for i = 1:d
    % distance matrix multiplied by theta
    d2 = d2 + 10.^(phi(i)).*(repmat(X1(:,i),[1,n2])-repmat(X2(:,i)',[n1,1])).^2;
end

% correlation matrix
if n1 ~= n2
    R = exp(- d2);
else
    nugget = n1/(10^12 - 1); % add nugget to make sure the matrix is well conditioned
    R = exp(-d2) + nugget*eye(n1);
end
end
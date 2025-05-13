function [y_pred, y_cov] = gpPredict(model, x_pred)
% This function makes predictions based on fitted GP model
% model - the fitted GP model containing the following parameters:
% -- params_star - optimal parameters found by MLE
% -- X - a matrix of input data, each row is a sample
% -- y - a vector of response
% -- dim_qual - an index array of qualitative variables
% -- d_lv - the dimension of latent space, 1 or 2
% -- levels - array of levels of each qualitative variable
% x_pred - m by d, m samples and d dim


% without qualitative variables
if numel(fieldnames(model)) == 5
    
    phi = model.phi; X = model.X; y = model.y;
    
    n = size(X,1);
    R_xpredx = computeR(phi,x_pred, X);
    R_xx = computeR(phi, X, X);
    R_xpred = computeR(phi, x_pred, x_pred);
    
    mu = (ones(1,n)*(R_xx\ones(n,1)))\(ones(1,n)*(R_xx\y)); 
    sigma = 1/n*(y-mu)'*(R_xx\(y-mu));
    y_pred = mu + R_xpredx*(R_xx\(y-mu));
    y_cov = sigma*(R_xpred - R_xpredx*(R_xx\R_xpredx'));

% with qualitative variables
elseif numel(fieldnames(model)) == 9
    
    phi = model.phi; z = model.z; X = model.X; y = model.y;
    dim_qual = model.dim_qual; d_lv = model.d_lv; levels = model.levels;
    
    [n,~] = size(X);
    d_qual = length(dim_qual);
    phi = [phi, zeros(1, d_lv*d_qual)];
    
    X1 = toLatent(X,dim_qual, z, d_lv, levels);
    x_pred1 = toLatent(x_pred, dim_qual, z,d_lv, levels);
    
    R_xpredx = computeR(phi,x_pred1, X1);
    R_xx = computeR(phi, X1, X1);
    R_xpred = computeR(phi, x_pred1, x_pred1);
    
    mu = (ones(1,n)*(R_xx\ones(n,1)))\(ones(1,n)*(R_xx\y)); 
    sigma = 1/n*(y-mu)'*(R_xx\(y-mu));

    y_pred = mu + R_xpredx*(R_xx\(y-mu));
    y_cov = sigma*(R_xpred - R_xpredx*(R_xx\R_xpredx'));
end
end

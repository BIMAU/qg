function [PC CovMat Var] = pca(G, dim)
%% PC-Analysis. [signals Var PC] = PCA(G). 
% Accepts an {m}x{n} matrix 
% containing {m} variables, measured {n} times.  
%
% It returns the principal components {PC} (combined they
% may perform a change of basis), the transformed data {signals} 
% and the variances {Var} in decreasing order.

[~, experiments] = size(G);
fprintf('Principal component analysis\n'); 
fprintf(' computing covariance matrix...\n'); tic;
mn     = mean(G,2);
GPCA   = G - repmat(mn,1,experiments);  % mean deviation form
CovMat = (GPCA * GPCA')/ (experiments - 1);
Y      = GPCA'/sqrt(experiments - 1);
fprintf(' computing covariance matrix... done (%fs)\n', toc);

fprintf(' computing svd...\n'); tic;
if nargin < 2
    [~, S, PC] = svd(Y);
else 
    [~, S, PC] = svds(Y, dim);
end
fprintf(' computing svd... done (%fs)\n', toc)

S = diag(S);
Var = S.*S;
end
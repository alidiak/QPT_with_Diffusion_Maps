function [U, S] = pca(X)

[m, n] = size(X);

covariancemat = (1/m)*((X).')*X;
[U,S] = svd(covariancemat);


end

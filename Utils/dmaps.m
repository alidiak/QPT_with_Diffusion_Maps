function [E, V, K] = dmaps(X,epsilon,n)
% X is a MxN sample matrix, each row correspond to a sample of N features
 
    M = size(X,1); % M is the number of samples
    Y = zeros(M,M); % Y is the Euclidean distance (squared) matrix
    
    for i = 1 : M
        Y(:,i) = sum(abs(X-X(i,:)).^2,2);
    end
    
    K = exp(-Y/epsilon); % Use the Gaussian kernel to generate the K matrix
 
    P = K./sum(K,2);    % Normalize K such that each row adds up to 1
    
    [V,E]=eigs(P,n); % Find the largest n eigenvalues and eigenvectors of P
    E = diag(E);
    
end
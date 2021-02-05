function [vects, evals, Y, P] = diffusionmaps(X, epsilon, alpha, n_comp)

% Diffusion map method (introduced by Coifman et. al.)

    N = size(X,1);
    N_samples = size(X,2);

    Y = zeros(N_samples,N_samples);  % Y is the Euclidean distance (squared) matrix

    for i = 1 : N_samples
        Y(:,i)= sum(((X(:,i)-X).^2),1)/N;
    end
    
    % Use the Gaussian kernel to generate the K matrix
    K = exp(-Y/epsilon); % Use the Gaussian kernel to generate the K matrix

    q_eps = sum(K,2); % useful normalization term
    
    K_alph = K./((q_eps.^alpha)*(q_eps.^alpha).'); % alpha altered K
    
    d_alph = sum(K_alph,2); % K_alpha normalization term
    
    P = K_alph./sqrt(d_alph*d_alph.'); % finally normalize/construct the transition prob matrix

    [vects, evals] = eigs(P, n_comp);
     
     evals=diag(evals);
    
end
function coef=MPS_psi(psi,s, evals)
% this function takes as input the matrix product state tensor that
% represents the wavefunction psi, along with a given spin
% configuration (s) and returns the MPS calculated coefficient of the 
% wavefunction for that spin configuration. 

% NOTE: this is all specific to how I extract the MPS in get_MPS_psi, the
% form of the tensor generated there is 
% psi= A(1:L, 1:matrix_rows,1:local_hilbert_dim,1:matrix_cols) 

d=max(size(evals));
% d = size(psi{1},2); % this is the local hilbert space dimension, 
% it is the same for each tensor/matrix.

coef=1;
for ii=1:max(size(psi)) % should be =L (lattice size)
    for jj=1:d
        if abs(evals(jj)-s(ii))<1e-4
            slice = jj;  % picks the slice/version of A depending on the spin val
        end
%         if s(ii+1)==evals(jj)
%             slice2=jj;
%         end
    end
  
    A=squeeze(psi{ii}(:,slice,:));
    %A2=squeeze(psi{ii+1}(:,slice2,:));
    if max(size(coef))==1 % fixes a small bug where squeeze transposes A if the singleton
        % dimenions are the first dimensions.
        A=A.'; %%%%%% MAKE SURE IT'S NOT DOING CONJUGATION!!! .' always
    end
    coef=coef*(A);
    
end

end
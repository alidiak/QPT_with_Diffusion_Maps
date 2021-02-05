function R=getR(Z,S)
% The N-spin basis (S1,S2,...SN) to linear basis conversion for spin S
% Each row of Z contains a N-spin basis.
% Z can have arbitrary rows for conversion all at once
% This is consistent with getZ.m
% For spin-1/2 make sure the values in Z are -1/2 & 1/2 instead of -1 & 1.

D=2*S+1; %Hilbert space dimension
N=size(Z,2); %Number of spins
R=ones(size(Z,1),1);
for m=1:N
    R=R+(S-Z(:,N+1-m))*D^(m-1);
end

end
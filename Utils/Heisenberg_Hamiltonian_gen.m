function [H,chi,H_2]=Heisenberg_Hamiltonian_gen(N,adjmat,Z,S, varargin)
%Lanczos diagonalization of H=sum_{<m n>}(X_m*X_n+Y_m*Y_n+Z_m*Z_n)
%Can restrict to spin-K subspace. Runtime and memory roughly scales as DK*N^2
%This is the shortest and easiest code, with slower sparse H generation

AKLT=false;
while ~isempty(varargin)
    switch lower(varargin{1})        
        case 'k'
            K=varargin{2};
        case 'aklt'
            AKLT=varargin{2};
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

% S=1/2; 
D=2*S+1;

if ~exist('Z','var')
     %Spin and Local Hilbert space dimension
    Z=getZ(N,S); %functional call puts literally zero cost on runtime and memory
end

ia=(1:D^N)'; %Linear index array
if exist('K','var')
    chi=ia(sum(Z,2)==K); %chi maps the indices from the projected to original space
else 
    chi=ia; % leave unchanged is condition on total spin is not entered
end
DK=length(chi); %projected space dimension

ia=(1:DK)'; %the new linear index array in the K subspace
Z=Z(chi,:); %the new Z in the K subspace

xi=zeros(D^N,1); xi(chi)=ia; %xi maps the indices from original to projected space

% Now get the off-diagonal Hamiltonian H^{+-}
M=sum(Z(:,1)~=S & Z(:,2)~=-S); %Total number of nonvanishing matrix elements for  S_1^+S_2^-, same for any S_m^+S_n^-

[row,col,val]=find(adjmat);
H=spalloc(DK,DK,N*(N-1)*M);

if AKLT
    H_2=spalloc(DK,DK,N*(N-1)*M);
else
    H_2=0;
end

for m=1:nnz(adjmat)
    ka=Z(:,row(m))~=S & Z(:,col(m))~=-S; %pick nonvanishing ket indices
    ja=xi(chi(ka)+D^(N-col(m))-D^(N-row(m))); %chi(ka) filters out invalid indices cases
    %Only the matrix element <ia(ka)|S_m^+.S_n^-|ja> will be nonzero,and it equals 2*S for S=1/2 and 1    
    Sval=sqrt(S*(S+1)-Z(ja,col(m)).*(Z(ja,col(m))+1)).*...
        sqrt(S*(S+1)-Z(ja,row(m)).*(Z(ja,row(m))-1));
    SpSm=sparse(ia(ka),ja,Sval./2,DK,DK,M); %the sparse matrix (S^+_m.S^-_n)/2
%     SpSm=sparse(ia(ka),ja,S,DK,DK,M); %the sparse matrix (S^+_m.S^-_n)/2
    SzSz=sparse(ia,ia,Z(:,row(m)).*Z(:,col(m)),DK,DK,DK); %the sparse matrix S^z_m.S^z_n
    SS=val(m)*(SpSm+SpSm'+SzSz); %The heisenberg interaction, just SpSm here as S=1/2 (already 1/2 factor)
    H=H+0.5*SS; %Heisenberg matrix (0.5 factor because adjacency matrix sum method doubles interactions)
    if AKLT
        H_2 = H_2+0.5*(SS*SS); % SS squared (for AKLT model) (adjacency mat sums over interactions twice)
    end
end


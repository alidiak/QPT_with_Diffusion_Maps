function matrix=kron_matrix_generator(op,D,N, bc, varargin)
% this function generates the Hamiltonian terms when they consist of a sum
% of local operators. The local operator should be input at op and the end
% of the sum/size of the system should be input as N. (for a local op 
% acting on two sites with open boundary conditions N should be Size-1. The
% op should also be entered as the kron product of the two operator 
% matrices). D is the local Hilbert space size.

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'spatial_var'
            spatial_var=varargin{2};
        case 'op_list'
            op_list=varargin{2};
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

% if exist('op_list',var)
%     for f=2:length(op_list)
%         op=kron(op_list{f-1},op_list{f});
%     end
% else
    op=sparse(op); % make the op sparse
% end

% matrix=zeros([D^N,D^N]);
matrix=spalloc(D^N,D^N, N*D^N);
nops=round(log(size(op,1))/log(D)); % number of sites the entered op is acting on

bc_term=(nops-1); % this is so the boundary conditions add up

% if strcmp(bc,'periodic') && nops>1 % DOESN'T WORK CONSISTENTLY!!!! DON'T USE PBC UNTIL FIXED
%     bc_term=bc_term-1;
%     error('PBC DOES NOT WORK CONSISTENTLY!')
% end

for j = 1:(N-bc_term)  % applies the operator to all sites
   if exist('spatial_var','var')
%        if spatial_var % creating a non even spatial variaiton (for breaking translational symmetry) 
%            w=j;
%        end
        w=spatial_var(j);
   else
       w=1;    
   end
   a=kron(speye(D^(j-1)),op); % this is the first half for a given j
   b=sparse(kron(a,speye(D^(N-j-(nops-1))))); % this is the second half for a given j
   matrix=matrix+w*b;     
end

if strcmp(bc,'periodic')
    a=kron(op_list{2},speye(D^(j-1)));
    b=sparse(kron(a,op_list{1}));
    matrix=matrix+w*b;
end


%    if (N-j-(nops-1))<0 % for conditions with periodic with nops>1
%         a=op; % the form when there is wrap around is always op1*op2*...op_nops*I^(N-nops)
%         b=sparse(kron(a,speye(D^(N-nops))));
%         matrix=matrix+w*b;
%    else
% matrix=2.*matrix;  % it messes up previous code, but this is only half
% the counted interaction terms. (adding the hermitian conjugate is often
% another way to account for the missing terms.)
%matrix=sparse(matrix);

end
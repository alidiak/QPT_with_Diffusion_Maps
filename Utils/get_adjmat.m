function adjmat=get_adjmat(lattice,per,L)
    
% Loading the Graph/Lattice data if not 1d chain
    if strcmp(lattice,'1dchain')
        v=ones((L-1),1);    adjmat = diag(v,1)+diag(v,-1); % adjacency matrix for obc 1dchain
    else
        edge_dir=['~/Documents/QML_Research/Edgemats/']; % Generated in mathematica script
        filedat = importdata([edge_dir '/' lattice per '_edges_' num2str(L) '_vertices.txt'],' ' ,3);
        edges = filedat.data(:,:);
        adjmat = sparse([edges(:,1) edges(:,2)],[edges(:,2) edges(:,1)],[edges(:,3) edges(:,3)],L,L);
    end

end
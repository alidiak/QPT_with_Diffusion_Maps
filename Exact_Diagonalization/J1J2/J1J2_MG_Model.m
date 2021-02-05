%% Init
clear ; close all; clc

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
addpath('../../Utils');

%% Loading adjacency matrix from a file
% Use mathematica notebook "Lattice_adjacency_mat_generator.nb" to generate
% adjacency matrices to impart model on. 

% per='periodic'; % enter '' if non-periodic and 'periodic' if periodic.
% N =20;
% lattice='triangular';
% edge_dir=['../../Edgemats'];
% filedat = importdata([edge_dir '/' lattice per '_edges_' num2str(N) '_vertices.txt'],' ' ,3);
% edges = filedat.data(:,:);
% adjmat = sparse([edges(:,1) edges(:,2)],[edges(:,2) edges(:,1)],[edges(:,3) edges(:,3)],N,N);
% G = graph(adjmat);
% nnfiledat = importdata([edge_dir '/' lattice per 'nextneighbor_edges_' num2str(N) '_vertices.txt'],' ' ,3);
% nnedges = nnfiledat.data(:,:);
% nnadjmat = logical(sparse([nnedges(:,1) nnedges(:,2)],[nnedges(:,2) nnedges(:,1)],...
%     [nnedges(:,3) nnedges(:,3)],N,N));

%% Creating an adjacency matrix for a simple 1d lattice

N=10; % choose system size 
save_dat=true;

v=ones((N-1),1); v2=ones((N-2),1);
adjmat = diag(v,1)+diag(v,-1);
nnadjmat=diag(v2,2)+diag(v2,-2);
%v(1:2:end)=-v(1:2:end); op_adjmat=diag(v,1)+diag(v,-1);

% periodic conditions
per='periodic';
lattice='1dchain';

if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
%     op_adjmat(1,end)=(-1)^N; op_adjmat(end,1)=(-1)^N;
    nnadjmat(1,end-1)=1;nnadjmat(end-1,1)=1;
    nnadjmat(2,end)=1;nnadjmat(end,2)=1;
end

G=graph(adjmat);% p=plot(G);

time_evol=false; dt=0.1; tmax=1;

degen_per_h=true; MH=true;

%% Exact Diagonalization Matrices

spin=1/2; Z=getZ(N,spin); Dim=2*spin+1; %Get the matrix elements for all {Z_m} which are diagonal

KK=0;
% This is the nearest neighbor heisenberg interaction Hamiltonian
options = {'K',KK};
[H1,chi]=Heisenberg_Hamiltonian_gen(N, adjmat,Z,spin,options{:}); % only enter K=0 if N is even
H2=Heisenberg_Hamiltonian_gen(N, nnadjmat,Z,spin,options{:}); 

DK=size(H1,1); % the reduced matrix size
D=Dim^N;
ia=(1:D)'; % Linear index array
ja= (1:DK)';
xi=zeros(Dim,1); xi(chi)=ja; %xi maps the indices from original to projected space
Zred=Z(chi,:); % Z in the reduced subspace

%Now constructing the total Sx,Sy,Sz operator (we can easily square them without constructing two-point correlators)
% Sx=spalloc(D,D,D);
% Sy=spalloc(D,D,D);
Sz=sparse(ja,ja,sum(Zred,2),DK,DK,DK); % should be all 0 sparse if K=0

% for m=1:N
%     %jj=chi+2*Zred(:,m)*2^(N-m);
%     %Sx=Sx+sparse(1:size(nonzeros(xi(jj))),nonzeros(xi(jj)),1,DK,DK,DK);
%     %Sy=Sy+sparse(1:size(nonzeros(xi(jj))),nonzeros(xi(jj)),1i*Zred(:,m),DK,DK,DK);
%     %sum(xi(jj)) % never seems to have any nonzero elements in the matrix
%     
%     kk=ia+2*Z(:,m)*2^(N-m); % normal way
%     Sx=Sx+sparse(ia,kk,1,D,D,D);
%     Sy=Sy+sparse(ia,kk,1i*Z(:,m),D,D,D);
% end

%% Option to add a symmetry breaking term to the dimer ground-state 
% this is the BAHC Hamiltonian which has the even or odd dimer as it's
% ground-state. Adding this to the J1-J2 Hamiltonian directly breaks the
% even/odd dimer ground-state symmetry.

J_symm_break=0.0025;

d=-1; % (-1) for an even dimer GS and (1) for an odd dimer GS
%  (breaks symmetry in favor toward each dimer respectively (if added term is positive))

v_b=ones((N-1),1); v_b(1:2:end)=-v_b(1:2:end);
k=ones((N-1),1)-d.*(v_b); adjmat_BAHC = diag(k,1)+diag(k,-1);
if strcmp(per,'periodic')
    adjmat_BAHC(1,end)=1-d*(-1)^N; adjmat_BAHC(end,1)=1-d*(-1)^N;
end

options={'K',0}; % K is the reduced subspace of Sz
[H_BAHC]=Heisenberg_Hamiltonian_gen(N, adjmat_BAHC,Z,spin,options{:});

%% Finding Energies and Creating Samples

J1= 1;  % nearest neighbor
J2= 0.75; % next nearest neighbor 1/2 is the MG model specifically

dJ = 0.05;
Jlist = 0.0:dJ:1.0;
% Jlist= [0:dJ:0.40 0.45:0.005:0.55 0.6:dJ:1.0];
% Jlist=0.45:0.01:0.55;
% Jlist=[0.5];
sample_size=10;

% creating an Sz term that increases in strength with lattice site. Should
% break translational symmetry.
% options={'spatial_var',true}; sz=[[1/2,0];[0,-1/2]];  
% spatial_vary_sz=kron_matrix_generator(sz,2,N,per,options{:}); 

% dimer_corr_r=zeros([floor(N/4),length(Jlist)]);
% for r=2:2:floor(N/2)
% creating the dimer correlation order parameter as an operator
if N<=20 && N>4
r=2;
l0=round((N+mod(r,2))/2); r2=2*round(r/2);
sx=0.5*[[0,1];[1,0]]; % op=sz; 
t1=l0-r2/2; t2=l0+r2/2;
a=kron(speye(Dim^(t1-1)),kron(sx,sx)); % first term corresponding to I*S^x_l0-r/2 S^x_l0-r/2+1
b=kron(a,speye(Dim^((t2-1)-(t1+1)))); % sandwich I term (between l0-r/2+2 and l0+r/2)
c=kron(b,kron(sx,sx)); % term S^x_l0+r/2 S^x_l0+r/2+1
pos_corr=kron(c,speye(Dim^(N-(t2+1)))); % I_N-l0+r/2+1 fills rest of the lattice points/hilb space
% second term t2 goes to previoud t2+1 (l0+r/2->l0+r/2+1)
t2=t2+1; 
b=kron(a,speye(Dim^((t2-1)-(t1+1)))); % sandwich I term (between l0-r/2+2 and l0+r/2)
c=kron(b,kron(sx,sx)); % term S^x_l0+r/2 S^x_l0+r/2+1
neg_corr=kron(c,speye(Dim^(N-(t2+1))));

% combining the correlators to make the order parameter.
Sx_dimer_corr=(1/spin^4)*(pos_corr-neg_corr);
else
    Sx_dimer_corr=0;
end

sz=[[1,0];[0,-1]]; szsz=kron(sz,sz);
l2=10; % second lattice site
sz1=kron(sz,speye(2^(l2-2))); sz2= kron(sz1,sz); 
% sz1=kron(speye(2^(N-2)),szsz); %a=kron(speye(2^(l2-1)),sz); sz2=kron(a,speye(2^(N-l2)));
Sz_spin_corr=kron(sz2,speye(2^(N-l2)));
% Sz_spin_corr=sz1;%+sz2;

% this is the model order parameter (changing the adjmat as it is
% alternating/staggered by lattice site in the order param
% HZZ=sparse(ia,ia,sum(((Z*op_adjmat).*Z),2),2^N,2^N,2^N);

samplepos = zeros(max(size(Jlist))*sample_size,1); % will have sample_size

E_gap = zeros(max(size(Jlist)),1);

% orderpar=zeros(1,max(size(Jlist)));;
lattice_Sz = zeros(N,max(size(Jlist)));
lattice_Sx = zeros(N,max(size(Jlist)));
sigmax = zeros(1,max(size(Jlist)));
groundstates = zeros(size(H1,1),max(size(Jlist)));
energy=zeros(1,max(size(Jlist)));
dimer_corr=zeros(1,length(Jlist));
spin_corr=zeros(1,length(Jlist));

for n=1:(max(size(Jlist)))
    
    %%%%%% CHOOSE HERE WHICH COUPLING J TO VARY! %%%%
    J2 = Jlist(n);
    
    if exist('H_BAHC','var')
        H=(J1*(H1)+J2*(H2))+J_symm_break*H_BAHC;
    else
        H=(J1*(H1)+J2*(H2)); % 2 for difference in evals between Carleo's codes 
    end
%     H=J2*(H2+H1);
    
    % For the ks testing (breaks translational symm)
%     J2=0.5; Bz=Jlist(n);
%     H=2*(J1*(H1)+J2*(H2))+Bz*spatial_vary_sz; 

    nevals = 2;
    [evecs,evals] = eigs(H, nevals, 'SA');
    fprintf(['Min. Energy for J2/J1=' num2str(J2/J1) ': ' num2str(evals(1,1)) '\n'])
    
    energy(n)=evals(1,1);
    
    % energy gap calc
    % uniq=uniquetol(diag(evals));
    diagvals = diag(evals);
    E_gap(n)= diagvals(2)-diagvals(1);
    
    GS=evecs(:,1);
    groundstates(:,n)=GS;
    GS_probs = abs(GS.*1).^2; % probabilities of each state
    
    lattice_Sz(:,n)=((abs(GS.*1).^2).')*Zred; 
    
%     sigmax(n) = (GS')*Sx(chi,chi)*GS; % average sigma X
    
%     dimer_corr(n)=GS'*Sx_dimer_corr(chi,chi)*GS;
    spin_corr(n)=GS'*Sz_spin_corr(chi,chi)*GS;
    %     orderpar(n)=(GS')*HZZ*GS;
    
    if J2==0.5
        dimer_probs=GS_probs>1e-8; % gives binary array
        dimer_Z=2*Zred(dimer_probs,:);
    end

    % Adding the probabilities so that each state can be represented by a
    % unique probability between 0 and 1. 
    for m=2:DK
        GS_probs(m) = GS_probs(m)+GS_probs(m-1);
    end

    % Sampling %
    for ii=1:sample_size
        a = rand;
        samplepos((n-1)*sample_size+ii,1) = sum(a>=GS_probs)+1;
        % finds the pos where a is greater than the lowest prob. This is the 
        % sample we choose by picking a random number a.
    end
end
% Plotting
fprintf('\n\n SigmaZ GS Expectation values \n')
fprintf('J2/J1= %.2f : %.3f \n', [(Jlist./J1)',sum(lattice_Sz,1)']')

savedir=['plots/' lattice per '/N_' num2str(N)];
mkdir(savedir)

% figure, plot(Jlist./J1,energy,'-o'), title('Energy vs. Ratio J2/J1')
% figname=[savedir '/Energy_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% 
% figure, plot(Jlist./J1,E_gap,'-o'), title('Energy Gap vs. Ratio J2/J1')
% figname=[savedir '/Energy_Gap_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% 
% figure, plot(Jlist./J1,dimer_corr,'-o'), title(['Dimer Correlation vs. Ratio J2/J1 for r=' num2str(r)])
% figname=[savedir '/Dimer_Corr_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% dimer_corr_r(r,:)=dimer_corr;
% end

setfig(16,'J_2',0,['<s_1^z s_{' num2str(l2) '}^z>'],0,'Mid-Lattice\ Spin\ Correlation','on')
plot(Jlist,spin_corr,'o-')

if exist('dimer_corr_r','var')
    figure; skip=2;cmap=colormap(parula(length(Jlist/skip)));
    for k=1:skip:length(Jlist)
        plot(log10(2:2:floor(N/2)),log10(dimer_corr_r(:,k)),'Color',cmap(k,:),...
            'DisplayName',['J2=' num2str(Jlist(k))])
        hold on
    end
end
% legend('FontSize',8); xlabel('Log_10(r)');ylabel('Log_10(C^x)'); title('Dimer Correlation vs r and J_2')
% figname=['plots/' lattice per '/N_' num2str(N) '/Dimer_corr_vs_r_and_J2_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% figure, plot(Jlist./J1,orderpar/(N*spin^2),'-o'), title('Order Param vs. Ratio J2/J1')
% ylabel('1/(LS^2) \Sigma_l (-1)^l \sigma_l^z \sigma_{l+1}^z'); xlabel('J2/J1');
% figname=['plots/' lattice per '/N_' num2str(N) '/Orderparam_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% 
% figure, plot(Jlist./J1,sum(lattice_Sz,1)','-o'), title('Avg SigmaZ vs. Ratio J2/J1')
% figname=['plots/' lattice per '/N_' num2str(N) '/SigmaZ_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

samples=Zred(samplepos,:); 

%save(['plots/1dchainperiodic/Finite_size/J1J2_spinhalf_N_' num2str(N) '_' ...
  %  num2str(J_symm_break) 'symm_breaking.mat']); 

%% loading previous run
% path = matlab.desktop.editor.getActiveFilename;
% cd(fileparts(path));
% N=24;
% load(['plots/1dchainperiodic/Finite_size/J1J2_spinhalf_N_' num2str(N) '.mat']); 
% load(['plots/1dchainperiodic/Finite_size/J1J2_spinhalf_N_24_0.0025symm_breaking.mat'])
% load(['BAHC_N' num2str(N) '_even_dimer_samples.mat']); % Adding a dimer state to the mix

% Comment out if just even or odd are desired
% odd_dimer_samples=dimer_samples;
% load(['BAHC_N' num2str(N) '_even_dimer_samples.mat']);
% even_dimer_samples=dimer_samples;
% dimer_samples=[even_dimer_samples(1:0.5*sample_size,:) ; ...
%     odd_dimer_samples((0.5*sample_size+1):end,:)];

% if exist('dimer_samples','var')
%     samples=[samples ; dimer_samples]; Jlist=[Jlist -0.1];
% end
% 
% sample_size=500; samplepos = zeros(max(size(Jlist))*sample_size,1); 
% for n=1:(max(size(Jlist)))
%    GS_probs = abs(groundstates(:,n)).^2; % probabilities of each state
%    for m=2:DK
%         GS_probs(m) = GS_probs(m)+GS_probs(m-1);
%     end
%     for ii=1:sample_size
%         a = rand;
%         samplepos((n-1)*sample_size+ii,1) = sum(a>=GS_probs)+1;
%     end
% end
% if spin==0.5
%     samples=2*Zred(samplepos,:); 
% else
%     samples=Zred(samplepos,:);
% end
% 
% % Number of states equal to the connection states /antiferromagnetic ground-states
% pt=-0.1;
% alt=0.5*ones([1,N]); alt(1:2:end)=-alt(1:2:end);
%   Nconnectors = sum(sum(alt==samples(((find(Jlist==pt)-1)*sample_size+1):(find(Jlist==pt))*sample_size,:),2)==N)...
%     +sum(sum(-alt==samples(((find(Jlist==pt)-1)*sample_size+1):(find(Jlist==pt))*sample_size,:),2)==N)
% connectors_loc=[find(sum(alt==samples(((find(Jlist==pt)-1)*sample_size+1):(find(Jlist==pt))*sample_size,:),2)==N) ; ...
%     find(sum(-alt==samples(((find(Jlist==pt)-1)*sample_size+1):(find(Jlist==pt))*sample_size,:),2)==N) ]; % where they occur

%% Other clustering methods
eps=2.85; minpts=3; % minpts rule of thumb, >= dim+1 or 2*dim where dim is dimension of data. min logical is 3
idx=dbscan(samples,eps,minpts);

num_clusters=zeros([length(Jlist),1]); N_outlier=zeros( [length(Jlist),1]);
for h=1:length(Jlist)
    num_clusters(h)=length(unique(idx((h-1)*sample_size+1:h*sample_size))); 
    N_outlier(h)=sum(idx((h-1)*sample_size+1:h*sample_size)==-1);
end
figure; yyaxis left; plot(Jlist,num_clusters,'.-'); hold on
% yyaxis right; plot(Jlist,N_outlier,'.-')
title(['DBSCAN predicted # of Clusters \epsilon=' num2str(eps)])
% legend({'# Clusters','# Outliers'})

%% DBSCAN over a range
eps_l=linspace(2.65,3.2,30); minpts=10;

tic
num_clusters=zeros([length(Jlist),length(eps_l)]); N_outlier=zeros( [length(Jlist),length(eps_l)]);
for ee=1:length(eps_l)
    idx=dbscan(samples,eps_l(ee),minpts);
    
    for h=1:length(Jlist)
        num_clusters(h,ee)=length(unique(idx((h-1)*sample_size+1:h*sample_size))); 
        N_outlier(h,ee)=sum(idx((h-1)*sample_size+1:h*sample_size)==-1);
    end
end
toc 

%% plot DBSCAN vs epsilon

plot_all_pt=0.49; plot_all_end=0.515; skip=6;     
kk=1;  lines=[]; 
setfig(16,'DBSCAN\ \epsilon',0,'Number\ of\ Clusters',0,['Number\ of\ Clusters\ vs\ radius\ minpts=' num2str(minpts)],'on');
while kk<=length(Jlist)
    kk 
    if kk==find(abs(Jlist-plot_all_pt)<1e-4)
        skprev=skip; skip=1;
    end
    if kk==find(abs(Jlist-plot_all_end)<1e-4)
        skip=skprev;
    end
    if kk==1 || mod(kk,skip)==0
        lines(end+1)=plot((eps_l(:)), (num_clusters(kk,:)),'.-','MarkerSize',10,'DisplayName',['J_2=' num2str(Jlist(kk))]);
        hold on;
    end
     kk=kk+1;
end
c_map=parula(length(lines)); %c_map(find(abs(v_list-0.5)<1e-4)/skip,:)=[1,0,0]; % makes J_2=0.5 red
c_map(ceil((length(lines))/2),:)=[1,0,0]; % assuming J_2=0.5 is in the center
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/DBSCAN_Nclusters_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Calculating the order parameter with the samples

% Czlist=[];
% for n=1:length(Jlist)
%     S=samples(((n-1)*sample_size+1):(n*sample_size),:);
%     for r=2:2:(N/2)
%         l0=round((N+mod(r,2))/2); %r=2*round(r/2);
%         if  (N-(l0+r/2+2))<0
%             break
%         else   
%             Cz=mean((1/(spin^4))*S(:,l0-r/2).*S(:,l0-r/2+1).*( ...
%                 S(:,l0+r/2).*S(:,l0+r/2+1) -S(:,l0+r/2+1).*S(:,l0+r/2+2)));
%             Czlist(end+1)=Cz;
%         end
%     end
% end
% rmax=length(Czlist)/length(Jlist);
% Czlist=reshape(Czlist,[rmax,length(Jlist)]);
% 
% figure; p=plot(log10(2:2:(N/2)),log10(abs(Czlist(:,2:2:end))),'.-'); 
% set(p,{'color'},num2cell(jet(length(p)),2))
% title('C_{dimer}'); legend(cellstr(num2str(Jlist(2:2:end),'J=%f'))); xlabel('log10(r)');
% 
% r=50;l0=round((N+mod(r,2))/2); r=2*round(r/2);
% Cz=(1/(spin^4))*samples(:,l0-r/2).*samples(:,l0-r/2+1).*( ...
%     samples(:,l0+r/2).*samples(:,l0+r/2+1) -samples(:,l0+r/2+1).*samples(:,l0+r/2+2));
% Czlist=zeros([length(Jlist),1]);
% for n=1:length(Jlist)
%     Czlist(n)=mean(Cz(((n-1)*sample_size+1):(n*sample_size)));
% end
% figure; plot(Jlist,abs(Czlist),'k-*'), title(['Cz dimer correlation vs J_2 for r=' num2str(r)])
% figname=[savedir '/Sample_Dimer_Corr_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
%% Testing probability comparison methods on symmetry transformations

xlate_by=1; options={'statistic','Cramer','plot','on'};
xlated_samples=[samples(:,end-(xlate_by-1):end)  samples(:,1:end-xlate_by)];
[ks_tests,ks_stats]=symm_kstesting(getR(samples/2,1/2),getR((xlated_samples)/2,1/2),Jlist,1/2,...
    ['translation\ by\ ' num2str(xlate_by)], 'J_2',options{:}); %'Bz\ of\ \Sigma_i i\sigma_i^z ',options{:});

%% Original Distance measure (Euclidian)
% eps_list=0.0025:0.0025:0.2; % normal method seems to work better at lower eps
degen_tolerance = 10^(-2.5);
% eps_list=logspace(-2.2,-1.5,20); 
eps_list=10^(-2.3);

% eps_list=logspace(-0.8,-0.2,20); % if unnormalized by L in Diff maps (more general)
options = {'degen_tol',degen_tolerance,'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',1,...
'fig_title', 'degenerate\ evals=1\ vs\ Symm\ Broken\ J_2\','var_name','J_2'};
[ndegenerate,secondeval,eval_avg,Ys]=diffmap_list(0.5.*samples,Jlist,eps_list,options{:});

%% Checking the dispersion/standard deviation of the distance matrix
if exist('Ys','var')

 % Histogramming the distances
distances=zeros([N/2+1,length(Jlist)]); i=1; dists=(0:2:N)/N;
for k=0:2:N
    distances(i,:)=squeeze(sum(sum(abs(N*Ys-k)<1e-4)));
    i=i+1;
end

std_dev=zeros([length(Jlist),1]);
for kk=1:length(Jlist)
%     std_dev(kk)=std(distances(:,kk)/sum(distances(:,kk)));
    std_dev(kk)=pdf_std(distances(:,kk)/sum(distances(:,kk)), dists.'); % have to do the pdf expect value when using histograms
% %     std_dev(kk)=std(nonzeros(d_hist)/sum(nonzeros(d_hist)));
%     y=reshape(Ys(:,:,kk),[sample_size*sample_size,1]);
%     pd=fitdist(y,'Normal'); 
%     std_dev(kk)=pd.sigma;
%     std_dev(kk)= sqrt(mean(y.^2)-mean(y)^2);
% %     std_dev(kk)=std2(y);
%     std_dev(kk)= sqrt(mean((y-mean(y)).^2));
%     std_dev(kk)=std(Ys(:,:,kk),1,'all');
end

figure=setfig(16,'J_2',0,'Standard\ Deviation',0,'Distance\ Distributions\ Standard\ Deviation','on'); 
plot(Jlist,std_dev,'.-')
figname=[savedir '/Distance_distribution_J1J2_std_dev_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
end

%% Testing features of the distance matrix

avg_dist=squeeze(mean(mean(Ys,1),2)); max_dist=squeeze(max(max(Ys)));
num_max=squeeze(sum(sum(0.5==Ys,1),2)); num_zero=squeeze(sum(sum(0.5==Ys,1),2));

tol=0.11; num_close_max=squeeze(sum(sum((Ys-(0.5-tol)>0),1),2)); 
figure; plot(Jlist,num_close_max,'-*'); 
title(['number distances near max with tol ' num2str(tol)]) 
figname=[savedir '/number_dist_near_max_' num2str(tol) '_' per '_' lattice '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% number unique
nunique=zeros([length(Jlist),1]);
for n=1:length(Jlist)
    nunique(n)=length(unique(samplepos(((n-1)*sample_size+1):(n*sample_size))));
end
figure('Visible','on'); plot(Jlist,nunique,'k-*'), ylabel('# unique samples'), xlabel('J_2') 
title('unique samples per J_2'); 
figname=[savedir '/number_unique_samples_' per '_' lattice '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% if save_dat
%     save(['plots/1dchainperiodic/Finite_size/J1J2_spinhalf_N_' num2str(N) '.mat']);
% end
%% Exploring log scaling of the 1st minus 2nd eval at J2=0.5 vs epsilon

setfig(16,'Log_{10}(\epsilon)',0,'Log_{10}(1-\lambda_1)',0,'Exponential\ scaling\ of\ 1-\lambda_1','on');
lines=[]; skip=4; 
for kk=2:skip:length(Jlist)
lines(end+1)=plot(log10(eps_list), log10(1-secondeval(:,kk)),'.-','DisplayName',['J_2=' num2str(Jlist(kk))]);
% xlabel('Log10(\epsilon)'); ylabel('Log10(1-\lambda_1)');
% title('Exponential scaling of 1-\lambda_1')
hold on
end
c_map=parula(length(lines)); c_map((find(Jlist==0.5))/skip,:)=[1,0,0]; % makes J_2=0.5 red
if exist('dimer_samples','var')
    c_map(find(Jlist==-0.1)/skip,:)=rgb('black'); % makes the added dimer samples black
end
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/1minuslambda2_exponentialscaling_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

plot_all_pt=0.49; plot_all_end=0.515; skip=4;     
kk=1;  lines=[]; 
setfig(16,'Log_{10}(\epsilon)',0,'Degeneracy',0,'Exponential\ Scaling\ of\ Degeneracy','on');
while kk<=length(Jlist)
    kk 
    if kk==find(abs(Jlist-plot_all_pt)<1e-4)
        skprev=skip; skip=1;
    end
    if kk==find(abs(Jlist-plot_all_end)<1e-4)
        skip=skprev;
    end
    if kk==1 || mod(kk,skip)==0
        lines(end+1)=plot((eps_list(:)), (ndegenerate(:,kk)),'.-','MarkerSize',10,'DisplayName',['J_2=' num2str(Jlist(kk))]);
        hold on;
    end
     kk=kk+1;
end
c_map=parula(length(lines)); %c_map(find(abs(v_list-0.5)<1e-4)/skip,:)=[1,0,0]; % makes J_2=0.5 red
c_map(ceil((length(lines))/2),:)=[1,0,0]; % assuming J_2=0.5 is in the center
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/ndegen_exponentialscaling_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'Log_{10}(\epsilon)',0,'Degeneracy',0,'Exponential\ Scaling\ of\ Degeneracy','on');
lines=[];
for kk=1:skip:length(Jlist)
lines(end+1)=plot(log10(eps_list), (ndegenerate(:,kk)),'.-','DisplayName',['J_2=' num2str(Jlist(kk))]);
% xlabel('Log10(\epsilon)'); ylabel('Log10(1-\lambda_1)');
% title('Exponential scaling of 1-\lambda_1')
hold on
end
c_map=parula(length(lines)); c_map((find(abs(Jlist-0.30)<1e-4)+3)/skip,:)=[1,0,0]; 
set(lines,{'color'},num2cell(c_map,2)); 
legend('FontSize',8)
figname=[savedir '/ndegen_exponentialscaling_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Finite sized scaling processing

cd '~/Documents/QML_Research/Exact_Diagonalization/J1J2'
dir= ['plots/1dchainperiodic/Finite_size/'];
file_list=glob([dir '**.mat']);
% fig3=setfig(16,'h',0,'Avg\ Samples',0,'Number\ of\ Unique\ Samples\ Avg','on');
N_list=zeros([length(file_list),1]); connectors_n=zeros([length(file_list),1]);
num_exact_point_min=0;
for ff = 1:length(file_list)
    file=file_list{ff}
    load(file); % load the data    
    N_list(ff)=N;

    if ff==1
        ndegenerate_n=zeros([size(ndegenerate,1),length(file_list)]);
        secondeval_n=zeros([size(secondeval,1),length(file_list)]);
%         second_eval_minJ2_diff=zeros([size(secondeval,1),length(file_list)]);
        second_eval_minJ2_diff=zeros([length(file_list),1]);
    end
    
    ndegenerate_n(:,ff)=ndegenerate(:,find(Jlist==0.5));
    secondeval_n(:,ff) = secondeval(:,find(Jlist==0.5));
    
    % finding the minimum J_2 besides (J_2=0.5)
    eps_min= find(abs(eps_list-10^(-2))<=1e-3); % regime where 0.5 stands out as minimum.
    sorted_J2_mins=sort(log10(1-secondeval(eps_min,:)));
    min_J2=find(sorted_J2_mins(1)==log10(1-secondeval(eps_min,:)));
    if min_J2== find(Jlist==0.5)
        num_exact_point_min=num_exact_point_min+1;
        min_J2=find(sorted_J2_mins(2)==log10(1-secondeval(eps_min,:)));
    end
    
    second_eval_minJ2_diff(ff)= log10(1-secondeval(eps_min,min_J2))-log10(1-secondeval(eps_min,find(Jlist==0.5)));
    
    % number of connecting antiferromagnetic GS states
    Nconnectors = sum(sum(alt==samples(((find(Jlist==0.5)-1)*sample_size+1):(find(Jlist==0.5))*sample_size,:),2)==N)...
    +sum(sum(-alt==samples(((find(Jlist==0.5)-1)*sample_size+1):(find(Jlist==0.5))*sample_size,:),2)==N);
%     connectors_loc=[find(sum(alt==samples((find(Jlist==0.5)*sample_size):(find(Jlist==0.5)+1)*sample_size,:),2)==N);
    connectors_n(ff)=Nconnectors;

    alt=0.5*ones([1,N]); alt(1:2:end)=-alt(1:2:end);

end

[N_list,p]=sort(N_list); secondeval_n=secondeval_n(:,p);ndegenerate_n=ndegenerate_n(:,p);
second_eval_minJ2_diff = second_eval_minJ2_diff(p);
connectors_n=connectors_n(p);

setfig(16,'N',0,'Log_{10}(1-\lambda_1)',0,'Finite\ size\ scaling\ of\ 1-\lambda_1','on');
for kk=1:3:(size(secondeval_n,1)-6)
    plot(N_list,log10(1-secondeval_n(kk,:)),'-*',...
        'DisplayName',['\epsilon=' num2str(eps_list(kk))]); hold on
end
legend;figname=[dir '/Finite_size_Scaling_secondeval'];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'N',0,'Number\ Connecting\ States',0,'Finite\ size\ scaling\ connecting\ states','on');
plot(N_list,connectors_n,'-*'); figname=[dir '/Connected_states_vs_N'];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'N',0,'Log_{10}(\lambda_{min}-\lambda_{J_2=0.5})',0,'Finite\ size\ scaling\ of\ min\ 1-\lambda_1\ to\ 1-\lambda_1\ at\ J=0.5\ difference','on');
for kk=1:3:(size(second_eval_minJ2_diff,1)-6)
    plot(N_list,(second_eval_minJ2_diff(:)),'-*',...
        'DisplayName',['\epsilon=' num2str(eps_list(kk))]); hold on
end
legend;figname=[dir '/Finite_size_Scaling_mineval_minus_J2eval'];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Exploration of the dimer state diffusion map behavior as a function of epsilon

N=12; binary_dimer=permn([0,1],N); % getting full distribution/# of possible states
% eps_list=logspace(-0.65,-0.2,40); 
% eps_list=logspace(-2.05,-1.5,40); 
eps_list=10^(-1.5);
degen_tol=10^(-2.5);
sample_size=1000;
% clear sample_size

% creating the actual dimer state
even_dim=zeros([2^N,2*N]); odd_dim=zeros([2^N,2*N]);
for jj=1:N
    s=binary_dimer(:,jj);
    dimer=s.*[0.5,-0.5]+abs(s-1).*[-0.5,0.5]; % 1 states are up, down - 0 states are down, up
    odd_dim(:,2*jj-1:(2*jj)) = dimer; 
    if jj==N
        even_dim(:,1)=dimer(:,end); even_dim(:,2*jj)=dimer(:,1);
    else
        even_dim(:,2*jj:2*jj+1) = dimer;
    end
end
% distance matrix of dimer like this and dist matrix of binary_dimer are
% equal within a factor of 2. (divide odd/even dimer by 2 and its equal).

% Distances in odd dimer vs even,  cross-terms.
xlated_dim=[odd_dim(:,end) odd_dim(:,1:end-1)];
cross_terms=zeros([size(odd_dim,1),size(odd_dim,1)]);
for kk=1:size(odd_dim,1)
    cross_terms(kk,:)=sum((xlated_dim-odd_dim(kk,:)).^2,2);
end
terms=[];
for k=0:2:2*N
    terms(end+1)=sum(sum(abs(cross_terms-k)<1e-4));
end
setfig(16,'Unnormalized\ Distance',0,'Number\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Single\ vs\ Combined\ Dimer\ Distance\ Cross-Terms','on');
plot(0:2:2*N,2.*terms/4,'o-');

combined_dim= [odd_dim ; even_dim];

alt=1*ones([1,2*N]); alt(1:2:end)=1-alt(1:2:end);
connectors_loc= sort([find(sum((alt==combined_dim),2)==2*N) ; ...
    find(sum((abs(1-alt)==combined_dim),2)==2*N)]);

% combined_dim_mod=zeros(size(combined_dim));
% for ii=length(connectors_loc)-1
% combined_dim_mod=[combined_dim_mod combined_dim(connectors_loc(ii)+1:(connectors_loc(ii+1)-1),:)];
% end

dim_degen = zeros([length(eps_list),1]); comb_degen=zeros([length(eps_list),1]);
for ee=1:length(eps_list)

    if   exist('sample_size','var') %&& (2^N)>=sample_size
        r1=randi(2^N,[sample_size,1]); r2=randi(2^N,[sample_size/2,1]);   
        r3=randi(2^N,[sample_size/2,1]);   
    else
        r1=1:2^N; r2=r1; r3=r1;
    end
    
    alpha=1;
    [~,vals,Y]= diffusionmaps(odd_dim(r1,:).',eps_list(ee),alpha,10);
    dim_degen(ee)=sum(abs(ones(size(vals))-vals)<degen_tol);
%     dim_degen(ee)= sum((1-vals).^2);

    [~,vals2,Y_comb]= diffusionmaps([odd_dim(r2,:) ; even_dim(r3,:)].',eps_list(ee),1,10);
    comb_degen(ee)=sum(abs(ones(size(vals2))-vals2)<degen_tol);
%     comb_degen(ee)= sum((1-vals2).^2);

end

setfig(16,'Log_{10}(\epsilon)',0,'Degeneracy',0,'Exponential\ Scaling\ of\ Degeneracy','on');
plot(log10(eps_list),log10(comb_degen),'.-'); hold on; plot(log10(eps_list),log10(dim_degen),'.-','MarkerSize',10);
% plot(1./eps_list,log10(comb_degen-2),'.-'); hold on; plot(1./eps_list,log10(2.*dim_degen),'.-');
legend({'Combined Dimers','Single Dimer'})
figname=['Toy_model/dimer_vs_combined_dimer_N' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% comparing the number of distances equal to n singlet flips
Analytic=[]; Observed=[]; Analytic_comb=[]; Observed_comb=[];
for k=0:N
    Analytic(end+1)=nchoosek(N,k)*2^N; 
    Analytic_comb(end+1)=2*nchoosek(N,k)*2^N+4*nchoosek(2*N,2*k);
    Observed(end+1)=sum(sum(abs(2*N*Y-k)<1e-4));
    Observed_comb(end+1)=sum(sum(abs(2*N*Y_comb-k)<1e-4)); 
%     Analytic(end+1)=nchoosek(N,k)/2^N; % normalized versions
%     Analytic_comb(end+1)=(2*nchoosek(N,k)*2^N+4*nchoosek(2*N,2*k))/(2^(2*N+2));
    % first term for each dimer state, second for 'cross terms'
end
setfig(16,'Unnormalized\ Distance',0,'Number\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Analytic\ vs\ Observed\ Distance\ Distribution','on');
plot(0:2:2*N,Analytic,'o-'); hold on; plot(0:2:2*N,Observed,'*-');
legend({'Analytic Distribution','Observed Distribution'})
figname=['Toy_model/Analytic_vs_Observed_dist' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'Unnormalized\ Distance',0,'Number\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Analytic\ vs\ Observed\ Distance\ Distribution','on');
plot(0:2:2*N,Analytic_comb,'o-'); hold on; plot(0:2:2*N,Observed_comb,'*-');
legend({'Analytic Distribution','Observed Distribution'})
figname=['Toy_model/Combined_Analytic_vs_Observed_dist' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'Unnormalized\ Distance',0,'Normalized\ Number\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Single\ vs\ Combined\ Dimer\ Distance\ Distributions','on');
plot(0:2:2*N,Observed./sum(Observed),'o-'); hold on; plot(0:2:2*N,Observed_comb./sum(Observed_comb),'*-');
legend({'Single Dimer','Combined Dimer'})
figname=['Toy_model/Normalized_Observed_Combined_vs_dimer_dist' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'Unnormalized\ Distance',0,'Normalized\ Number\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Single\ vs\ Combined\ Dimer\ Distance\ Distributions\ Analytic','on');
plot(0:2:2*N,Analytic,'o-'); hold on; plot(0:2:2*N,Analytic_comb,'*-');
legend({'Single Dimer','Combined Dimer'})
figname=['Toy_model/Normalized_Observed_Combined_vs_dimer_dist' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
%% Plotting the analytical curves for high N and as a function of N
% n_list=6:40; rel_dif=[];
% for N=n_list
N=40; 
Analytic=[]; Analytic_comb=[]; 
for k=0:N
    Analytic(end+1)=nchoosek(N,k)/2^N;
    Analytic_comb(end+1)=(2*nchoosek(N,k)*2^N+4*nchoosek(2*N,2*k))/(2^(2*N+2));
end

dists=0:2:2*N; dists=dists/(2*N); % account for normalization
single_sig=pdf_std(Analytic, dists); % function simply performs a prob based expectation value
comb_sig=pdf_std(Analytic_comb, dists); % and sqrt( (y.^2)-y.^2)
rel_dif(end+1)=(single_sig-comb_sig)/comb_sig;

setfig(16,'Distance\ d',0,'Number\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Analytical\ Single\ vs\ Combined\ Dimer\ Distance\ Distributions','on');
plot(dists,Analytic,'o-'); hold on; plot(dists,Analytic_comb,'*-'); hold on;
legend({['Single Dimer, std=' num2str(single_sig)],['Combined Dimer, std=' num2str(comb_sig)]})
figname=['Toy_model/Normalized_Analytic_Combined_vs_dimer_dist' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

 setfig(16,'Distance\ \frac{d}{L}',0,'Ratio\ of\ Matrix\ Elements\ Equal\ to\ Distance',...
    0,'Analytical\ Single\ Over\ Combined\ Dimer\ Distance\ Distributions','on');
plot((0:2:2*N)./(2*N),Analytic./Analytic_comb,'.-','MarkerSize',10); 
figname=['Toy_model/Normalized_Dimer_Ratio' num2str(2*N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% end
% setfig(16,'Distance\ d',0,'Standard\ Deviation\ Relative\ Difference\',...
%     0,'Analytical\ Single\ vs\ Combined\ Dimer\ Standard\ Deviation\ Relative\ Difference','on');
%  plot(n_list,rel_dif,'.-')
%% difference as a function of N
N_max=500; % max N for which no error about computational error is raised
analytic_diff=[]; N_list=10:10:N_max;
for nn=N_list
    Analytic=[]; Analytic_comb=[];
    for k=0:nn
        Analytic(end+1)=nchoosek(nn,k)/2^nn;
        Analytic_comb(end+1)=(2*nchoosek(nn,k)*2^nn+4*nchoosek(2*nn,2*k))/(2^(2*nn+2));
    end
    analytic_diff(end+1)=sum(abs(Analytic-Analytic_comb));
end
setfig(16,'L',0,'Total\ Distribution\ Difference',...
    0,'Analytical\ Single\ Minus\ Combined\ Dimer\ Distributions','on');
plot(N_list,analytic_diff,'o-'); figname=['Toy_model/Normalized_Analytic_Combined_minus_dimer_dist'];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% fitting an exponential decay function 
% g = fittype('a+b*exp(-c*x)'); f0=fit(log10(eps_list).',dim_degen,g); f1=fit(log10(eps_list).',comb_degen,g);
% f0=fit(log10(eps_list).',dim_degen,'cubicinterp','breaks',3); f1=fit(log10(eps_list).',comb_degen,'cubicinterp');
% ft=fittype('piecewiseexp(x,a,c,d,e,k)');
% % a is constant for first half of x, c is starting point of exp decay, d is
% % exp decay coef, k is peicewise starting point, e is exp decay rate
% f0=fit(log10(eps_list).',dim_degen,ft,'StartPoint',[sample_size,sample_size,1,-2.4,1]);
% f1=fit(log10(eps_list).',comb_degen,ft);
% plot(f0,'.-'); plot(f1,'.-') 

%% PCA

K=2; dpoint_per_v=sample_size;
options = {'kmeans_dim', 3, 'symm', 'spherical'};
PCA_kmeans(samples,dpoint_per_v,K, Jlist, adjmat, savedir, options{:});

%% Trying the Autoencoder on this model (input have to be real or imag)

K=2; dpoint_per_v=sample_size; latent_dim=2;
options = {'apply_transferfunc', true, 'symm', 'circular'};
Autoenc_kmeans(samples,dpoint_per_v,latent_dim,K, Jlist, adjmat, savedir, options{:});

%% For error analysis and averaging over samples
path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
% if ~exist('groundstates','var')


% nlist=8:2:24;
n_runs=90;  % for 1000 samples, the runs take about 1 minute each
% cmap=parula(length(nlist)); lines=[];
% n_ndegen=zeros([length(nlist),n_runs,length(Jlist)]);
% for nn=1:length(nlist)
% N=nlist(nn);
% N=6;
%N=24;
% load(['plots/1dchainperiodic/Finite_size/J1J2_spinhalf_N_' num2str(N) '.mat']); 
load(['plots/1dchainperiodic/Finite_size/J1J2_spinhalf_N_24_0symm_breaking.mat'])
% end
% print(['begining N=' num2str(N) 'run' ])

eps_list=0.02;
% fig1=setfig(16,'J_2',0,'Degeneracy',0,['Max\ Eigenvalue\ Degeneracy\ for\ \epsilon=' num2str(eps_list)],'on');
% choose the value of epsilon to do error analysis for
ndegen_tot=zeros([n_runs,length(Jlist)]);
parpool(3)
parfor run=1:n_runs
    tic
    sample_size=1000; samplepos = zeros(max(size(Jlist))*sample_size,1); 
    for n=1:(max(size(Jlist)))
       GS_probs = abs(groundstates(:,n)).^2; % probabilities of each state
       for m=2:DK
            GS_probs(m) = GS_probs(m)+GS_probs(m-1);
        end
        for ii=1:sample_size
            a = rand;
            samplepos((n-1)*sample_size+ii,1) = sum(a>=GS_probs)+1;
        end
    end
    if spin==0.5
        samples=2*Zred(samplepos,:); 
    else
        samples=Zred(samplepos,:);
    end

    degen_tolerance = 10^(-2.5);

    options = {'degen_tol',degen_tolerance,'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',1,...
    'fig_title', 'degenerate\ evals=1\ vs\ J_2\','var_name','J_2'};
    [ndegenerate,~,~,~]=diffmap_list(0.5.*samples,Jlist,eps_list,options{:});

    toc

    ndegen_tot(run,:)=squeeze(ndegenerate);
    
    fprintf(['finishing run # ' num2str(run) '\n'])
    
end

delete(gcp('nocreate'));

% print(['finishing N=' num2str(N) 'run' ])

std_err=sqrt(sum((ndegen_tot-mean(ndegen_tot,1)).^2,1)/n_runs)/sqrt(n_runs);

% n_ndegen(nn,:,:)= ndegen_tot;
% sb0=ndegen_tot;
% figure; 
% lines(end+1)=plot(Jlist,mean(ndegen_tot,1)/mean(ndegen_tot(:,1),1),...
%     '.-', 'DisplayName',['N=' num2str(N)]);
% plot(Jlist,mean(ndegen_tot,1),'.-', 'DisplayName',...
%     ['Symm breaking=' num2str(J_symm_break)])
hold on
% e=errorbar(Jlist,mean(ndegen_tot,1),std_err,'.-');
% ylim([40,520])
% end % end N loop

% set(lines,{'color'},num2cell(cmap,2));
% legend('show')

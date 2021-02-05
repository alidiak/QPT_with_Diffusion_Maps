%% Init
clear ; close all; clc

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
addpath('../../Utils');

%% Creating an adjacency matrix for a simple 1d lattice

N=6; % choose system size

v=ones((N-1),1);
adjmat = diag(v,1)+diag(v,-1);

% periodic conditions
per='periodic';
%per='';
lattice='1dchain';

if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
end

G=graph(adjmat);
% p=plot(G);

time_evol=false; dt=0.1; tmax=1;

degen_per_h=true; MH=false;

%% Loading adjacency matrix from a file
% Use mathematica notebook "Lattice_adjacency_mat_generator.nb" to generate
% adjacency matrices to impart model on. 

per='periodic'; % enter '' if non-periodic and 'periodic' if periodic.
N =16;
lattice='checker';
filedat = importdata(['../../Edgemats/' lattice per '_edges_' num2str(N) '_vertices.txt'],' ' ,3);
edges = filedat.data(:,:);
% N = max(max(edges));
adjmat = sparse([edges(:,1) edges(:,2)],[edges(:,2) edges(:,1)],[edges(:,3) edges(:,3)],N,N);
G = graph(adjmat); plot(G);
time_evol=false;

%% Exact Diagonalization
Z=getZ(N,1/2); %Get the matrix elements for all {Z_m} which are diagonal

ia=(1:2^N)'; %Linear index array

%Now constructing the total Sx,Sy,Sz operator (we can easily square them without constructing two-point correlators)
Sx=spalloc(2^N,2^N,N*2^N);
% Sx_interact=spalloc(2^N,2^N,N*2^N);Z_interact= Z*adjmat;
Sy=spalloc(2^N,2^N,N*2^N);
Sz=sparse(ia,ia,sum(Z,2),2^N,2^N,2^N);

for m=1:N
    Sx=Sx+sparse(ia,ia+2*Z(:,m)*2^(N-m),1,2^N,2^N,2^N);
    % Sx_interact=Sx_interact+sparse(ia,ia+2*Z_interact(:,m)*2^(N-m),1,2^N,2^N,2^N);
    Sy=Sy+sparse(ia,ia+2*Z(:,m)*2^(N-m),1i*Z(:,m),2^N,2^N,2^N);
end

HZZ=sparse(ia,ia,sum(((Z*adjmat).*Z),2),2^N,2^N,2^N);

%% Testing my periodic kron_matrix_generator code

% sigmaz=[[1,0];[0,-1]]; op=kron(sigmaz,sigmaz);
% op_list={};
% for i=1:2
%     op_list{end+1}=sigmaz;
% end
% options={'op_list',op_list};
% J_HZZ= kron_matrix_generator(op,2,N, per,options{:});
% 
% sum(J_HZZ==2.*HZZ)-2^N % should return all 0's and no negatives
%%
Z=2*Z; % we work in the eigenbasis where spins are =+-1
%Generate the Ising ZZ interactions: sum_{m<n}lambda*J(m,n)*Z_m*Z_n (diagonal matrix)

% To break spatial symm for inversion and translation symm tests
% options={'spatial_var',true}; sz=[[1/2,0];[0,-1/2]];  
% spatial_vary_sz=kron_matrix_generator(sz,2,N,per,options{:}); 

dh = 0.1;
% hlist = 0.00:dh:2.0;
hlist=[1];
sample_size=500;

J= 0;  % 1 antiferromagnetic
hz = -0.00;% Small symmetry breaking term
Bx=0.0;

samplepos = zeros(max(size(hlist))*sample_size,1); % will have sample_size
% number of samples for each change in the variable h.

E_gap = zeros(max(size(hlist)),1);

lattice_Sz = zeros(N,max(size(hlist)));
lattice_Sx = zeros(N,max(size(hlist)));
sigmax = zeros(1,max(size(hlist)));
groundstates = zeros(2^N,max(size(hlist)));
energy=zeros(1,max(size(hlist)));

if time_evol
    t=0:dt:tmax; 
    sigmax_t=zeros(max(size(t)),max(size(hlist)));
    sigmaz_t=zeros(max(size(t)),max(size(hlist)));
    psi= zeros(max(size(t)),max(size(hlist)),2^N);
    %samplepos=zeros(max(size(t)),max(size(hlist))*sample_size);
    ED_tsamples=zeros(max(size(t)),max(size(hlist))*sample_size,N);
end 

for n=1:(max(size(hlist)))
    h = hlist(n);

    nevals = 2^N;
    H=-J*HZZ+0.5*h*Sx;
    %H=J*HZZ-Bx*Sx-h*spatial_vary_sz;
    %H=((2*J).*HZZ-0.5*Sx+h*2*Sz);
    %H=h*2*Sx;
    [evecs,evals] = eigs(H+H', nevals, 'SA');
    fprintf(['Min. Energy for h=' num2str(h) ': ' num2str(evals(1,1)) '\n'])
    
    energy(n)=evals(1,1);
    
    % energy gap calc
    % uniq=uniquetol(diag(evals));
    diagvals = diag(evals);
    E_gap(n)= diagvals(2)-diagvals(1);
    
    GS=evecs(:,1);
    groundstates(:,n)=GS;
    GS_probs = abs(GS.*1).^2; % probabilities of each state
    
    lattice_Sz(:,n)=((abs(GS.*1).^2).')*Z; 
    % blanko(n) = (GS.')*Sz*Sz*GS;
    
    %lattice_Sz(:,n)=(GS_probs.')*Z;
    % lattice_Sz(:,n)=(GS.')*Z; 
    % lattice_Sx(:,n)=GS_probs*(Sx*Z); % WRONG
       
    sigmax(n) = (GS')*Sx*GS; % average sigma X
    
    %fprintf(['SigmaZ GS Expectation value for h=' num2str(h) ': ' ...
     % num2str((GS')*Sz*GS) '\n' num2str(sum(lattice_Sz(:,n))) '\n']) %gives same result as method below
    
    % Adding the probabilities so that each state can be represented by a
    % unique probability between 0 and 1. 
    if ~time_evol
        for m=2:2^N
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
    
    % Exact Time Evolution:
    %if time_evol
    %    psi(1,n,:)= ones(1,2^N)/norm(ones(1,2^N)); % equal superposition of all states (sigmaX state)
    %    for jj=2:max(size(t))
    %        psit=psi(jj-1,n,:);
    %        psi(jj,n,:)= exp(-1i*H*squeeze(psi(jj-1,n,:))*dt);
    %    end
    %end
    
    if time_evol
        % Enter here the initial state
        %init_psi=ones(2^N,1); %equal superpos of all states
        init_psi=zeros([2^N,1]); init_psi(1) =1; %init_psi(end)=1; % cat-state, all spin up/down 
        init_psi=init_psi/norm(init_psi);
        for tt=1:max(size(t))
            psi(tt,n,:)=expm(1i*H*t(tt))*init_psi; 
            % psi(tt,n,:)=evecs*(exp(-1i*t(tt).*diag(evals)).*pinv(evecs)*(init_psi));
            ps=squeeze(psi(tt,n,:));
            sigmax_t(tt,n) =ps'*Sx*ps;
            sigmaz_t(tt,n)= ps'*Sz*ps;
            %sigmaz_t(tt,n)=sum(((abs(ps).^2).')*Z);
            
            % sampling
            psi_probs = abs(ps.*1).^2; % probabilities of each state
            for m=2:2^N
                psi_probs(m) = psi_probs(m)+psi_probs(m-1);
            end
            for ii=1:sample_size
                a = rand;
                pos = sum(a>=psi_probs)+1;
                samplepos(tt,(n-1)*sample_size+ii)=pos;
            end
        end
    end
end

if time_evol
    figure;  hold on
    for jj=1:(max(size(hlist)))
        h=hlist(jj);
        %plot(t,sum((ps*2*Sx)'.*ps',1),'o-','DisplayName',['h=' num2str(h)])
        plot(t,sigmax_t(:,jj),'-o','DisplayName',['h=' num2str(h)])
    end
    title('\sigma_x vs t'); xlabel('t'); ylabel('Sum \sigma^x');  hold off
    legend show; figure;  hold on
    for jj=1:(max(size(hlist)))
        h=hlist(jj);
        plot(t,sigmaz_t(:,jj),'o-','DisplayName',['h=' num2str(h)])
    end
    title('\sigma_z vs t'); xlabel('t'); ylabel('Sum \sigma^z');  hold off
    legend show
end

fprintf('\n\n SigmaZ GS Expectation values \n')
fprintf('h= %.2f : %.3f \n', [hlist',sum(lattice_Sz,1)']')

mkdir(['plots/' lattice per])
mkdir(['plots/' lattice per '/N_' num2str(N)])
if time_evol
    savedir0=(['plots/' lattice per '/N_' num2str(N) '/Time_Evo']);
    mkdir(savedir0);
else
    savedir0=(['plots/' lattice per '/N_' num2str(N) ]);
end

figure, plot(hlist,energy,'-o'), title('Energy vs. Transverse Field h')
savefig([savedir0 '/Energy_' per '_' lattice '_N_' num2str(N)])

figure, plot(hlist,E_gap,'-o'), title('Energy Gap vs. Transverse Field h')
savefig([savedir0 '/Energy_Gap_' per '_' lattice '_N_' num2str(N)])

figure, plot(hlist,sigmax,'-o'), title('Avg \sigma_x vs. Transverse Field h')
savefig([savedir0 '/SigmaX_' per '_' lattice '_N_' num2str(N)])

figure, plot(hlist,sum(lattice_Sz,1)','-o'), title('Avg \sigma_z vs. Transverse Field h')
savefig([savedir0 '/SigmaZ_' per '_' lattice '_N_' num2str(N)])

ED_samples = Z(samplepos(:),:); % generate all of the samples 

% two point correlator
d=distances(G);
[max_dist,ind]=max(d(:)); % finds the index of the furthest distance spots
% the column (floor(ind/N)+1) and row (ind) corresponding to the max dist.
row = rem(ind,N);
if row==0
   row =N;
end
two_point_corr = zeros([N^2,max(size(hlist))]);
for ii=0:(N-1)
    for jj=1:N
        % ii*N+jjwhy does matlab look so weird on linux
        two_point_corr(ii*N+jj,:)=lattice_Sz(ii+1,:).*lattice_Sz(jj,:);
    end
end
imaging_mat=reshape(two_point_corr,[N,N,max(size(hlist))]);
%two_point_corr = lattice_Sz(floor(ind/N)+1,:).*lattice_Sz(row,:); 
%figure, plot(hlist,two_point_corr,'-o'), title('Two Point Corr. vs. Transverse Field h')
%savefig(['plots/' lattice per 'N_' num2str(N) '/TwoPointCorr_' per '_Triangular_N_' num2str(N)])

%% RBM vs exact

% latt_name= ' periodic\ checkerboard\';
% % Triangular
% % RBM_energies=[-15.99897204, -16.55107248, -17.16353367, -17.8234492 ,       -18.5350453 , -19.2911498 , -20.10235252, -20.95091572,-21.86320315, -22.8238222 , -23.8400391 , -24.91286207,  -26.04052092, -27.21938481, -28.44083025, -29.71014847,-31.02046384, -32.35507881, -33.71940399, -35.10133036,       -36.50929159];
% RBM_energies = [-24.99856086, -25.02768697, -25.11051164, -25.25085099,  -25.43197074, -25.71570172, -26.04496388, -26.44765942,-26.91792488, -27.46615451, -28.07736657, -28.76173727,-30.30467504, -31.16793777, -32.08170623, -33.0541639 ,-34.06998884, -46.18905123, -36.27426976, -37.45565626,-38.69064916, -39.97011992, -41.28955124, -42.64408854,-44.02568664, -45.43277489, -46.86088599, -53.70287508,-49.770118  , -51.24377708];
% setfig(16,'h',0,'Energy',0,['Frustrated\ 2D\' latt_name ' TFIM'],'on'); 
% % setfig(16,'h',0,'Avg\ log(1-\lambda_2)',0,'Lambda\ measures','off');
% plot(hlist,energy,'-o','DisplayName','ED Energy'); hold on
% plot(hlist,RBM_energies,'-*','DisplayName','RBM Energy');
% legend()
% 
% setfig(14,'h',0,'Energy\ Difference',0,['Frustrated\ 2D\' latt_name ' TFIM\ \Delta\ E'],'on'); 
% plot(hlist,(RBM_energies-energy),'k-*')

%% Metropolis Hastings Sampling of ED groundstate
if MH
burnin=200; chain_size=1000; hilbspace=[1,-1];
num_starts=100; eff_size=chain_size/num_starts;

MHchains=zeros([max(size(hlist))*chain_size,N]); 
temp_chain=zeros([num_starts*eff_size,N]);
for ii=1:max(size(hlist))
    GS=groundstates(:,ii);
    psi = {GS,Z}; % have to combine the samples like this, input Z so MH 
    % alg can identify which coefficient based on a sample.

    for ll=1:num_starts % adding extra burn ins/starting points (10 to be exact)
        options={'psi',psi}; 
        burnchain=MetropolisHastings(burnin,hilbspace,options{:});

        options={'s0', burnchain(end,:),'psi',psi};
        temp_chain((ll-1)*eff_size+1:(ll*eff_size),:)=...
            MetropolisHastings(eff_size,hilbspace,options{:});
    end 
    MHchains(((ii-1)*chain_size+1):(ii*chain_size),:)=temp_chain(:,:);
end
end

%% Testing probability comparison methods on symmetry transformations

% [ks_tests,ks_stats]=symm_kstesting(ED_samples/2,-ED_samples/2,hlist,1/2,'spin flip','B_z');

options={'statistic','Cramer'};
[ks_tests,ks_stats]=symm_kstesting(samplepos,getR(flip(ED_samples,2)/2,1/2),hlist,1/2,...
    'inversion','Bz of \Sigma_i i\sigma_i^z ',options{:});

%% creating a histogram comparison of MH sampling versus direct sampling
if MH
MHsample_pos=getR(MHchains/2,1/2);
% for ii=1:max(size(MHchains))
%     [q,idx]=ismember(MHchains(ii,:),Z,'rows');
%     MHsample_pos(ii)=idx;
% end

step=5; figure; count=1;
for kk =1:step:max(size(hlist))
    subplot((round(max(size(hlist))/5)+1),2,count)
    %histogram(sum(ED_samples(((kk-1)*sample_size+1):(kk*sample_size),:),2)); 
    %title(['Sum Direct Sampling Histogram for h=' num2str(hlist(kk))])
    histogram(samplepos(((kk-1)*sample_size+1):(kk*sample_size)),50); 
    title(['Direct Sampling Evec Histogram for h=' num2str(hlist(kk))])
    subplot((round(max(size(hlist))/5)+1),2,count+1)
    histogram(MHsample_pos(((kk-1)*chain_size+1):(kk*chain_size)),50); 
    title(['MH Sampling Evec Histogram for h=' num2str(hlist(kk))])
    %histogram(sum(MHchains(((kk-1)*chain_size+1):(kk*chain_size),:),2)); 
    %title(['Sum MH Sampling Histogram for h=' num2str(hlist(kk))])
    count=count+2;
end
savefig([savedir0 '/histogram_evecs_MH_vs_DS_' ...
    per '_' lattice '_N_' num2str(N) '_' num2str(num_starts) '_MH_restarts'])
end
%% Diffusion Map

% to test MH on ED generated samples
if MH
    samples=MHchains;
    savestr='MH';
else
    savestr='';
    samples=ED_samples;
end
if ~time_evol
    t=1; dt=0; ndegenerate=zeros([max(size(hlist)),max(size(t))]);
end

labels=zeros([1,size(samples,1)]); 
for ii=1:max(size(hlist))
    labels(((ii-1)*sample_size+1):((ii)*sample_size))=hlist(ii);
end

mkdir([savedir0 '/diffusionmap/'])
save_dir1=[savedir0 '/diffusionmap/range_' num2str(max(size(hlist)))];
mkdir(save_dir1)
Y_avg=zeros([length(hlist),1]); Y_uniq=zeros([length(hlist),1]); 
Ys=zeros([sample_size,sample_size,length(hlist)]); 
for eps_epsstar=1:0.5:6
if degen_per_h

    nevals=30;
    ndegenerate=[]; degen_tol=1e-4;

    for kk =1:max(size(hlist)) % applying diffusion maps to each Jx
        
        S=samples(((kk-1)*sample_size+1):((kk)*sample_size),:);
        % randomize samples and labels
        k=randperm(max(size(S)));
        S=S(k,:);
        % labels=labels(k);

        epsilon=eps_epsstar*(2*pi/sample_size);
        n_comp=5; alpha=0; % so that diffusion map is A_ll' 
        [vect2, vals2,Ys(:,:,kk)]= diffusionmaps(S',epsilon,alpha,n_comp);

        Y_avg(kk)=mean(mean(Y)); Y_uniq(kk)=length(unique(Y));
                
        ndegenerate(end+1)=sum(abs(ones(size(vals2))-vals2)<degen_tol);
      
    end

    save_dir=[save_dir1 '/degenerate_evals'];
    mkdir(save_dir)

    figure, plot(hlist, ndegenerate,'-*'), xlabel('Jx'),ylabel('# of evals=1')
    title(['number of degenerate max evals vs Z Ising coupling Jx for \epsilon/\epsilon* =' num2str(eps_epsstar)])
    figname=[save_dir '/numdegenevals_' per '_N_' num2str(N) '_epsepsstar_' num2str(eps_epsstar) savestr];
    saveas(gcf,[figname '.png'],'png') %savefig([figname '.fig']); 

else 
    
    % randomize samples and labels
    k=randperm(max(size(labels)));
    samples=samples(k,:);
    labels=labels(k);

    nevals=30;
    epsilon=eps_epsstar*(2*pi/sample_size);
    n_comp=5; alpha=0; % so that diffusion map is A_ll' 
    [vect2, vals2]= diffusionmaps(samples',epsilon,alpha,n_comp);

    save_dir=[save_dir1 '/' num2str(hlist(1)) '_' num2str(hlist(2))];
    mkdir(save_dir)

    figure
    plot(1:nevals,vals2(1:nevals),'o')
    title(['Diffusion Map evals of Ising Lattice with N = ' num2str(N) ' and eps/eps*=' ...
        num2str(eps_epsstar)])
    savefig([save_dir '/numevals_per_' per '_N_' num2str(N) '_epsepsstar_' num2str(eps_epsstar) '.fig'])
    figure; 
    gscatter(vect2(:,2),vect2(:,3),labels,autumn(max(size(hlist))),'*',3)%,'rkgb','*',5) %the 1st/0th component is always constant ]
    title(['Diffusion Map evecs 2 and 3 eps/eps*='  num2str(eps_epsstar)])
    xlabel('vector2'); ylabel('vector3')
    savefig([save_dir '/per_' per '_N_' num2str(N) '_epsepsstar_' num2str(eps_epsstar) '.fig'])
end 
end

%% PCA
%samples=samples-mean(samples,1);  % subtracting the mean first
[U,S] = svd((samples.')*samples);

figure, ax1 = subplot(3,1,1); ax2 = subplot(3,1,2); ax3=subplot(3,1,3);
plot(ax1,1:N,diag(S)/norm(diag(S)),'-o')
title(ax1,['PCA Evals for ' per lattice ' Lattice N= ' num2str(N)])
ylabel(ax1,'Ordered Evals')

plot(ax2,1:N,U(:,1),'-*'), 
title(ax2,['Largest Evec for ' per 'Triangular Lattice N= ' num2str(N)])
xlabel(ax2,'Lattice Site'), ylabel(ax2,'Lattice Site Evec Val.')

plot(ax3,1:N,U(:,2),'-*'), 
title(ax3,['2nd Largest Evec for ' per 'Triangular Lattice N= ' num2str(N)])
xlabel(ax3,'Lattice Site'), ylabel(ax3,'Lattice Site Evec Val.')
savefig([savedir0 '/PCA_Results_' per '_' lattice '_N_' num2str(N)])

figure, plot(hlist,(U(:,1)'*lattice_Sz)','-o'), title('Sum_i (Evec1_i*\sigma_i^Z) vs. Transverse Field h')
savefig([savedir0 '/PCA_SigmaZ_Projection' per '_' lattice '_N_' num2str(N)])

figure, plot(hlist,abs(mean(squeeze(reshape(samples*U(:,1),...
    [sample_size,size(hlist)])))),'-o','DisplayName',['Evec' num2str(1)])
hold on
plot(hlist,abs(mean(squeeze(reshape(samples*U(:,2),...
    [sample_size,size(hlist)])))),'-o','DisplayName',['Evec' num2str(2)])
hold off, legend('show')
title('Data projection of evecs as a function of h')
savefig([savedir0 '/PCA_data_Projection' per '_' lattice '_N_' num2str(N)])

% Plotting it on the non-periodic latttice regardless of if it is periodic
% or not
if ~strcmp(lattice,'1dchain')
nonper = importdata(['../Edgemats/' lattice '_edges_' num2str(N) '_vertices.txt'],' ' ,3);
edges = nonper.data(:,:);
adjmat = sparse([edges(:,1) edges(:,2)],[edges(:,2) edges(:,1)],[edges(:,3) edges(:,3)],N,N);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
G = graph(adjmat);
p=plot(G);
G.Nodes.NodeColors =U(:,1);
p.NodeCData = G.Nodes.NodeColors;
colorbar
title('Colormap for largest evec')
savefig([savedir0 '/Graph_' per '_' lattice '_N_' num2str(N) '_1stEvec'])

figure
p2=plot(G);
G.Nodes.NodeColors = U(:,2);
p2.NodeCData = G.Nodes.NodeColors;
colorbar
title('Colormap for 2nd largest evec')
savefig([savedir0 '/Graph_' per '_' lattice '_N_' num2str(N) '_2ndEvec'])

dpoints_per_v=500; % number of points to plot per change in variable
skip=sample_size/dpoints_per_v; % the number of points to be skipped in samples
ncolors=size(samples,1)/sample_size; % number of colors in dataset - 1 per variable change
C=repmat(1:ncolors,sample_size/skip,1); % for mapping each array to a color
figure
y1=samples(1:skip:end,:)*U(:,1);y2=samples(1:skip:end,:)*U(:,2);
y3=samples(1:skip:end,:)*U(:,3); % projections onto PCA subspace
s=repmat(5*ones(ncolors,1),sample_size/skip,1); % dot size (all same here)
scatter3(y1,y2,y3,s(:),C(:)) 
xlabel('Evec1 Projection'); ylabel('Evec2 Projection'), zlabel('Evec3 Projection')
title('Data projection on PCA for evecs 123')
colormap(jet(ncolors))
c=colorbar('Ticks',1:2:ncolors,'TickLabels',hlist(1:2:end));
c.Label.String='Transverse field Jx';
savefig([savedir0 '/PCA_Projection_evecs'])


%% checking MPS psi functions
% if strcmp(lattice,'1dchain')
%     [all_psi,v_list]=get_MPS_psi('TFIM',N); % generates MPS psi from .log data
% 
%     %evals=unique(Z);
%     evals=[1;-1];
%     mps_groundstates= zeros([2^N,max(size(v_list))]);
%     for ii=1:max(size(v_list))
%         for jj=1:2^N
%             s = Z(jj,:);
%             mps_groundstates(jj,ii)=MPS_psi(all_psi{ii},s,evals); 
%         end
%     end
% end

%% Testing probability comparison methods
% compareto='10';
% 
% if max(hlist)==10
%     samplepos_B10=samplepos; save("N"+num2str(N)+"_B10_samples.mat",'samplepos_B10')
% elseif min(hlist)==0
%     samplepos_B0=samplepos; save("N"+num2str(N)+"_B0_samples.mat",'samplepos_B0')
% elseif min(hlist)==0.1
%     samplepos_B01=samplepos; save("N"+num2str(N)+"_B01_samples.mat",'samplepos_B01')
% end
% if strcmp(compareto,'10')
%     load("N"+num2str(N)+"_B10_samples.mat"); cmprto=samplepos_B10;
%     % this is the dist I will compare to the others
% elseif strcmp(compareto,'0')
%     load("N"+num2str(N)+"_B0_samples.mat"); cmprto=samplepos_B0;
% elseif strcmp(compareto,'0.1')
%     load("N"+num2str(N)+"_B01_samples.mat"); cmprto=samplepos_B01;
% end

% apply symmetry here
% ED_symm_samples=
% 
% KLdivergence=zeros([max(size(hlist)),1]); total_var=zeros([max(size(hlist)),1]); 
% ks_tests=zeros([max(size(hlist)),1]); ks_hypothesis=zeros([max(size(hlist)),1]); 
% ks_stat=zeros([max(size(hlist)),1]); 
% t_tests=zeros([max(size(hlist)),1]); t_hypothesis=zeros([max(size(hlist)),1]); 
% t_stat=zeros([max(size(hlist)),1]); t_ci=zeros([max(size(hlist)),1]); 
% kuip_stat=zeros([max(size(hlist)),1]);  cramer_stat=zeros([max(size(hlist)),1]); options={'cramer',true};
% for kk=1:max(size(hlist))
%     samplepos_h=samplepos(((kk-1)*sample_size+1):(kk*sample_size));
%     KLdivergence(kk)=KLDiv(samplepos_h.',cmprto.');
%     kuip_stat(kk)=kuipertest2(samplepos_h,cmprto,500,false);
%     cramer_stat(kk)=kuipertest2(samplepos_h,cmprto,500,false,options{:});
%     [ks_hypothesis(kk),ks_tests(kk),ks_stat(kk)]=kstest2(samplepos_h,cmprto);
%     [t_hypothesis(kk),t_tests(kk),ci,stat]=ttest2(samplepos_h,cmprto);
%     t_stat(kk)=stat.tstat; t_ci(kk)=mean(ci);
%     total_var(kk)=total_variation_dist(samplepos_h.',cmprto.');
% end

% subplot(4,1,1);plot(hlist,KLdivergence,'-*');  title(['Kullback-Leiber Divergence comparing against ' compareto])
% figure; subplot(4,1,1);plot(hlist,t_tests,'-*');  title(['T-Test p-value comparing against ' compareto])
% subplot(4,1,2); plot(hlist,log(abs(t_stat)),'-*'); title(['T-Test statistic comparing against ' compareto])
% figure; subplot(4,1,1);plot(hlist,kuip_stat,'-*');  title(['Kuiper Statistic comparing against B_x=' compareto])
% subplot(4,1,2);plot(hlist,cramer_stat,'-*');  title(['Cramer Statistic comparing against B_x=' compareto])
% subplot(4,1,3); plot(hlist,log(ks_tests),'-*'); title(['kolmogorov-smirnov test p value comparing against B_x=' compareto])
% subplot(4,1,4); plot(hlist,ks_stat,'-*'); title(['kolmogorov-smirnov statistic comparing against B_x=' compareto])
% % subplot(4,1,4); plot(hlist,total_var,'-*'); title(['Total variation distance of probability measures comparing against ' compareto])
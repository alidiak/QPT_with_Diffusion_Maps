%% Init
clear ; close all; clc

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
addpath('../../Utils');

%% Creating an adjacency matrix for a simple 1d lattice

N=12; % choose system size
spin=1/2; % choose spin
num_seeds=20; % choose the number of random seeds/runs to average over
ensemble='long_time';  % choose an ensemble
state='neel'; % pick neel for neel order state, half for half up, second half down initial state

v=ones((N-1),1); v2=ones((N-2),1);
adjmat = diag(v,1)+diag(v,-1);

% periodic conditions
per='';
lattice='1dchain';

if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
end

G=graph(adjmat);% p=plot(G);

degen_per_h=true; MH=true;

%% Exact Diagonalization Matrices 
% warning('off','all')

Z=getZ(N,spin); Dim=2*spin+1; D=Dim^N; %Get the matrix elements for all {Z_m} which are diagonal


ia=(1:D).'; % Linear index array
 if exist('KK','var')
    chi = (sum(Z,2)==0); Zred=Z(chi,:);
 else
    chi=ones([D,1]); Zred=Z(ia,:);
end
DK=sum(chi); % the reduced matrix size
ja= (1:DK)';
xi=zeros(D,1); xi(chi)=ja; %xi maps the indices from original to projected space


%Now constructing the total Sx,Sy,Sz operator (we can easily square them without constructing two-point correlators)
Sx=spalloc(D,D,N*D);
Sy=spalloc(D,D,N*D);

for m=1:N
    Sx=Sx+sparse(ia,ia+Dim*Z(:,m)*Dim^(N-m),1,D,D,D);
    Sy=Sy+sparse(ia,ia+Dim*Z(:,m)*Dim^(N-m),1i*Z(:,m),D,D,D);
end

% Sx=Sx(ki,ki); Sy=Sy(ki,ki); % goes to zero when KK=0, and the Hamiltonian looses its diagonal elements

% HZZ=sparse(ja,ja,sum(((Zred*adjmat).*Zred),2),DK,DK,DK);
% Zred=2*Zred; % we work in the eigenbasis where spins are =+-1
%Generate the Ising ZZ interactions: sum_{m<n}lambda*J(m,n)*Z_m*Z_n (diagonal matrix)

%Now constructing the total Sx,Sy,Sz operator (we can easily square them without constructing two-point correlators)
Sz=sparse(ja,ja,sum(Zred,2),DK,DK,DK); % should be all 0 sparse if K=0

%% Finding Energies and Creating Samples

J= 1.0;  % nearest neighbor Heisenberg interaction (positive for anti-ferro)
hx=0.6; % using same values as in paper
J2=0.3; % next nearest neighbor inter

sigmaz=[[1,0];[0,-1]]; sigmax=[[0,1];[1,0]];

% uncomment below for a flipped version corresponding to measurements in the sx direction
J2_HXX= J2.*kron_matrix_generator(kron(sigmax,kron(eye(Dim),sigmax)),Dim,N, per); % checked for simple N=3, seems to work

dh = 0.25;
% hlist =0.0:dh:8.0; % here h is the randomization in J1 (the nearest neighbor interaction term)
% hlist=logspace(-1,1.1); % matching what is done in the paper
hlist= 0.1:dh:10.1;
% hlist=[0];

sample_size=500;
deps=2*3.1250e-04;
% eps_list=0.00025:deps:0.020; % normal method seems to work better at lower eps
% eps_list=logspace(-2.4,-1.8,40);

eps_list=0.005; 

savedir=['plots/' ensemble '_ensemble/' lattice per '/' state '/N_' num2str(N) ];
mkdir(savedir)

E_gap_avg = zeros(max(size(hlist)),1);
energy_avg=zeros(max(size(hlist)),1);
ndegenerate_avg=zeros(length(eps_list),max(size(hlist)));
secondeval_avg=zeros([length(eps_list),length(hlist)]);
eval_avg_avg=zeros([length(eps_list),length(hlist)]);
Kpred_avg=zeros(max(size(hlist)),1);
PCA_pred_avg=zeros(max(size(hlist)),1);
nunique_avg=zeros([length(hlist),1]);
dist_to_initial_state = zeros([sample_size*length(hlist),1]);

% for saving and recreating purposes (so I don't have to rerun the full ED again)
wavefunctions=zeros([num_seeds, length(hlist), DK]);
samplepos_tot=zeros([num_seeds, length(hlist)*sample_size]);

if strcmp(state,'neel') % creating the neel state, alternating up and down or vise versa
    mix=0.5*ones([1,N]); mix(2:2:end)=-mix(2:2:end); 
else  % creating an initial state with half down and half up
    mix=zeros([1,N]);
    mix(1:(N/2))=0.5; mix((N/2+1):N)=-0.5;
end
if sum(mix)==0
    init_psi=(N==sum(Zred==mix,2)); 
else
    init_psi=(N==sum(Z==mix,2)); 
end
init_state=find(init_psi); % just to keep track of which state was found

parpool(5)
parfor seed=1:num_seeds
    
warning('off','all')

fprintf('Beginning  run %.2f out of %.3f \n', [seed,num_seeds]')

samplepos = zeros(max(size(hlist))*sample_size,1); % will have sample_size

E_gap = zeros(max(size(hlist)),1);

lattice_Sz = zeros(N,max(size(hlist)));
lattice_Sx = zeros(N,max(size(hlist)));
states = zeros(max(size(hlist)),DK);
energy=zeros(1,max(size(hlist)));

for n=1:(max(size(hlist)))
    
    h = hlist(n); % here h is the range of the uniformly distributed random variables on Sz
    J_lattice= J - h + (2*h)*rand(N,1); % array of uniformly drawn random #s applied to lattice sites
    
    % Was for comparing to other TFIM code 
%      hx = hlist(n); 
%     J_lattice=ones(size(J_lattice));
    options={'spatial_var',J_lattice};
    % Normal version
%    J_HZZ= kron_matrix_generator(kron(sigmaz,sigmaz),Dim,N, per, options{:}); 
    % Flipped so we can measure in sx direction
    J_HXX= kron_matrix_generator(kron(sigmax,sigmax),Dim,N, per, options{:}); 

    H_int=J_HXX+J2_HXX+hx*Sz;
    H=(H_int+H_int')./2;
    
%     nevals = 2;
    [evecs,evals] = eigs(H, DK, 'SA');
    fprintf(['Min. Energy for h=' num2str(h) ': ' num2str(evals(1,1)) '\n'])
    
    energy(n)=evals(1,1);
    
    % energy gap calc
    diagvals = diag(evals);
    E_gap(n)= diagvals(2)-diagvals(1);
    
    % For the diagonal ensemble, define an initial state - at t->inf this state will still be prevalent 
%     init_psi=zeros([1,length(chi)]); init_psi(2) =1; %init_psi(end)=1; % cat-state, all spin up/down 

%     init_psi=init_psi/norm(init_psi);
    
    % The alternative long time evolution of initial psi
    if strcmp(ensemble,'long_time')
        t=10000;
%         psi_E  = expm(1i*H*t)*init_psi; % note that there are many far more efficient ways to do this      
        psi_E  = (evecs*diag(exp(1i*diagvals*t))/evecs)*init_psi;  % using the already diagonalized H will be a bit faster and memory efficient
    elseif strcmp(ensemble,'diag')
    % Project the initial state onto the energy basis i.e. get all the coefficients c_0, c_1,... c_(chi). 
        psi_E=evecs*init_psi;
    end
    
    states(n,:)=psi_E;
    
%     wavefunction(seed,n,:)=psi_E;
    state_probs=abs(psi_E).^2; % same as diag(psi_E.' * psi_E)

    for m=2:DK
        state_probs(m) = state_probs(m)+state_probs(m-1);
    end
    for ii=1:sample_size
    
    % Sampling %
    a = rand;
    samplepos((n-1)*sample_size+ii,1) = sum(a>=state_probs)+1;
    % finds the pos where a is greater than the lowest prob. This is the 
    % sample we choose by picking a random number a.
    end
end

wavefunctions(seed,:,:)=states;
samplepos_tot(seed,:)=samplepos;

samples=Zred(samplepos,:); 

%% Original Distance measure (Euclidian)
options = {'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',1,...
'fig_title', 'number of degenerate max evals vs h'};
[ndegenerate,secondeval,eval_avg]=diffmap_list(samples,hlist,eps_list,options{:});

%% number unique
nunique=zeros([length(hlist),1]);
for n=1:length(hlist)
    nunique(n)=length(unique(samplepos(((n-1)*sample_size+1):(n*sample_size))));
end
setfig(16,'J_2',0,'#\ unique\ samples',0,'Exponential\ Scaling\ of\ Degeneracy','off');
plot(hlist,nunique,'k-*'); figname=[savedir '/number_unique_samples_' per '_' lattice '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
%% PCA

K=2; dpoint_per_v=sample_size;
options = {'kmeans_dim', 3, 'symm', 'spherical','no_plot',true};
[Kpred,PCA_pred]=PCA_kmeans(samples,dpoint_per_v,K, hlist, adjmat, savedir, options{:});

%% averaging all the important vals

E_gap_avg = E_gap_avg+E_gap;
energy_avg=energy_avg+energy;
ndegenerate_avg=ndegenerate_avg+ndegenerate;
secondeval_avg=secondeval_avg+secondeval;
eval_avg_avg=eval_avg_avg+eval_avg;
Kpred_avg=Kpred_avg+Kpred;
PCA_pred_avg=PCA_pred_avg+PCA_pred;
nunique_avg=nunique_avg+nunique;
dist_to_initial_state=dist_to_initial_state+sum((mix-samples).^2,2);
% if spacings
% level_spacing_avg=level_spacing_avg+level_spacings;
% end
fprintf('Run # %.2f out of %.3f done \n', [seed,num_seeds]')

end
 
delete(gcp('nocreate'));

% saving the data from the run for potential finite size scaling 
save(['plots/' ensemble '_ensemble/N_' num2str(N) state '.mat']);

%% Loading previous run
path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));

N=12; state='neel'; 
load(['plots/long_time_ensemble/N_' num2str(N) state '_J1rand.mat']);

%% Rerunning diffusion maps on the sample data with different parameters
degen_tolerance=10^(-2.5); alpha=1;
eps_list=logspace(-1.75,-1.48,20);
% eps_list=10^(-2);

savedir=['plots/' ensemble '_ensemble/' lattice per '/' state '/tol_' num2str(degen_tolerance)...
    '/N_' num2str(N) 'alpha_' num2str(alpha) ]; mkdir(savedir)
options = {'degen_tol',degen_tolerance,'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',alpha,...
'fig_title', 'number of degenerate max evals vs h'};
ndegenerate_avg=zeros(length(eps_list),max(size(hlist))); secondeval_avg=zeros([length(eps_list),length(hlist)]);
distance_avgs=zeros([N/2+1,length(hlist)]);
% evals_tot_avg=zeros([length(eps_list),length(hlist), sample_size]);
parpool(4)
parfor seed=1:num_seeds
    warning('off','MATLAB:MKDIR:DirectoryExists')
    samples=Zred(samplepos_tot(seed,:),:);
    
    [ndegenerate,secondeval,~,Ys]=diffmap_list(samples,hlist,eps_list,options{:});
    
    ndegenerate_avg=ndegenerate_avg+ndegenerate; 
    secondeval_avg=secondeval_avg+secondeval;
    
    % for histogramming the distance 
    distances=zeros([N/2+1,length(hlist)]); i=1;
    for k=0:2:N
      distances(i,:)=squeeze(sum(sum(abs(N*Ys-k)<1e-4))); % BE CAREFUL HERE ABOUT NORMALIZATION USED IN DIFFMAPS.m
      i=i+1;
    end
    distance_avgs=distance_avgs+distances;
%     evals_tot_avg=evals_tot_avg+evals;
end
delete(gcp('nocreate'));

%% Degeneracies versus epsilon 

setfig(16,'Log_{10}(\epsilon)',0,'Log_{10}(1-\lambda_1)',0,'Exponential\ scaling\ of\ 1-\lambda_1','on');
lines=[]; skip=2; 
for kk=1:skip:length(hlist)
lines(end+1)=plot(log10(eps_list), log10(1-secondeval_avg(:,kk)/num_seeds),'.-','DisplayName',['\delta J=' num2str(hlist(kk))]);
hold on
end
c_map=parula(length(lines)); % c_map(find(hlist==0.5)/skip,:)=[1,0,0]; % makes J_2=0.5 red
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/1minuslambda2_exponentialscaling_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

lines=[]; 
setfig(16,'Log_{10}(\epsilon)',0,'Log_{10}(Degeneracy)',0,'Exponential\ Scaling\ of\ Degeneracy','on');
for kk=1:skip:length(hlist)
        lines(end+1)=plot(log10(eps_list), log10(ndegenerate_avg(:,kk)/num_seeds),'.-','MarkerSize',10,'DisplayName',['\delta J=' num2str(hlist(kk))]);
        hold on;
end
c_map=parula(length(lines)); %c_map(find(abs(v_list-0.5)<1e-4)/skip,:)=[1,0,0]; % makes J_2=0.5 red
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/ndegen_exponentialscaling_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% plotting as 3D

% [X,Y] = meshgrid(hlist,eps_list);
% ax1=figure; surf(X,Y,log10(ndegenerate_avg/num_seeds));
% xlabel('h'); ylabel('\epsilon'); zlabel('Log_{10}(Degeneracy)')

%% Histogram of the distance of the samples to the initial state (and just distances present in general)

figure=setfig(16,'h',0,'States\ Equal\ Distance\ d',0,['Inter-Sample\ Distances;\ State\ ' state],'on'); 
for k=1:(N/2+1)
    plot(hlist,distance_avgs(k,:)/num_seeds,'-o','DisplayName',['d=' num2str((k-1)*2) ])
    hold on
end
legend; figname=[savedir '/distances_vs_h_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

figure=setfig(16,'Distance',0,'States\ Equal\ Distance\ d',0,['Inter-Sample\ Distances;\ State\ ' state],'on'); 
for kk=1:10:length(hlist)
    plot(0:2:N,distance_avgs(:,kk)/num_seeds,'-o','DisplayName',['h=' num2str(hlist(kk)) ])
    hold on
end
legend(); figname=[savedir '/Distance_distribution_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

figure=setfig(16,'Average\ Distance',0,'Counts',0,['Sample\ Distances\ to\ Initial\ State:\ ' state],'on'); 
for kk=1:10:length(hlist)
    histogram(dist_to_initial_state((kk-1)*sample_size+1:kk*sample_size)./20,'FaceAlpha',0.5,...
      'DisplayName',['h=' num2str(hlist(kk))]);
    hold on
end
legend(); figname=[savedir '/Hist_dist_to_init_state_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Plotting averages


figure; plot(hlist,nunique_avg.'/num_seeds,'k.-'); title('Number of Unique Samples Avg')
saveas(gcf,[savedir '/avg_unique_samples.png'],'png')

for ee =1:length(eps_list)
    epsilon=eps_list(ee);
    eps_epsstar=(sample_size/(2*pi))*epsilon;
    fig1=setfig(16,'h',0,'Degeneracy',0,['Avg\ Degeneracy\ for\ \epsilon=' num2str(epsilon)],'off'); 
    plot(hlist,ndegenerate_avg(ee,:)/num_seeds,'.-')
    figname=[savedir '/diffusionmap/degenerate_evals/Ndegen_Evals_epsepsstar_' num2str(eps_epsstar)];
    savefig(fig1,[figname '.fig']); saveas(fig1,[figname '.png'],'png')
    
end

figure; plot(hlist,energy_avg.'/num_seeds,'k.-'); title('Energy avg'); saveas(gcf,[savedir '/energy.png'],'png')
figure; plot(hlist,E_gap_avg.'/num_seeds,'k.-'); title('Energy Gap'); saveas(gcf,[savedir '/E_gap.png'],'png')
figure; plot(hlist,Kpred_avg.'/num_seeds,'k.-'); title('K-means Prediction avg'); saveas(gcf,[savedir '/PCA_Kmeans.png'],'png')
figure; plot(hlist,PCA_pred_avg.'/num_seeds,'k.-'); title('PCA prediction avg'); saveas(gcf,[savedir '/PCApred.png'],'png')

% Plotting overlayed numdegen lines from different epsilon regimes
to_25=0.0025/deps; 
figure; plot(hlist, [ndegenerate_avg(1,:)./mean(ndegenerate_avg(1,:)); ndegenerate_avg(1+to_25,:)./mean(ndegenerate_avg(1+to_25,:))...
    ; ndegenerate_avg(1+2*to_25,:)./mean(ndegenerate_avg(1+2*to_25,:))],'.-')
legend({['\epsilon=' num2str(eps_list(1))],['\epsilon=' num2str(eps_list(1+to_25))],['\epsilon=' num2str(eps_list(1+2*to_25))]})
title('number degenerate overlayed for different \epsilon regimes'); 
saveas(gcf,[savedir '/NdegenComb.png'],'png');savefig([savedir '/NdegenComb.fig'])

%% Histogram of the distance of the samples to the initial state

figure=setfig(16,'Average\ Distance',0,'Counts',0,['Sample\ Distances\ to\ Initial\ State:\ ' state],'on'); 
for kk=1:20:length(hlist)
    histogram(dist_to_initial_state((kk-1)*sample_size+1:kk*sample_size)./20,'FaceAlpha',0.5,...
      'DisplayName',['h=' num2str(hlist(kk))]);
    hold on
end
legend(); figname=[savedir '/Hist_dist_to_init_state_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Trying the Autoencoder on this model (input have to be real or imag)

% K=2; dpoint_per_v=sample_size; latent_dim=2;
% options = {'apply_transferfunc', true, 'symm', 'circular'};
% Autoenc_kmeans(samples,dpoint_per_v,latent_dim,K, hlist, adjmat, savedir, options{:});

%% For error analysis 
path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
if ~exist('wavefunctions','var')
    N=16; state='neel'; % _local
    load(['plots/long_time_ensemble/N_' num2str(N) state '.mat']);
end

n_runs=50;  
eps_list=0.031; 
degen_tolerance = 10^(-2.5);

options = {'degen_tol',degen_tolerance,'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',1,...
'fig_title', 'degenerate\ evals=1\ vs\ J_2\','var_name','J_2'};
ndegen_tot=zeros([n_runs,num_seeds,length(hlist)]); ndegen_seeds=zeros([num_seeds,length(hlist)]);
% parpool(3)
parpool(4)
for run=1:n_runs
    warning('off','MATLAB:MKDIR:DirectoryExists')
    tic
    sample_size=500;
    samplepos_tot = zeros(num_seeds, max(size(hlist))*sample_size); % will have sample_size
    states=zeros(length(hlist),DK);
    for seed=1:num_seeds
        for hh=1:length(hlist)
            psi_E=wavefunctions(seed,hh,:);
            state_probs=abs(psi_E).^2;
            for m=2:DK
                state_probs(m) = state_probs(m)+state_probs(m-1);
            end
            for ii=1:sample_size
            a = rand;
            samplepos_tot(seed,(hh-1)*sample_size+ii,1) = sum(a>=state_probs)+1;
            end
        end
    end
    parfor seed=1:num_seeds
        warning('off','MATLAB:MKDIR:DirectoryExists')
    
        samples=Zred(samplepos_tot(seed,:),:);
    
        [ndegenerate,~,~,~,~,~]=diffmap_list(samples,hlist,eps_list,options{:});
%     evals,Ys,Ls,diffmap_vecs
    
        ndegen_seeds(seed,:)=squeeze(ndegenerate);
        
    end
    fprintf(['run # ' num2str(run) ' finished: '])
    toc

    ndegen_tot(run,:,:)=ndegen_seeds;

end

delete(gcp('nocreate'));

ndegen_tot= reshape(ndegen_tot,[n_runs*num_seeds,length(hlist)]);

% std_err=sqrt(sum((ndegen_tot-mean(ndegen_tot,1)).^2,1)/n_runs)/sqrt(n_runs);
std_err= std(ndegen_tot,1)/sqrt(size(ndegen_tot,1));
figure; plot(hlist,std_err,'.-')

fig1=setfig(16,'J_2',0,'Degeneracy',0,['Max\ Eigenvalue\ Degeneracy\ for\ \epsilon=' num2str(eps_list)],'on');
e=errorbar(hlist,mean(ndegen_tot,1),std_err,'.-');


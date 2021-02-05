%% Init
clear ; close all; clc

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
addpath('../../Utils');

%% Creating an adjacency matrix for a simple 1d lattice

N=12; % choose system size
clocknum=3; % The n in the clock state references, should = the Z symmetry
% subspacesz=0; % constraints to reduce the subspace (and hopefully improve accuracy)
% subspacesig=0;

gridsize=41;

dtheta=pi/(clocknum*gridsize); % change in the angles
theta_list=0:dtheta:pi/clocknum; % due to symmetry, pi/3 is the max unique value for the angles
% theta_list=pi/4;
phi_list=theta_list;

sample_size=10; % MH_chain size and largest dimension of diffusion maps
MH=false; % set to true if MH sampling is desired
num_starts=100; 
% eps_epsstar_list =2.0:0.25:6.5; eps_list=eps_epsstar_list*(2*pi/sample_size); % set the variation of epsilon/epsilon*
% eps_list=logspace(-1.24,-1.14,40); eps_epsstar_list=eps_list/(2*pi/chain_size);
eps_epsstar_list=[]; % set to empty if just ks and autoencoder etc are needed

df = 1/(gridsize-1); 
flist = 0.1:df:1.0; % vary over J-f

ndegenerate=zeros([max(size(eps_epsstar_list)),max(size(theta_list)),max(size(flist))]);
clusterstdv_interclusterdist_ratio=zeros([max(size(eps_epsstar_list)),...
                            max(size(theta_list)),max(size(flist))]);

PCAkpred2=zeros([max(size(theta_list)),max(size(flist))]);
PCAkpred3=zeros([max(size(theta_list)),max(size(flist))]);                   
autoenckpred2=zeros([max(size(theta_list)),max(size(flist))]);
autoenckpred3=zeros([max(size(theta_list)),max(size(flist))]);
PCA_pred =zeros([max(size(theta_list)),max(size(flist))]);
autoenc_pred=zeros([max(size(theta_list)),max(size(flist))]);

ks_stats=zeros([max(size(theta_list)),max(size(flist))]); 
ks_p=zeros([max(size(theta_list)),max(size(flist))]);
                        
per='';

v=ones((N-1),1);
adjmat = diag(v,1)+diag(v,-1);

% periodic conditions
if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
end 

lattice='1dchain';

%%%%%%%%%%% DEFINING THE ELEMENTS OF THE CHIRAL CLOCK MODEL: %%%%%%%%%%%%%%
Jsig=0.00; % small symmetry breaking term
omega=exp(2i*pi/clocknum); % the omega term which depends on n
omega_list=[]; spins=[];
for ii=0:clocknum-1
    omega_list(end+1)=omega^(ii);
    spins(end+1)=(clocknum-1)/2-ii;
end

sigma = diag(omega_list); % sigma is the val of the angle
tau = diag(ones([clocknum-1,1]),-1); tau(1,end)=1;

D=clocknum;
if mod(clocknum,2)==1
    S=(D-1)/2; Z=getZ(N,S); 
    if exist('subspacesz','var')
        subspace=(sum(Z,2)==subspacesz);
        ia=(1:D^N)'; chi=ia(subspace);
        Zred=Z(chi,:);
%         sigmas=Zred;
%         for jj=1:clocknum
%             sigmas(sigmas==spins(jj))=omega^(jj-1);
%         end
    else
        chi=(1:D^N)';
    end
end

% Changing permuatations to the eigenvalues of the angles/sigma
% if mod(clocknum,2)==0
sigmas=permn(1:clocknum, N);
for jj=1:clocknum
    sigmas(sigmas==jj)=omega_list(jj);
end
% end

% if exist('subspacesig','var')
%     s0_condition=(sum(sigmas,2)==subspacesig);
% end

%G=graph(adjmat);
%p=plot(G);
samplepos_total=zeros([max(size(theta_list)),max(size(flist))*sample_size]);
EE=zeros([max(size(theta_list)),max(size(flist))]);

for tt=1:max(size(theta_list))
tic
theta = theta_list(tt); % interaction term phase - quantifies chirality when phi=0
phi = theta;%phi_list(tt); % 'rotating' term phase, 
% at both theta and phi= +-pi/3 it is antiferromagnetic!

% Normal version
% sum_tauj = kron_matrix_generator(tau*exp(-1i*phi)+tau'*exp(1i*phi),D,N, per);
% sum_sigmajj1 = kron_matrix_generator(kron(sigma',sigma)*exp(-1i*theta)+...
% kron(sigma,sigma')*exp(1i*theta),D,N,per);
% 
% sum_sigmaj= kron_matrix_generator(sigma+sigma',D,N,per); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tau/swapped version
sum_tauj = kron_matrix_generator(sigma*exp(-1i*phi)+sigma'*exp(1i*phi),D,N, per);
sum_sigmajj1 = kron_matrix_generator(kron(tau',tau)*exp(-1i*theta)+...
kron(tau,tau')*exp(1i*theta),D,N,per);

sum_sigmaj= kron_matrix_generator(tau+tau',D,N,per); 

%% Exact Diagonalization

samplepos = zeros(max(size(flist))*sample_size,1); % will have sample_size
% number of samples for each change in the variable h.

E_gap = zeros(max(size(flist)),1);

lattice_sigma = zeros(N,max(size(flist)));
avg_tau= zeros(1,max(size(flist)));
avg_sigma= zeros(1,max(size(flist)));
groundstates = zeros(length(chi),max(size(flist)));
energy=zeros(1,max(size(flist)));

for n=1:(max(size(flist)))
    f = flist(n);
    J=(1-f); %negative here or theta=phi=pi/3 should make it antiferromagnetic like
    
    nevals = 4;

    H= -f*(sum_tauj) - J*(sum_sigmajj1)+Jsig*sum_sigmaj;
    [evecs,evals] = eigs(H+H', nevals, 'sr'); % make sure it is H+H' (need H.c.)
    fprintf(['Min. Energy for f=' num2str(f) ': ' num2str(evals(1,1)) '\n'])
    
    % evals=real(evals); % will be anyway as real observable
    energy(n)=min(diag(evals));
    
    diagvals = diag(evals);
    E_gap(n)= diagvals(2)-diagvals(1);
    
    GS=evecs(:,1);
    
    % Calculate entanglement entropy/spectrum for the subsystem made of first N2 spins
    N2=floor(N/2);
    V0=reshape(GS,D^N2,D^(N-N2));
    ES=eig(V0*V0'); %ES=eigs(V0*V0',10,'la');
    ES=sort(nonzeros(ES),'descend'); %Entanglement spectrum (sorted and remove zeros)
    EE(tt,n)=-sum(ES.*log(ES)); %Entanglement entropy
    
    groundstates(:,n)=GS;
    GS_probs = abs(GS.*1).^2; % probabilities of each state
    
    lattice_sigma(:,n) = ((abs(GS.*1).^2).')*sigmas; 
    avg_sigma(n) = (GS')*sum_sigmaj*GS;
    avg_tau(n)= (GS')*sum_tauj*GS;
    
    % Adding the probabilities so that each state can be represented by a
    % unique probability between 0 and 1. 
    for m=2:D^N
        GS_probs(m) = GS_probs(m)+GS_probs(m-1);
    end

    % reducing to subspace
    if exist('subspace','var')
        GS_probs=GS_probs(subspace);
    end
    % Sampling %
    for ii=1:sample_size
        a = rand;
        samplepos((n-1)*sample_size+ii,1) = sum(a>=GS_probs)+1;
        % finds the pos where a is greater than the lowest prob. This is the 
        % sample we choose by picking a random number a.
    end
    
end
samplepos_total(tt,:)=samplepos; % keeping track of all the samples
% fprintf('\n\n Sigma GS Expectation values \n')
% fprintf('h= %.2f : %.3f \n', [flist',sum(lattice_sigma,1)']')

savedir0= ['Plots/Z_' num2str(clocknum) '/' lattice per '/N_' num2str(N) '/PhaseMap'];
mkdir(savedir0)

% figure, plot(flist,energy,'-o'), title('Energy vs. f with J=1-f')
% savefig([savedir0 '/Energy_' per '_' lattice ...
%     '_N_' num2str(N) '_theta_' num2str(theta) '_phi_' num2str(phi) '.fig'])
% 
% figure, plot(flist,E_gap,'-o'), title('Energy Gap vs. f with J=1-f')
% savefig([savedir0 '/Energy_Gap_' per '_' lattice ...
%     '_N_' num2str(N) '_theta_' num2str(theta) '_phi_' num2str(phi) '.fig'])
% 
% figure, plot(flist,avg_tau,'-o'), title('Avg Tau vs. f with J=1-f')
% savefig([savedir0 '/Avg_tau_' per '_' lattice ...
%     '_N_' num2str(N) '_theta_' num2str(theta) '_phi_' num2str(phi) '.fig'])
% 
% figure, plot(flist,avg_sigma,'-o'), title('Avg Sigma vs. f with J=1-f')
% savefig([savedir0 '/Sigma_' per '_' lattice ...
%     '_N_' num2str(N) '_theta_' num2str(theta) '_phi_' num2str(phi) '.fig'])

% figure; G=graph(adjmat);
% p=plot(G);
% G.Nodes.NodeColors = angle(lattice_sigma(:,1));
% p.NodeCData = G.Nodes.NodeColors;
% colorbar
% savefig([savedir0 '/Sigma_' per '_' lattice ...
%     '_N_' num2str(N) '_theta_' num2str(theta) '_phi_' num2str(phi) '.fig'])

%% Metropolis Hastings Sampling of ED groundstate
if MH
burnin=200; chain_size=sample_size; hilbspace=omega_list;
eff_size=chain_size/num_starts;
rot= 2*pi/clocknum; % how much the samples are rotated

MHchains=zeros([max(size(flist))*chain_size,N]); 
temp_chain=zeros([num_starts*eff_size,N]);
for ii=1:max(size(flist))
    GS=groundstates(:,ii);
    psi = {GS,sigmas}; % have to combine the samples like this, input omegas
    % so that MH alg can identify which coefficient based on the sample.

    for ll=1:num_starts % adding extra burn ins/starting points (10 to be exact)
        options={'psi',psi,'angle',rot}; 
        burnchain=MetropolisHastings(burnin,hilbspace,options{:});

        options={'s0', burnchain(end,:),'psi',psi,'angle',rot};
        temp_chain((ll-1)*eff_size+1:(ll*eff_size),:)=...
            MetropolisHastings(eff_size,hilbspace,options{:});
    end 
    MHchains(((ii-1)*chain_size+1):(ii*chain_size),:)=temp_chain(:,:);
end
end
%% Diffusion Map

if MH
    samples=MHchains;
    savestr='MH';
else
    samples = sigmas(samplepos,:); % generate all of the samples
    savestr='';
end

nevals=30;

% labels=zeros([1,size(samples,1)]);
% for ii=1:max(size(flist))
%     labels(((ii-1)*sample_size+1):((ii)*sample_size))=flist(ii);
% end
if ~isempty(eps_epsstar_list)
poolobj=parpool(5); % allocates N cpu workers
degen_tol=10^(-2.5); 
nfiles=max(size(flist));
% Ys=zeros(sample_size,sample_size,nfiles);

eps_epsstar_list=(sample_size/(2*pi))*eps_list;
% ndegenerate=zeros([max(size(eps_epsstar_list)),nfiles]);

n_comp=5; alpha=1; % so that diffusion map is A_ll' 
parfor ee=1:max(size(eps_epsstar_list))

%     skip=5;count=1; % figure
    for ff =1:nfiles % applying diffusion maps to each Jx

        S=samples(((ff-1)*sample_size+1):((ff)*sample_size),:);
        % randomize samples and labels
        k=randperm(max(size(S)));
        S=S(k,:);
        % labels=labels(k);

        epsilon= eps_epsstar_list(ee)*(2*pi/sample_size);
        [vect2, vals2,Y,LL]= diffusionmaps(S',epsilon,alpha,n_comp);
        
        K=sum(abs(ones(size(vals2))-vals2)<degen_tol);
        ndegenerate(ee,tt,ff)=K;
%         ndegenerate(ee,ff)=K;

%         [idx,centroids]=kmeans(vect2,K, 'Replicates',10,'Start','plus'); 
        % ,'EmpytAction','error') %returns error if a cluster has no points
        % plus is the default and uses a k-means++ algorithm (starts the centroids
        % with probabilities proportional to the distance from other centroids
        % (more complex than that)). 'sample', 'uniform', 'cluster' other options
        % calculating the average within cluster standard deviation
        
%         avg_cluster_dev=0; avg_intercluster_dist=0;
%         for k=1:K
            % vect2(idx==k,:) returns the samples belonging to cluster k, then, the 
            % sqrt of the average of the distance squared is the standard deviation.
%             cluster_dev = sqrt(mean((vect2(idx==k,:)-centroids(k,:)).^2,1));
%             avg_cluster_dev=avg_cluster_dev+mean(cluster_dev);

            % Calculating the average intercluster distance
%             for k2=1:K
%                 Dcl=sqrt(sum((centroids(k,:)-centroids(k2,:)).^2));
%                 avg_intercluster_dist=avg_intercluster_dist+Dcl;
%             end
%         end
%         
%         clusterstdv_interclusterdist_ratio(ee,tt,ff)= 2*K*(K-1)*avg_cluster_dev/avg_intercluster_dist;
        
%         if mod(kk,skip)==1
%             subplot(ceil(max(size(flist))/(2*skip)),floor(max(size(flist))/(2*skip)),count)
%             scatter3(vect2(:,2),vect2(:,3),vect2(:,4))
%             count=count+1;
%             title(['f=' num2str(flist(kk))])
%             xlabel('vector2'), ylabel('vector3'), zlabel('vector3')
%         end
        
        % vect2(:,2).'*S; %this performs the mapping from the reduced
        %%feature space to the physical lattice space 

    end
     
%     sgtitle(['Diffusion Map evecs 2, 3, and 4 for eps/eps*='  num2str(eps_epsstar)])
%     figname=[save_dir '/Evec_Projections_N_' num2str(N) '_epsepsstar_'...
%     num2str(eps_epsstar) '_theta_' num2str(theta) '_phi_' num2str(phi)];
%     savefig([figname '.fig']), saveas(gcf,figname,'png')
    
%     figure, plot(flist, ndegenerate,'-*'), xlabel('f'),ylabel('# of evals=1')
%     title(['# degenerate max evals vs transverse field strength f for \epsilon/\epsilon* =' num2str(eps_epsstar)])
%     figname=[save_dir '/numdegenevals_' per '_N_' num2str(N) '_epsepsstar_' num2str(eps_epsstar) ...
%         '_theta_' num2str(theta) '_phi_' num2str(phi) savestr];
%     saveas(gcf,figname,'png')
%     savefig([figname '.fig'])

%     figure, plot(flist, dens_list, '-o','DisplayName','density of K')
%     hold on, plot(flist, disp_list,'-*','DisplayName','disparity of K')
%     hold on, plot(flist, cc_list,'-^','DisplayName','clustering coef. of K')
%     legend('show'),title('network measures on distance matrix K')
%     savefig([save_dir '/net_measures_' per '_N_' num2str(N) '_epsepsstar_' ...
%         num2str(eps_epsstar) '_theta_' num2str(theta) '_phi_' num2str(phi) savestr '.fig'])

%     fprintf(['Completion percentage=' num2str(100*(ee/max(size(eps_epsstar_list)))) '\n'])
end % for eps_epsstar loop
end % for if statement if eps_epsstar is empty
delete(gcp('nocreate'))
%delete(poolobj)
toc

%% PCA
K=3; dpoint_per_v=sample_size;
options = {'kmeans_dim', 3, 'symm', 'spherical','no_plot',true};
% [PCAkpred3(tt,:),~]=PCA_kmeans(imag(samples),dpoint_per_v,K, flist, adjmat, savedir0, options{:});
[PCAkpred2(tt,:),PCA_pred(tt,:)]=PCA_kmeans(imag(samples),dpoint_per_v,2, flist, adjmat, savedir0, options{:});
%% Trying the Autoencoder on this model 

% K=3; dpoint_per_v=sample_size; latent_dim=2;
% options = {'apply_transferfunc', true, 'symm', 'circular','no_plot',true};
% % autoenckpred3(tt,:)=Autoenc_kmeans(imag(samples),dpoint_per_v,latent_dim,K, flist, adjmat, savedir0, options{:});
% [autoenckpred2(tt,:),autoenc_pred(tt,:)]=Autoenc_kmeans(imag(samples),dpoint_per_v,latent_dim,2, flist, adjmat, savedir0, options{:});

%% KS testing

rot_samples=omega.*samples;
samples(round(samples,4)==round(omega,4))=0;
samples(round(samples,4)==round(omega^2,4))=-1; 
rot_samples(round(rot_samples,4)==round(omega,4))=0;
rot_samples(round(rot_samples,4)==round(omega^2,4))=-1; 
rot_samplepos=getR(real(rot_samples),1); samplepos=getR(real(samples),1);

options={'statistic','Cramer','plot',false};
[ks_p(tt,:),ks_stats(tt,:)]=symm_kstesting(samplepos,rot_samplepos,flist,1,...
    '\frac{2i\pi}{3}\ rotation', ['f\ and\ \theta=' theta ',\ \phi=' theta], options{:});

fprintf(['Completion percentage=' num2str(100*(tt/length(theta_list))) '\n'])

end % for theta loop

% save('tau_sampled_clock_L10.mat')

map_plot(real(EE.'),flist,'EE')

%% plotting the results
for ee=1:max(size(eps_epsstar_list))
    % transposed and flipped to match data presentation of paper
    figname=[savedir0 '/numdegen_phase_map' per '_N_' num2str(N) '_epsepsstar_'...
        num2str(eps_epsstar_list(ee)) savestr];
    map_plot(squeeze(ndegenerate(ee,:,:)).',flist,['Degeneracy\ for\ \frac{\epsilon}{\epsilon^*}='...
        num2str(eps_epsstar_list(ee))],figname)

%     figname=[savedir0 '/cluster_measures_ratio_' per '_N_' num2str(N) '_epsepsstar_'...
%         num2str(eps_epsstar_list(ee)) savestr];
%     map_plot(squeeze(clusterstdv_interclusterdist_ratio(ee,:,:)).',flist,...
%         [' \frac{2n \sigma}{D}\ ratio\ for\ \frac{\epsilon}{\epsilon^*}=' num2str(eps_epsstar_list(ee))],figname)
end

figname=[savedir0 '/KS_pval_phase_map_N_' num2str(N)];
map_plot(ks_p.',flist,'phi=theta\ K-S\ P\ Value',figname)
figname=[savedir0 '/Cramer_phase_map_N_' num2str(N)];
map_plot(ks_stats.',flist,'phi=theta\ Cramer\ Statistic',figname)

figname=[savedir0 '/PCA_phase_map_N_' num2str(N)];
map_plot(PCA_pred.',flist,'phi=theta\ PCA\ Phase\ Results',figname)
figname=[savedir0 '/Autoenc_phase_map_N_' num2str(N)];
map_plot(autoenc_pred.',flist,'phi=theta\ Autoencoder\ Phase\ Results',figname)

% one method to get rid of cluster=1 vs 2 assignment
for ii=1:size(PCAkpred2,1)
PCAkpred2(ii,:)=PCAkpred2(ii,:)-(max(PCAkpred2(ii,:))+min(PCAkpred2(ii,:)))/2;
PCAkpred3(ii,:)=PCAkpred3(ii,:)-(max(PCAkpred3(ii,:))+min(PCAkpred3(ii,:)))/2;
autoenckpred2(ii,:)=autoenckpred2(ii,:)-mean(autoenckpred2(ii,:));
autoenckpred3(ii,:)=autoenckpred3(ii,:)-(max(autoenckpred3(ii,:))+min(autoenckpred3(ii,:)))/2;
end

for ii=1:size(PCAkpred2,1)
    if PCAkpred2(ii,end)>0
        PCAkpred2(ii,:)=-PCAkpred2(ii,:);
    end
    if PCAkpred3(ii,end)>0
        PCAkpred3(ii,:)=-PCAkpred3(ii,:);
    end
    if autoenckpred2(ii,end)>0
        autoenckpred2(ii,:)=-autoenckpred2(ii,:);
    end
    if autoenckpred3(ii,end)>0
        autoenckpred3(ii,:)=-autoenckpred3(ii,:);
    end
end

figname=[savedir0 '/PCAk2means_phase_map_N_' num2str(N)];
map_plot(PCAkpred2.',flist,'phi=theta\ PCA\ K=2\ Phase\ Prediction',figname)
figname=[savedir0 '/PCAk3means_phase_map_N_' num2str(N)];
map_plot(PCAkpred3.',flist,'phi=theta\ PCA\ K=3\ Phase\ Prediction',figname)
figname=[savedir0 '/Autoenck2means_phase_map_N_' num2str(N)];
map_plot(autoenckpred2.',flist,'phi=theta\ Autoencoder\ K=2\ Phase\ Prediction',figname)
figname=[savedir0 '/Autoenck3means_phase_map_N_' num2str(N)];
map_plot(autoenckpred3.',flist,'phi=theta\ Autoencoder\ K=3\ Phase\ Prediction',figname)

pause

%% Training the Autoencoder and PCA on the full Phase Map

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
load("N12_phase_map.mat");
nfiles=length(flist)*length(theta_list);

if size(samplepos_total,2)>size(samplepos_total,1) % means should be flipped for (:) to access correctly
    samplepos_total=samplepos_total.';
end
samples=sigmas(samplepos_total(:),:);

% calculating and plotting the number of unique samples per parameter as
% well as other observables
nunique=zeros([nfiles,1]);
sigma_exp=zeros([nfiles,1]);
sigma_sq=zeros([nfiles,1]);
corr=zeros([nfiles,1]);
for n=1:nfiles
    nunique(n)=length(unique(samplepos_total(((n-1)*sample_size+1):(n*sample_size))));
    sigma_exp(n)=mean(samples(((n-1)*sample_size+1):(n*sample_size)));
    sigma_sq(n)= mean(samples(((n-1)*sample_size+1):(n*sample_size)).^2);
    corr(n)=mean(samples((((n-1)*sample_size+1):(n*sample_size)),1).*...
        samples((((n-1)*sample_size+1):(n*sample_size)),round(N/2)));
end
std_dev=sqrt(sigma_sq-sigma_exp.^2);

std_dev= reshape(std_dev,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
figname=[savedir0 '/imag_Sigma_std_phase_map_N_' num2str(N)];
map_plot(imag(std_dev),flist,'\phi=\theta\ imag(\tau\ Std)',figname)

sigma_exp= reshape(sigma_exp,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
figname=[savedir0 '/imag_Sigma_exp_phase_map_N_' num2str(N)];
map_plot(imag(sigma_exp),flist,'\phi=\theta\ imag(\tau\ expectation\ value)',figname)

corr= reshape(corr,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
figname=[savedir0 '/imag_Sigma_corr_phase_map_N_' num2str(N)];
map_plot(imag(corr),flist,'\phi=\theta\ imag(\tau\ 1,\ to\ \frac{L}{2}\ Correlation)',figname)

nunique= reshape(nunique,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
figname=[savedir0 '/Numunique_phase_map_N_' num2str(N)];
map_plot(nunique,flist,'\phi=\theta\ Number\ of\ unique\ samples',figname)

% applying KS/Cramer
omega=exp(2i*pi/3); % Testing how it does with a different neq 2pi/3 rotation
rot_samples=omega.*samples; altsamples=samples;
altsamples(round(altsamples,4)==round(omega,4))=0;
altsamples(round(altsamples,4)==round(omega^2,4))=-1; 
rot_samples(round(rot_samples,4)==round(omega,4))=0;
rot_samples(round(rot_samples,4)==round(omega^2,4))=-1; 

options={'statistic','Cramer','plot',false};
[ks_p,ks_stats]=symm_kstesting(getR(real(altsamples),1),getR(real(rot_samples),1),1:nfiles,1,...
    '\frac{i\pi}{4}\ rotation', ['f\ and\ \theta=' theta ',\ \phi=' theta], options{:});

ks_stats= reshape(ks_stats,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
figname=[savedir0 '/Cramer_phase_map_N_' num2str(N) '_pi2_rot'];
map_plot(ks_stats,flist,'\phi=\theta\ Cramer\ Statistic',figname)

tic % timing PCA
K=3;
options = {'kmeans_dim', 3, 'symm', 'spherical','no_plot',false}; dpoint_per_v=sample_size;
[K2PCA_pred,PCA_pred]=PCA_kmeans(imag(samples),dpoint_per_v,K, ...
    1:(size(samples,1)/sample_size), adjmat, savedir0, options{:});
toc

K2PCA_pred= reshape(K2PCA_pred,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
PCA_pred= reshape(PCA_pred,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);

% PCA trained on entire phase map results
figname=[savedir0 '/PCA_k' num2str(K) 'means_phase_map_N_' num2str(N)];
map_plot(K2PCA_pred,flist,['\phi=\theta\ PCA\ K=' num2str(K) '\ Phase\ Prediction'],figname)
figname=[savedir0 '/PCA_phase_map_N_' num2str(N)];
map_plot(abs(PCA_pred),flist,'\phi=\theta\ PCA\ Phase\ Results',figname)

%% Rerun diffusion maps
% eps_list=logspace(-1.4,-0.8,20); 
eps_list=10^(-1.0);
degen_tolerance = 10^(-2.5);
% eps_list=logspace(-0.8,-0.2,40); % if unnormalized by L in Diff maps (more general)
options = {'degen_tol',degen_tolerance,'plot',false,'alpha',1};
[ndegenerate,secondeval,evals,Ys]=diffmap_list(samples,1:nfiles,eps_list,options{:});

% eps_list=eps_epsstar_list*(2*pi/sample_size);
% for ee=1:length(eps_list)
%     Ndegen= reshape(squeeze(ndegenerate(ee,:,:)).',[max(size(uniquetol(flist,1e-4))),...
%     max(size(uniquetol(theta_list,1e-4)))]);
%     figname=[savedir0 '/Ndegen_epsilon_' num2str(eps_list(ee))];
%     map_plot(Ndegen,flist,'\phi=\theta\ Degeneracy',figname)
% end

save("N12_phase_map_Ys.mat","Ys")

%% Checking the dispersion/standard deviation of the distance matrix 
if ~exist('Ys','var')
    load("N12_phase_map_Ys.mat","Ys")
end

if ~exist('distances','var')
uniqs=uniquetol(unique(Ys),1e-4); 
distances=zeros([length(uniqs),nfiles]); 
for k=1:length(uniqs)
    distances(k,:)=squeeze(sum(sum(abs(Ys-uniqs(k))<1e-4)));
end
end

std_dev=zeros([nfiles,1]);
for kk=1:nfiles 
%     std_dev(kk)=std(distances(:,kk));
%     d_hist=histcounts(Ys(:,:,kk));
%     std_dev(kk)=std(nonzeros(d_hist)/sum(nonzeros(d_hist)));
%     y=reshape(Ys(:,:,kk),[sample_size*sample_size,1]);
%     pd=fitdist(y,'Normal');
%     std_dev(kk)=pd.sigma;
    std_dev(kk)=pdf_std(distances(:,kk)/sum(distances(:,kk)), uniqs.');
%     std_dev(kk)= sqrt(mean(y.^2)-mean(y)^2);
%     std_dev(kk)=std(Ys(:,:,kk),1,'all');
end

std_dev= reshape(std_dev,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);

figname=[savedir0 '/Clock3_distance_dist_stddev_N_' num2str(N)];
map_plot(std_dev,flist,'\phi=\theta\ Distance\ Distribution\ Stardard\ Deviation',figname)
   
%% autoencoder
tic 
K=3;
latent_dim=3; options = {'apply_transferfunc', true, 'symm', 'circular','no_plot',true};
[K2autoenc,autoenc_pred,autoenc_pred2,autoenc]=Autoenc_kmeans(imag(samples),dpoint_per_v,latent_dim,...
    K, 1:(size(samples,1)/sample_size), adjmat, savedir0, options{:});
toc

K2autoenc= reshape(K2autoenc,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
autoenc_pred=reshape(autoenc_pred,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
autoenc_pred2=reshape(autoenc_pred2,[max(size(uniquetol(flist,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);

% Autoenc trained on entire phase map results
figname=[savedir0 '/Autoenc_k' num2str(K) 'means_phase_map_N_' num2str(N)];
map_plot(K2autoenc,flist,['\phi=\theta\ Autoencoder\ K=' num2str(K) '\ Phase\ Prediction'],figname)
figname=[savedir0 '/Autoenc_phase_map_N_' num2str(N)];
map_plot(abs(autoenc_pred),flist,'\phi=\theta\ Autoencoder\ Phase\ Projection\ 1\ Results',figname)
figname=[savedir0 '/Autoenc_proj2_phase_map_N_' num2str(N)];
map_plot(abs(autoenc_pred2),flist,'\phi=\theta\ Autoencoder\ Projection\ 2\ Phase\ Results',figname)

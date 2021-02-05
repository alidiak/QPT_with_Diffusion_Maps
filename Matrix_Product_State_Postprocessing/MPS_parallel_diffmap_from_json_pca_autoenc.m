%% Init
clc; close all; clear;

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
addpath('../Utils');

%% Input the system parameters

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));

model='J1J2'; 
% valid args include: TFIM, AF-TFIM, XXZ, CMod, clock_n3,
% clock_n3_phase_map, J1J2, 
N =100; % lattice/system size
% date='_2020-03-12'; % only use this if looking at the Clock model
date='_50resamples2';
%% Load the data from .json files (these dictionaries are way easier to work with)
% 30x faster just to load the data, also already calculated MH chain samples in 
% parallel so it is way faster.

mps_dir=['/home/alex/Documents/QML_Research/NetKet/Output/MPS_Results/'...
    model '/' model '_OUTPUTS' num2str(N) date];
% out_files=glob([mps_dir '/Clock**.json']);  % this function gets a list of all .json files in the directory
out_files=glob([mps_dir '/**.json']);  % this function gets a list of all .json files in the directory
nfiles=max(size(out_files));

%% Extracting the data so it is transparent for parallelization

% matlab parallel workers require 'transparency' which means that 
% you cannot create, delete, modify, access, or query variables if
% you do not explicitly specify these variables in the code.
dat1=jsondecode(fileread(out_files{1}));
chain_size=size(dat1.MHchain_real,1); 
MHchains=zeros([nfiles,chain_size,N]); 
v_list=zeros([nfiles,1]);theta_list=zeros([nfiles,1]);phi_list=zeros([nfiles,1]);
conv_list=zeros([1,nfiles]); bondentropy=zeros([1,nfiles]); energy_list=zeros([1,nfiles]);
ks_p=zeros([1,nfiles]); ks_stats=zeros([1,nfiles]);ks_static=zeros([1,nfiles]);
samples=zeros([nfiles*chain_size,N]);
tic
for kk=1:nfiles
    jsondat=jsondecode(fileread(out_files{kk}));
    if isfield(jsondat,'n_excited_states') 
    nexcited=jsondat.n_excited_states; 
    else 
        nexcited=0;
    end
    if nexcited>0
        % just taking the ground state (see script of same name with nexc to use excited states) 
        MHchains(kk,:,:)=jsondat.MHchain_real(:,:,1)+1i*jsondat.MHchain_imag(:,:,1); 
        energy_list(kk)=jsondat.energy(1);
    else
        MHchains(kk,:,:)=jsondat.MHchain_real+1i*jsondat.MHchain_imag;
        energy_list(kk)=jsondat.energy;
    end
    if isfield(jsondat,'restarts')
        restarts=jsondat.restarts;
    end
    v_list(kk)=jsondat.g; 
    if isfield(jsondat,'theta_')
        theta_list(kk)=jsondat.theta_;  phi_list(kk)=jsondat.phi;
        bondentropy(kk)=jsondat.BondEntropy(round(N/2));
        omega=exp(2j*pi/3); chains=squeeze(MHchains(kk,:,:));
        samples((kk-1)*chain_size+1:kk*chain_size,:)=chains;
        rot_samples=chains*omega; 
        chains(round(chains,4)==round(omega,4))=0; 
        chains(round(chains,4)==round(omega^2,4))=-1;
        rot_samples(round(rot_samples,4)==round(omega,4))=0;
        rot_samples(round(rot_samples,4)==round(omega^2,4))=-1; 
        rot_samplepos=getR(real(rot_samples),1); samplepos=getR(real(chains),1);
        [~,ks_p(kk),ks_static(kk)]=kstest2(samplepos,rot_samplepos);options={'Cramer',true};
        ks_stats(kk)=kuipertest2(samplepos,rot_samplepos,2000,false,options{:});
    end
    conv_list(kk)=jsondat.Converged;
end
toc

% Making sure the output is organized correctly (glob/list of files is
% weird with putting files like 0.9 behind 0.9* (ex. 0.925)
[v_list,ind] = sort(v_list); 
bondentropy=bondentropy(ind); theta_list=theta_list(ind);
energy_list=energy_list(ind); conv_list=conv_list(ind); 
MHchains=MHchains(ind,:,:); phi_list=phi_list(ind);
ks_p=ks_p(ind); ks_stats=ks_stats(ind); ks_static=ks_static(ind);

ntheta=numel(uniquetol(theta_list,1e-4));
for jj=1:numel(v_list)/ntheta
    [theta_list((jj-1)*ntheta+1 : jj*ntheta),I]=sort(theta_list((jj-1)*ntheta+1 : jj*ntheta));
    chng=bondentropy((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    bondentropy((jj-1)*ntheta+1 : jj*ntheta)=chng;
    chng=energy_list((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    energy_list((jj-1)*ntheta+1 : jj*ntheta)=chng;
    chng=MHchains((jj-1)*ntheta+1 : jj*ntheta,:,:); chng=chng(I,:,:);
    MHchains((jj-1)*ntheta+1 : jj*ntheta,:,:)=chng;
    chng=phi_list((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    phi_list((jj-1)*ntheta+1 : jj*ntheta)=chng;
    chng=conv_list((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    conv_list((jj-1)*ntheta+1 : jj*ntheta)=chng;
    chng=ks_stats((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    ks_stats((jj-1)*ntheta+1 : jj*ntheta)=chng;
    chng=ks_p((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    ks_p((jj-1)*ntheta+1 : jj*ntheta)=chng;
    chng=ks_static((jj-1)*ntheta+1 : jj*ntheta); chng=chng(I);
    ks_static((jj-1)*ntheta+1 : jj*ntheta)=chng;
end

save_dir=[model '/plots/N_' num2str(N) '/' date];
mkdir(save_dir)

%% Applying k-means PCA and autoencoder on the full phase diagram

% adjacency matrix
v=ones((N-1),1); adjmat = diag(v,1)+diag(v,-1); per=''; lattice='1dchain';
if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
end

tic % timing PCA
options = {'kmeans_dim', 3, 'symm', 'spherical','no_plot',true}; dpoint_per_v=chain_size;
[K2PCA_pred,PCA_pred]=PCA_kmeans(imag(samples),dpoint_per_v,2, 1:nfiles, adjmat, save_dir, options{:});
toc

K2PCA_pred= reshape(K2PCA_pred,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
PCA_pred= reshape(PCA_pred,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
% PCA on entire phase plot results
figname=[save_dir '/PCA_k2means_phase_map_N_' num2str(N)];
map_plot(K2PCA_pred.',v_list,'phi=theta\ PCA\ K=2\ Phase\ Prediction',figname)
figname=[save_dir '/PCA_phase_map_N_' num2str(N)];
map_plot(abs(PCA_pred).',v_list,'phi=theta\ PCA\ Phase\ Results',figname)

tic % timing autoencoder
latent_dim=2; options = {'apply_transferfunc', true, 'symm', 'circular','no_plot',true};
[K2autoenc,autoenc_pred,autoenc_pred2]=Autoenc_kmeans(imag(samples),dpoint_per_v,latent_dim,...
    2, 1:nfiles, adjmat, save_dir, options{:});
toc

K2autoenc= reshape(K2autoenc,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
autoenc_pred=reshape(autoenc_pred,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
autoenc_pred2=reshape(autoenc_pred2,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);

% Autoenc on entire phase results
figname=[save_dir '/Autoenc_k2means_phase_map_N_' num2str(N)];
map_plot(abs(K2autoenc).',v_list,'phi=theta\ Autoencoder\ K=2\ Phase\ Prediction',figname)
figname=[save_dir '/Autoenc_phase_map_N_' num2str(N)];
map_plot(abs(autoenc_pred).',v_list,'phi=theta\ Autoencoder\ Phase\ Projection\ 1\ Results',figname)
figname=[save_dir '/Autoenc_proj2_phase_map_N_' num2str(N)];
map_plot(abs(autoenc_pred).',v_list,'phi=theta\ Autoencoder\ Projection\ 2\ Phase\ Results',figname)



%% parallel computation of number of Diffusion Map degeneracies
chain_size=500;
% eps_epsstar_list=1.0:0.25:5.75; % this range is generally good for samples=500
eps_list=logspace(-1.4,-1.2,20); eps_epsstar_list=eps_list/(2*pi/chain_size);
% eps_epsstar_list=1.0:0.5:12.5; 
eps_list=eps_epsstar_list*(2*pi/chain_size);
ndegenerate=zeros([max(size(eps_epsstar_list)),nfiles]);
n_comp=5;  % this is almost always the maximum degeneracy 
nreplicas=10; % number of distinct cluster initiations to use (picks the single best fit)
clusterstdv_interclusterdist_ratio=zeros([max(size(eps_epsstar_list)),nfiles]);
degen_tol=10^(-2.5); 
% diffmap_vecs=zeros([max(size(eps_epsstar_list)),nfiles,250,3end]); % 3 because plotting in max 3 dims
diffmap_vecs=zeros([max(size(eps_epsstar_list)),nfiles,chain_size,3]); % 3 because plotting in max 3 dims
poolobj=parpool(4); % allocates N cpu workers
parfor ee=1:max(size(eps_epsstar_list))% does this in parallelK
    
    tic 
    for kk =1:nfiles % applying diffusion maps to each change in variable        
        %S=MHchains(((kk-1)*chain_size+1):((kk)*chain_size),:);
        %S=jsondat.MHchain_real+1i*jsondat.MHchain_imag;
        S=squeeze(MHchains(kk,1:500,:));
        % randomize samples
        k=randperm(max(size(S)));
        S=S(k,:);

        epsilon=eps_epsstar_list(ee)*(2*pi/chain_size);
%         epsilon=eps_list(ee);
        alpha=1; 
        [vect2, vals2]= diffusionmaps(0.5.*S',epsilon,alpha,n_comp);

        
        K=sum(abs(ones(size(vals2))-vals2)<degen_tol);
        ndegenerate(ee,kk)=K;
        diffmap_vecs(ee,kk,:,:)=[vect2(:,1) vect2(:,2) vect2(:,3)];
        
%         [idx,centroids]=kmeans(vect2,K, 'Replicates',nreplicas,'Start','plus'); 
        % ,'EmpytAction','error') %returns error if a cluster has no points
        % plus is the default and uses a k-means++ algorithm (starts the centroids
        % with probabilities proportional to the distance from other centroids
        % (more complex than that)). 'sample', 'uniform', 'cluster' other options
        % calculating the average within cluster standard deviation
        
%         avg_cluster_dev=0; avg_intercluster_dist=0;
%         for k=1:K
%             % vect2(idx==k,:) returns the samples belonging to cluster k, then, the 
%             % sqrt of the average of the distance squared is the standard deviation.
%             cluster_dev = sqrt(mean((vect2(idx==k,:)-centroids(k,:)).^2,1));
%             avg_cluster_dev=avg_cluster_dev+mean(cluster_dev);
% 
%             % Calculating the average intercluster distance
%             for k2=1:K
%                 D=sqrt(sum((centroids(k,:)-centroids(k2,:)).^2));
%                 avg_intercluster_dist=avg_intercluster_dist+D;
%             end
%         end
        
%         clusterstdv_interclusterdist_ratio(ee,kk)= 2*K*(K-1)*avg_cluster_dev/avg_intercluster_dist;
        %cluster_dev_mat(ee,kk)=avg_cluster_dev/K;
%         if K==1
%             interclust_dist_mat(ee,kk)=0;
%         else
%             interclust_dist_mat(ee,kk)=(1/(K*(K-1)))*avg_intercluster_dist;   
                %possibly 0.5 times this as distances are counted twice
%             
%         end               
    end

    toc

end

delete(gcp('nocreate'));

%%   Plotting the diffusion map projections
save_dir=[model '/plots/N_' num2str(N) '/' date];
mkdir(save_dir)

IC=find((round(v_list,3)==0.525).*(theta_list==pi/3)); % IC phase conditions
Topo=find((round(v_list,3)==0.025).*(theta_list==pi/6)); % Topological phase conditions
Trivial=find((round(v_list,3)==1.025).*(theta_list==pi/6)); % Trivial phase conditions
skip=5; count=1; 
for ee=1:max(size(eps_epsstar_list))
    if mod(ee-1,skip)==0
        subplot(ceil(max(size(eps_epsstar_list))/(skip)),3,count);
        ICvecs=squeeze(diffmap_vecs(ee,IC,:,1:3));
        scatter3(ICvecs(:,1),ICvecs(:,2),ICvecs(:,3))
        title(['IC Phase Evecs \epsilon/\epsilon *=' num2str(eps_epsstar_list(ee))])

        %subplot(ceil(max(size(eps_epsstar_list))/(2*skip)),3*floor(max(size(eps_epsstar_list))/(2*skip)),count+1);
        subplot(ceil(max(size(eps_epsstar_list))/(skip)),3,count+1); 
        Topovecs=squeeze(diffmap_vecs(ee,Topo,:,1:3));
        scatter3(Topovecs(:,1),Topovecs(:,2),Topovecs(:,3))
        title(['Topo Phase Evecs \epsilon/\epsilon *=' num2str(eps_epsstar_list(ee))])

        subplot(ceil(max(size(eps_epsstar_list))/(skip)),3,count+2);
        Trivialvecs=squeeze(diffmap_vecs(ee,Trivial,:,1:3));
        scatter3(Trivialvecs(:,1),Trivialvecs(:,2),Trivialvecs(:,3))
        title(['Trivial Phase Evecs \epsilon/\epsilon *=' num2str(eps_epsstar_list(ee))])
        
        xlabel('vector2'), ylabel('vector3'), zlabel('vector4')
        sgtitle(['Diffusion Map evecs 2, 3, and 4 for IC, Topological and Trivial phases'])

        count=count+3;
    end     
    
end
figname=[save_dir '/DiffMap_EvecProjections_N_' num2str(N)];
saveas(gca,[figname '.png'],'png'), savefig([figname '.fig'])
%% plotting the results

save_dir=[model '/plots/N_' num2str(N) '/' date];
mkdir(save_dir)

if length(v_list)==nfiles % means it's 1D data
     setfig(16,'J_2',0,'Energy',0,[model '\ MPS\ Energy'],'on'); plot(v_list,energy_list,'.-');
     figname=[save_dir '/Energy_N_' num2str(N)]; saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
     mkdir([save_dir '/diffusionmaps'])
     for ee=1:max(size(eps_epsstar_list))
        figname=[save_dir '/diffusionmaps/numde gen_phase_map_N_' num2str(N) '_epsepsstar_' ...
            num2str(eps_epsstar_list(ee))];
        setfig(16,'J_2',0,' \#\ of\ Clusters',0,[model '\ MPS\ \#\ of\ Clusters\  \epsilon=' ...
            num2str(eps_epsstar_list(ee)*(2*pi/chain_size))],'on'); 
        plot(v_list,ndegenerate(ee,:),'.-'); saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
    end
end

[vals, ia]=uniquetol(v_list,1e-4); v_grid=max(size(vals));
while numel(energy_list)<(v_grid*max(size(uniquetol(theta_list,1e-4)))-1)
%     if numel(energy_list)==(v_grid*max(size(unique(theta_list))))
%         break
%     end
for ii=1:v_grid-1
    if ia(ii+1)-ia(ii)~=v_grid
        kk=1;
        for theta=unique(theta_list).'
            if sum(theta==theta_list(ia(ii):(ia(ii+1)-1)))==0 % finds missing theta
                bondentropy=[bondentropy(1:(ia(ii)+kk)) bondentropy(ia(ii)+kk-1) bondentropy((ia(ii)+kk)+1:end)];
                energy_list=[energy_list(1:(ia(ii)+kk)) energy_list(ia(ii)+kk-1) energy_list((ia(ii)+kk)+1:end)];
                conv_list=[conv_list(1:(ia(ii)+kk)) 0 conv_list((ia(ii)+kk)+1:end)];
                ndegenerate=[ndegenerate(:,1:(ia(ii)+kk)) ndegenerate(:,ia(ii)+kk-1) ndegenerate(:,(ia(ii)+kk+1):end)];
                clusterstdv_interclusterdist_ratio=[clusterstdv_interclusterdist_ratio(:,1:(ia(ii)+kk)) ...
                    clusterstdv_interclusterdist_ratio(:,ia(ii)+kk-1) clusterstdv_interclusterdist_ratio(:,(ia(ii)+kk+1):end)];
            end
        end
        [vals, ia]=unique(v_list); 
    end
end
end

bondentropy(bondentropy<0)=0; bondentropy(bondentropy>1.6)=0;% fixing errored values
bondentropy=reshape(bondentropy,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);

ks_p=reshape(ks_p,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
ks_stats=reshape(ks_stats,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
ks_static=reshape(ks_static,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
energy_list=reshape(energy_list,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);
conv_list=reshape(conv_list,[max(size(uniquetol(theta_list,1e-4))),...
    max(size(uniquetol(v_list,1e-4)))]);

if max(size(unique(phi_list)))>1
    str='\theta';
else
    str=num2str(uniquetol(phi_list,1e-3));
end

% Plotting the central cut bond entropy
figname=[save_dir '/BondEE_phase_map_N_' num2str(N)];
map_plot(bondentropy.',v_list,['phi=' str '\ central\ cut\ Bond\ Entropy'],figname)

% Plotting the energy
figname=[save_dir '/Energy_phase_map_N_' num2str(N)];
map_plot(energy_list.',v_list,['phi=' str '\ Energy'],figname)

% Plotting the ks pval
figname=[save_dir '/ks_p_phase_map_N_' num2str(N)];
map_plot(ks_p.',v_list,['phi=' str '\ KS\ P-value'],figname)
% Plotting the ks stat
figname=[save_dir '/ks_stat_phase_map_N_' num2str(N)];
map_plot(ks_static.',v_list,['phi=' str '\ KS\ Statistic'],figname)
% Plotting the Cramer stat
figname=[save_dir '/Cramer_stats_phase_map_N_' num2str(N)];
map_plot(ks_stats.',v_list,['phi=' str '\ Cramer-von-Mises\ Statistic'],figname)


% Plotting the number of degenerate diff map evals and kmeans cluster results
for ee=1:max(size(eps_epsstar_list))
    % Plotting Ndegeneratefigname=[save_dir '/Convlist_phase_map_N_' num2str(N)];
    ndeg=reshape(ndegenerate(ee,:),[max(size(uniquetol(theta_list,1e-4))),max(size(uniquetol(v_list,1e-4)))]);
    figname=[save_dir '/numdegen_phase_map_N_' num2str(N) '_epsepsstar_'...
        num2str(eps_epsstar_list(ee))];
    map_plot(squeeze(ndeg.'),v_list, ['phi=' str '\ degeneracy\ for\ \epsilon='... %['phi=' str '\ degeneracy\ for\ \frac{\epsilon}{\epsilon^*}='...
        num2str(eps_epsstar_list(ee)*(2*pi/chain_size))],figname)

    % Plotting ratio of avg cluster standard deviation to avg intercluster distance
%     ratio=reshape(clusterstdv_interclusterdist_ratio(ee,:),[max(size(uniquetol(theta_list,1e-4))),...
%         max(size(uniquetol(v_list,1e-4)))]);
%     figname=[save_dir '/cluster_measures_ratio_N_' num2str(N) '_epsepsstar_'...
%         num2str(eps_epsstar_list(ee))];
%     map_plot(flip(squeeze(ratio.')),v_list,['phi=' str ' 2*n*\sigma/D ratio for \epsilon/\epsilon*='...
%         num2str(eps_epsstar_list(ee))],figname)
end

% Plotting the convergence list
figname=[save_dir '/Convlist_phase_map_N_' num2str(N)];
map_plot(conv_list.',v_list,['\phi = ' str '\ Convergence\ List'],figname)

%% PCA and autoenc sample prep (taking a slice of theta) 

theta='0.942477796077'; % tried 0.942477796077, 0.785398163397, 0.837758040957
% Just going to look at an interesting slice for now
save_dir=[model '/plots/N_' num2str(N) '/' date];
mps_dir=['/home/alex/Documents/QML_Research/NetKet/Output/MPS_Results/'...
    model '/' model '_OUTPUTS' num2str(N) date];
out_files=glob([mps_dir '/theta_' theta '/**.json']); nfiles=max(size(out_files));
dat1=jsondecode(fileread(out_files{1})); chain_size=size(dat1.MHchain_real,1);
MHchains=zeros([nfiles,chain_size,N]);
v_list2=zeros([nfiles,1]);theta_list2=zeros([nfiles,1]);phi_list2=zeros([nfiles,1]);
for kk=1:nfiles
    jsondat=jsondecode(fileread(out_files{kk}));
    if isfield(jsondat,'n_excited_states') 
    nexcited=jsondat.n_excited_states; 
    else 
        nexcited=0;
    end
    if nexcited>0
        MHchains(kk,:,:)=jsondat.MHchain_real(:,:,1)+1i*jsondat.MHchain_imag(:,:,1); 
    else
        MHchains(kk,:,:)=jsondat.MHchain_real+1i*jsondat.MHchain_imag;
    end
    if isfield(jsondat,'restarts')
        restarts=jsondat.restarts;
    end
    v_list2(kk)=jsondat.g; theta_list2(kk)=jsondat.theta_;  phi_list2(kk)=jsondat.phi;
end

[v_list2,ind] = sort(v_list2); MHchains=MHchains(ind,:,:);

chains=zeros([length(v_list2)*chain_size,N]);
for jj=1:length(v_list2)
    chains((jj-1)*chain_size+1:jj*chain_size,:)=squeeze(MHchains(jj,:,:));
end

% adjacency matrix
v=ones((N-1),1); adjmat = diag(v,1)+diag(v,-1);
per=''; lattice='1dchain';
if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
end
omega=exp(2j*pi/3); % applying Z3 symmetry by rotating the samples 
rot_samples=omega.*chains; %rotating

%% Testing the number of unique samples

altsamples=samples;
altsamples(round(altsamples,4)==round(omega,4))=0;
altsamples(round(altsamples,4)==round(omega^2,4))=-1; 

sample_size=size(MHchains,2);
samplepos_total=getR(real(altsamples),1);

nunique=zeros([nfiles,1]);
for n=1:nfiles
    nunique(n)=length(unique(samplepos_total(((n-1)*sample_size+1):(n*sample_size))));
end

nunique= reshape(nunique,[  max(size(uniquetol(v_list,1e-4))),...
    max(size(uniquetol(theta_list,1e-4)))]);
figname=[save_dir '/Numunique_phase_map_N_' num2str(N)];
map_plot(nunique.',v_list,'\phi=\theta\ Number\ of\ unique\ samples',figname)

        
%% PCA on slice

% Have to map the samples to real values or only use imaginary
% chains(round(chains,4)==round(omega,4))=0; 
% chains(round(chains,4)==round(omega^2,4))=-1; chains=real(chains);

K=3; dpoint_per_v=500;
options = {'kmeans_dim', 3, 'symm', 'spherical'};
kpred=PCA_kmeans(imag(chains),dpoint_per_v,K, v_list2, adjmat, save_dir, options{:});

%% Trying the Autoencoder on slice

K=3; dpoint_per_v=500; latent_dim=2;
options = {'apply_transferfunc', true, 'symm', 'circular'};
kpred=Autoenc_kmeans(imag(chains),dpoint_per_v,latent_dim,K, v_list2, adjmat, save_dir, options{:});

%% Testing probability comparison method K-S on symmetry transformations on slice

% Changing the sigma basis to a spin=1 (-1,0,1) basis for the linear index mapping
chains(round(chains,4)==round(omega,4))=0; chains(round(chains,4)==round(omega^2,4))=-1; 
rot_samples(round(rot_samples,4)==round(omega,4))=0;
rot_samples(round(rot_samples,4)==round(omega^2,4))=3-1; 
rot_samplepos=getR(real(rot_samples),1); samplepos=getR(real(chains),1);

options={'statistic','Cramer'};
[ks_tests,ks_stats]=symm_kstesting(samplepos,rot_samplepos,v_list2,1,'\frac{2i\pi}{3}\ rotation',...
     ['f\ and\ \theta=' theta ',\ \phi=' theta], options{:});
 
 %% Specific for J1J2 processing
  
 if strcmp(model,'J1J2')
    plot_all_pt=0.48; plot_all_end=0.53; skip=4;     
    kk=1;  lines=[]; 
    setfig(16,'Log_{10}(\epsilon)',0,'Degeneracy',0,'Exponential\ Scaling\ of\ Degeneracy','on');
    while kk<=length(v_list)
        kk
        if kk==find(abs(v_list-plot_all_pt)<1e-4)
            skprev=skip; skip=1;
        end
        if kk==find(abs(v_list-plot_all_end)<1e-4)
            skip=skprev;
        end
        if mod(kk,skip)==0
            lines(end+1)=plot(log10(eps_list), (ndegenerate(:,kk)),'.-','DisplayName',['J_2=' num2str(v_list(kk))]);
            hold on;
        end
         kk=kk+1;
    end
    c_map=parula(length(lines)); %c_map(find(abs(v_list-0.5)<1e-4)/skip,:)=[1,0,0]; % makes J_2=0.5 red
    c_map(length(lines)/2,:)=[1,0,0]; % assuming J_2=0.5 is in the center
    set(lines,{'color'},num2cell(c_map,2));
    legend('FontSize',8)
    figname=[save_dir '/ndegen_exponentialscaling_N_' num2str(N)];
    saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
 end
 
 %% Averaging the diffusion map results over sample sets and doing error analysis

n_runs=size(MHchains,2)/dat1.chain_size;  
sample_size=size(MHchains,2)/n_runs;
 
degen_tolerance = 10^(-2.5);
eps_list=10^(-1.3);% choose the value of epsilon to do error analysis for
ndegen_tot=zeros([n_runs,length(v_list)]);

parpool(4)
parfor run=1:n_runs
    tic

    S=MHchains(:,(run-1)*sample_size+1:(run*sample_size),:);
    samples=squeeze(S(1,:,:));
    for ii=2:length(v_list)
        samples=[samples ; squeeze(S(ii,:,:))];
    end
    options = {'degen_tol',degen_tolerance,'save_dir', save_dir, 'fig_name', 'Ndegen_Evals','alpha',1,...
    'fig_title', 'degenerate\ evals=1\ vs\ J_2\','var_name','J_2'};
    [ndegenerate,~,~,~]=diffmap_list(0.5.*samples,v_list,eps_list,options{:});
    % 0.5 for the distance normalization in this case
    toc

    fprintf(['run ' num2str(run) ' finished of  ' num2str(n_runs) '\n'])
    
    ndegen_tot(run,:)=squeeze(ndegenerate);

end

delete(gcp('nocreate'));

std_err=sqrt(sum((ndegen_tot-mean(ndegen_tot,1)).^2,1)/n_runs)/sqrt(n_runs);

figure; plot(v_list,std(ndegen_tot,1)/sqrt(n_runs),'.-')

fig1=setfig(16,'J_2',0,'Degeneracy',0,['Max\ Eigenvalue\ Degeneracy\ for\ \epsilon=' num2str(eps_list)],'on');
e=errorbar(v_list,mean(ndegen_tot,1),std_err,'.-');


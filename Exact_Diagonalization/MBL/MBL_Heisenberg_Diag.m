%% Init
clear ; close all; clc

path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));
addpath('../../Utils');

%% Creating an adjacency matrix for a simple 1d lattice

N=4; % choose system size
spin=1/2; % choose spin
% KK=0; % choose sz subspace to work in
num_seeds=20; % choose the number of random seeds/runs to average over
ensemble='long_time'; 
state='neel'; % pick neel for neel order state, leave blank if want first half up, second half down initial state
finite_size=false;
degen_tolerance = 1e-4; % default is 1e-4
%ensemble='diag';

v=ones((N-1),1); v2=ones((N-2),1);
adjmat = diag(v,1)+diag(v,-1);

% periodic conditions
per='periodic';
lattice='1dchain';

if strcmp(per,'periodic')
    adjmat(1,end)=1; adjmat(end,1)=1;
end

G=graph(adjmat);% p=plot(G);

degen_per_h=true; MH=false;

% warning('off','all')

%% Exact Diagonalization Matriceswarning('off','all')

Z=getZ(N,spin); Dim=2*spin+1; %Get the matrix elements for all {Z_m} which are diagonal

% This is the nearest neighbor heisenberg interaction Hamiltonian
if exist('KK','var')
[H1,chi]=Heisenberg_Hamiltonian_gen(N, adjmat,Z,spin,KK); % only enter K=0 if N is even
else
 [H1,chi]=Heisenberg_Hamiltonian_gen(N, adjmat,Z,spin); 
end

DK=size(H1,1); % the reduced matrix size
D=Dim^N;
ia=(1:D); % Linear index array
ja= (1:DK)';
xi=zeros(D,1); xi(chi)=ja; %xi maps the indices from original to projected space
Zred=Z(chi,:); % Z in the reduced subspace

%Now constructing the total Sx,Sy,Sz operator (we can easily square them without constructing two-point correlators)
Sz=sparse(ja,ja,sum(Zred,2),DK,DK,DK); % should be all 0 sparse if K=0

%% Finding Energies and Creating Samples

J= 1.0;  % nearest neighbor Heisenberg interaction

dh = 10;
hlist=[0.01:dh:300.01];
% hlist =[0.0:dh:(2-dh) 2:(dh/2):4 (4+dh):dh:10];
% hlist=[10];

sample_size=500;
deps=2*3.1250e-04;
% eps_list=0.0075:deps:0.020; % normal method seems to work better at lower eps
% eps_list=logspace(-2.2,-1.2,30);
eps_list=[0.0089]; 

savedir=['plots/' ensemble '_ensemble/' lattice per '/' state '/tol_' num2str(degen_tolerance)...
    '/N_' num2str(N) ]; mkdir(savedir)

E_gap_avg = zeros(max(size(hlist)),1);
energy_avg=zeros(max(size(hlist)),1);
energy_density=zeros(max(size(hlist)),1);
ndegenerate_avg=zeros(length(eps_list),max(size(hlist)));
secondeval_avg=zeros([length(eps_list),length(hlist)]);
eval_avg_avg=zeros([length(eps_list),length(hlist)]);
Kpred_avg=zeros(max(size(hlist)),1);
PCA_pred_avg=zeros(max(size(hlist)),1);
nunique_avg=zeros([length(hlist),1]);
init_state_overlap_avg=zeros([1,length(hlist)]);

% Histogramming
dist_to_initial_state = zeros([sample_size*length(hlist),1]);
distance_avgs=zeros([N/2+1,length(hlist)]);

% for saving and recreating purposes (so I don't have to rerun the full ED again)
wavefunctions=zeros([num_seeds, length(hlist), DK]);
samplepos_tot=zeros([num_seeds, length(hlist)*sample_size]);

nevals = length(chi);

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
        
warning('off','MATLAB:MKDIR:DirectoryExists')

fprintf('Beginning  run %.2f out of %.3f \n', [seed,num_seeds]')

samplepos = zeros(max(size(hlist))*sample_size,1); % will have sample_size

E_gap = zeros(max(size(hlist)),1);

lattice_Sz = zeros(N,max(size(hlist)));
lattice_Sx = zeros(N,max(size(hlist)));
sigmax = zeros(1,max(size(hlist)));
states = zeros(max(size(hlist)),DK);
energy=zeros(1,max(size(hlist)));
init_state_overlap=zeros(1,max(size(hlist)));
e_dens=zeros(1,length(hlist));

for n=1:(max(size(hlist)))
    
    h = hlist(n); % here h is the range of the uniformly distributed random variables on Sz
    h_lattice= -h + (2*h)*rand(N,1); % array of uniformly drawn random #s applied to lattice sites
    
    Sz_h=sparse(ja,ja,Zred*h_lattice,DK,DK,DK);
    H=(J*H1+Sz_h);
    
%     Sx=spalloc(D,D,N*D);
%     Sy=spalloc(D,D,N*D);
%     for m=1:N
%         Sx=Sx+sparse(ia,ia+(Dim*Z(:,m)*2^(N-m)).',h_lattice(m)*1,D,D,D);
% %         Sy=Sy+sparse(ia,ia+(Dim*Z(:,m)*2^(N-m)).',h_lattice(m)*1i*Z(:,m),D,D,D);
%     end
    
%     H=(J*H1+(Sx));%(Sx+Sy)*0.5);  
    
%     nevals = 2;
    [evecs,evals] = eigs(H, DK, 'SA');
    fprintf(['Min. Energy for h=' num2str(h) ': ' num2str(evals(1,1)) '\n'])
    
    energy(n)=evals(1,1);
    
    e_dens(n)=(init_psi'*(H)*init_psi-min(diag(evals)))/(max(diag(evals))-min(diag(evals)));
    
    % energy gap calc
    % uniq=uniquetol(diag(evals));
    diagvals = diag(evals);
    E_gap(n)= diagvals(2)-diagvals(1);
    
%     state=evecs(:,1); % 1 is GS
%     state_probs = abs(state).^2; % probabilities of each state
    
    % This simulation assumes infinite temperature which means all states
    % are equally probable, therefore, for each sample we will take not
    % only a random state from the eigenvector, but also a random eigenvector.
    
    % For the diagonal ensemble, define an initial state - at t->inf this state will still be prevalent 
%     init_psi=zeros([1,length(chi)]); init_psi(2) =1; %init_psi(end)=1; % cat-state, all spin up/down 

%     init_psi=init_psi/norm(init_psi);
    
    
    % The alternative long time evolution of initial psi
    if strcmp(ensemble,'long_time')
    t=10000;
%     psi_E  = expm(1i*H*t)*init_psi; % note that there are many far more efficient ways to do this: 
       %  use diagonalized H or treat as differential equation
       
    psi_E  = (evecs*diag(exp(1i*diagvals*t))/evecs)*init_psi;  % using the already diagonalized H will be a bit faster and memory efficient
       
    elseif strcmp(ensemble,'diag')
    % Project the initial state onto the energy basis (i.e. get all the
    % coefficients c_0, c_1,... c_(chi). 
        psi_E=evecs*init_psi;
    end
    
%     tot=0;  % check the recreated state
%     for ii=1:length(chi)
%     tot=tot+psi_E(ii)*evecs(:,ii);
%     end
    
    states(n,:)=psi_E;

    init_state_overlap(n)= abs(init_psi'*psi_E).^2; % state overlap (measure of how MBL it is?)
    state_probs=abs(psi_E).^2; % same as diag(psi_E.' * psi_E)

    % equivalent way
%     p=init_psi.'*init_psi;  % init state projector
%     p_tild=evecs.'*p*evecs; % projection onto E 'basis'
%     state_probs=diag(p_tild); % diag ensemble probabilities
    
    for m=2:DK
        state_probs(m) = state_probs(m)+state_probs(m-1);
    end
    for ii=1:sample_size
        
    % In the MBL paper, they only look at the middle 1/3 of states
%     E_state= randi(round(DK/3));
%     state_list(ii,n)=round(DK/3)+E_state; % keeping track of which states were used
  
    % Adding the probabilities so that each state can be represented by a
    % unique probability between 0 and 1. 
    
    % Sampling %
    a = rand;
    samplepos((n-1)*sample_size+ii,1) = sum(a>=state_probs)+1;
    % finds the pos where a is greater than the lowest prob. This is the 
    % sample we choose by picking a random number a.
    end
end

wavefunctions(seed,:,:)=states;
samplepos_tot(seed,:)=samplepos;

% Plotting
% fprintf('\n\n SigmaZ GS Expectation values \n')
% fprintf('h= %.2f : %.3f \n', [(hlist)',sum(lattice_Sz,1)']')

% 
% figure, plot(hlist,energy,'-o'), title('Energy vs. h')
% figname=[savedir '/Energy_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% 
% figure, plot(hlist,E_gap,'-o'), title('Energy Gap vs. h')
% figname=[savedir '/Energy_Gap_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% 
% figure, plot(hlist,dimer_corr,'-o'), title(['Dimer Correlation vs. h for r=' num2str(r)])
% figname=[savedir '/Dimer_Corr_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% 
% % figure, plot(hlist,orderpar/(N*spin^2),'-o'), title('Order Param vs. h')
% % ylabel('1/(LS^2) \Sigma_l (-1)^l \sigma_l^z \sigma_{l+1}^z'); xlabel('h');
% % figname=['plots/' lattice per '/N_' num2str(N) '/Orderparam_' per '_' lattice '_N_' num2str(N)];
% % saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
% samplepos_tot(seed,:)=samplepos;
% figure, plot(hlist,sum(lattice_Sz,1)','-o'), title('Avg SigmaZ vs. h')
% figname=['plots/' lattice per '/N_' num2str(N) '/SigmaZ_' per '_' lattice '_N_' num2str(N)];
% saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% Checking the found states
% h_lattice 
% gs=Z(find(abs(abs(evecs(:,1)).^2-1)<1e-1),:) % ground state majority
% long_time_state=Z(find(max(abs(psi_E).^2)==abs(psi_E).^2),:) 
% sum(psi_E.'*Sz*psi_E)

samples=Zred(samplepos,:); 

%% Testing probability comparison methods on symmetry transformations

xlate_by=1; options={'statistic','Cramer','plot',false};
xlated_samples=[samples(:,end-(xlate_by-1):end)  samples(:,1:end-xlate_by)];
[ks_tests,ks_stats]=symm_kstesting(getR(samples/2,1/2),getR((xlated_samples)/2,1/2),hlist,1/2,...
    ['translation\ by\ ' num2str(xlate_by)], 'J_2',options{:}); %'Bz\ of\ \Sigma_i i\sigma_i^z ',options{:});

%% Distance Measure 2
% eps_list=0.01:0.01:0.02; % this alternate dist seems to work best at higher eps
% options = {'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',1,...
% 'fig_title', 'number of degenerate max evals vs h','dist','staggered'};
% ndegenerate=diffmap_list(samples,hlist,eps_list,options{:});

%% Original Distance measure (Euclidian)
options = {'degen_tol',degen_tolerance,'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',1,...
'fig_title', 'number of degenerate max evals vs h'};
[ndegenerate,secondeval,eval_avg, Ys]=diffmap_list(samples,hlist,eps_list,options{:});

if length(unique(Ys))~=(N/2+1)
    fprintf('Number of unique distances different than predicted')
end

% for histogramming the distance 
distances=zeros([N/2+1,length(hlist)]); i=1;
for k=0:2:N
    distances(i,:)=squeeze(sum(sum(abs(2*N*Ys-k)<1e-4)));
    i=i+1;
end
distance_avgs=distance_avgs+distances;

%% number unique
nunique=zeros([length(hlist),1]);
for n=1:length(hlist)
    nunique(n)=length(unique(samplepos(((n-1)*sample_size+1):(n*sample_size))));
end
figure('Visible','off'); plot(hlist,nunique,'k-*'), ylabel('# unique samples'), xlabel('J_2') 
title('unique samples per J_2'); 
figname=[savedir '/number_unique_samples_' per '_' lattice '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
%% PCA

K=2; dpoint_per_v=sample_size;
options = {'kmeans_dim', 3, 'symm', 'spherical','no_plot',true};
[Kpred,PCA_pred]=PCA_kmeans(samples,dpoint_per_v,K, hlist, adjmat, savedir, options{:});

%% averaging all the important vals

E_gap_avg = E_gap_avg+E_gap;
energy_avg=energy_avg+energy;
energy_density=energy_density+e_dens;
ndegenerate_avg=ndegenerate_avg+ndegenerate;
secondeval_avg=secondeval_avg+secondeval;
eval_avg_avg=eval_avg_avg+eval_avg;
Kpred_avg=Kpred_avg+Kpred;
PCA_pred_avg=PCA_pred_avg+PCA_pred;
nunique_avg=nunique_avg+nunique;
init_state_overlap_avg=init_state_overlap_avg+init_state_overlap;
dist_to_initial_state=dist_to_initial_state+sum((mix-samples).^2,2);
% if spacings
% level_spacing_avg=level_spacing_avg+level_spacings;
% end
fprintf('Run # %.2f out of %.3f done \n', [seed,num_seeds]')

end
 
delete(gcp('nocreate'));

save(['plots/' ensemble '_ensemble/N_' num2str(N) state '_SX.mat']);

%% resampling 
sample_size=500;

samplepos_tot = zeros(num_seeds, max(size(hlist))*sample_size); % will have sample_size
states=zeros(length(hlist),DK);
tic
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
toc

%% Loading previous run
path = matlab.desktop.editor.getActiveFilename;
cd(fileparts(path));

N=14; state='neel'; % _local
% load(['plots/long_time_ensemble/N_' num2str(N) state '.mat']);
load(['plots/long_time_ensemble/N_14neel_large_h.mat']);


%% Rerunning diffusion maps on the sample data with different parameters
degen_tolerance=10^(-2.5); alpha=1; n_comp=20;
eps_list=logspace(-1.9,-1.52,20);
% eps_list=10^(-1.75);
eps_list=0.014;

savedir=['plots/' ensemble '_ensemble/' lattice per '/' state '/tol_' num2str(degen_tolerance)...
    '/N_' num2str(N) 'alpha_' num2str(alpha) ]; mkdir(savedir)
options = {'degen_tol',degen_tolerance,'save_dir', savedir, 'fig_name', 'Ndegen_Evals','alpha',alpha,...
'fig_title', 'number of degenerate max evals vs h','ncomp',n_comp};
ndegenerate_avg=zeros(length(eps_list),max(size(hlist))); secondeval_avg=zeros([length(eps_list),length(hlist)]);
distance_avgs=zeros([N/2+1,length(hlist)]);
% evals_tot_avg=zeros([length(eps_list),length(hlist), sample_size]);
% Ls_avg=zeros([sample_size,sample_size,length(hlist),length(eps_list)]);
Ys_avg=zeros([sample_size,sample_size,length(hlist)]);
std_avg=zeros([1,length(hlist)]);
Ls_dist_avg=zeros([N/2+1,length(hlist),length(eps_list)]);

if length(eps_list)==1
ndegen_tot=zeros(num_seeds,length(hlist));
end

% diffmap_vecs_avg=zeros([sample_size,n_comp,length(hlist),length(eps_list)]);
% parpool(4)
for seed=1:num_seeds
    warning('off','MATLAB:MKDIR:DirectoryExists')
    
    samples=Zred(samplepos_tot(seed,:),:);
    
    fprintf(['run ' num2str(seed) ' started \n'])
    
%     tic
    [ndegenerate,secondeval,evals,Ys,~,~]=diffmap_list(samples,hlist,eps_list,options{:});
%     evals,Ys,Ls,diffmap_vecs
%     toc
    
    ndegenerate_avg=ndegenerate_avg+ndegenerate; 
    secondeval_avg=secondeval_avg+secondeval;
    
    if length(eps_list)==1
        ndegen_tot(seed,:)=squeeze(ndegenerate);
    end
    
%     for histogramming the distance (and avging Ls as a func of distance)
%     uniqs=uniquetol(unique(Ys),1e-4); 
%     distances=zeros([length(uniqs),length(hlist)]);
%     for k=1:length(uniqs)
%         distances(k,:)=squeeze(sum(sum(abs(Ys-uniqs(k))<1e-4)));
%     end
    distances=zeros([N/2+1,length(hlist)]); dists=(0:2:N)/N;
    Ls_d=zeros([N/2+1,length(hlist),length(eps_list)]);
    i=1;
    for k=0:2:N
      bool_dist=abs(N*Ys-k)<1e-4;
      distances(i,:)=squeeze(sum(sum(bool_dist)));
      % Averaging the term Ls as a function of distance 
%       Ls_d2=zeros([N/2+1,length(hlist)]);
%       for ee=1:length(eps_list)
%           for hh=1:length(hlist)
%                 Ls2=squeeze(Ls(:,:,hh,ee));
%                 Ls_d(i,hh,ee)=mean(mean(Ls2(bool_dist(:,:,hh)),1),2);
%           end
%       end
      i=i+1;
    end
    distance_avgs=distance_avgs+distances;
%     Ls_dist_avg=Ls_dist_avg+Ls_d;
%     Ls_avg=Ls_avg+Ls;
    
    for kk=1:length(hlist) % only works if not parfor loop
%         y=reshape(Ys(:,:,kk),[sample_size*sample_size,1]);
%         pd=fitdist(y,'Normal');
%         std_avg(kk)=std_avg(kk)+pd.sigma;
%         std_avg(kk)= std_avg(kk)+sqrt(mean(y.^2)-mean(y)^2);
%         std_avg(kk)=std_avg(kk)+std(distances(:,kk)/sum(distances(:,kk)));
%           std_avg(kk)=std_avg(kk)+pdf_std(distances(:,kk)/sum(distances(:,kk)), dists.');
        std_avg(kk)=std_avg(kk)+std(Ys(:,:,kk),1,'all');
    end
    Ys_avg=Ys_avg+Ys;
    
%     diffmap_vecs_avg=diffmap_vecs_avg+diffmap_vecs;
%     evals_tot_avg=evals_tot_avg+evals;
    
    fprintf(['run ' num2str(seed) ' completed \n'])
    
end
delete(gcp('nocreate'));

% Standard deviation or Dispersion calculation
std_dev=zeros([length(hlist),1]);
for kk=1:length(hlist)
%     std_dev(kk)=std(distance_avgs(:,kk)/num_seeds);
    std_dev(kk)= std(N*Ys_avg(:,:,kk)/num_seeds,0,'all');
end
setfig(16,'h',0,'Standard\ Deviation',0,'Distance\ Distributions\ Standard\ Deviation','on'); 
plot(hlist,std_avg/num_seeds,'.-')
% plot(hlist,std_dev,'.-')
figname=[savedir '/Distance_distribution_MBL_Heis_std_dev_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% save('N_10_evals_tot_avg.mat','evals_tot_avg')
% save('N_10_Ls_avg.mat','Ls_avg')
% load('')
%% Graphing connectivity

G=graph(Ys(:,:,2)); figure; plot(G); title(['Distance Graph h=' num2str(hlist(2))]);
G2=graph(Ys(:,:,end-1)); figure; plot(G2); title(['Distance Graph h=' num2str(hlist(end-1))]);

%%
% h_ind=1;
% K=exp(-Ys_avg(:,:,h_ind)/eps_list(end)); 
% Dn=sum(K,2);
% L=K./(Dn*Dn.');
% 
%  distances=zeros([N/2+1,length(hlist)]); Ls_d=zeros([N/2+1,1]);
% i=1;
% for k=0:2:N
%   distances(i)=squeeze(sum(sum(abs(2*N*Ys(:,:,h_ind)-k)<1e-4)));
%   % Averaging the term Ls as a function of distance 
%   Ls_d(i)=squeeze(mean(mean((abs(2*N*Ys(:,:,h_ind)-k)<1e-4).*L,1),2)); 
%   i=i+1;
% end
% 
% figure;  plot(0:2:N,Ls_d(:),'.-')
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

setfig(16,'Log_{10}(\epsilon)',0,'Log_{10}(Degeneracy)',0,'Exponential\ Scaling\ of\ Degeneracy','on');
lines=[]; skip=2;
for kk=1:skip:length(hlist)
        lines(end+1)=plot(log10(eps_list(:)), log10(ndegenerate_avg(:,kk)/num_seeds),'.-','MarkerSize',10,'DisplayName',['\delta J=' num2str(hlist(kk))]);
        hold on;
end
c_map=parula(length(lines)); %c_map(find(abs(v_list-0.5)<1e-4)/skip,:)=[1,0,0]; % makes J_2=0.5 red
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/ndegen_exponentialscaling_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
hold off;

% Converting to and plotting as 3D
figure; [X,Y] = meshgrid(hlist,eps_list);
surf(X,Y,log10(ndegenerate_avg/num_seeds));
xlabel('h'); ylabel('\epsilon'); zlabel('Log_{10}(Degeneracy)')
figname=[savedir '/3D_Degeneracy_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

figure; surf(X,Y,1-secondeval_avg/num_seeds,'FaceAlpha',0.75);
xlabel('h'); ylabel('\epsilon'); zlabel('Log_{10}(Degeneracy)')
figname=[savedir '/3D_1_minus_lambda_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% eps=repmat(eps_list,1,length(hlist)); hs=repmat(hlist,length(eps_list),1);
% figure; scatter3(log10(eps),hs(:),log10(ndegenerate_avg(:)/num_seeds))

%%% Other stuff we were checkin

% samples=Zred(samplepos_tot,:);
% samples(1:2:end)=-samples(1:2:end);
% 
% avg_dist=[];
% staggered_sz=[];
% for kk=1:length(hlist)
%     staggered_sz(end+1)=mean(samples((kk-1)*sample_size*num_seeds+1:kk*sample_size*num_seeds));
%     avg_dist(end+1)=mean(dist_to_initial_state((kk-1)*sample_size+1:kk*sample_size)./20);
% end
% 
% figure; plot(hlist,abs(staggered_sz),'o-')
% figure; plot(hlist,abs(avg_dist),'o-')


%% Plotting the average L_ij to better understand the effect of normalization

epsnum=1;
setfig(16,'Distance,\ d',0,'Avg.\ L_{ij}=d',0,['Average\ L_{ij}=d\ vs\ d\ \epsilon=' ...
    num2str(eps_list(epsnum)) '\ \alpha=' num2str(alpha)],'on');
% average of L_ij wrt the distance for different disorders h
for hnum=1:20:41
plot(0:2:N,log10(Ls_dist_avg(:,hnum,epsnum)./num_seeds...
    ),'.-','DisplayName',['h=' num2str(hlist(hnum))])
hold on
end
legend; figname=[savedir '/Average_N_' num2str(N) '_eps_' num2str(eps_list(epsnum))];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

setfig(16,'Distance,\ d',0,'log\ of\ Effective\ Probability',0,['Log_{10}(\frac{L_{ij}=d}{N_d}) vs\ d\ \epsilon=' ...
    num2str(eps_list(epsnum)) '\ \alpha=' num2str(alpha)],'on');
% average of L_ij wrt the distance for different disorders h
for hnum=1:20:41
plot(0:2:N,log10(Ls_dist_avg(:,hnum,epsnum)./distance_avgs(:,hnum)),'.-','DisplayName',['h=' num2str(hlist(hnum))])
hold on
end
legend; figname=[savedir '/Effective_prob_LvNd_' num2str(N) '_eps_' num2str(eps_list(epsnum))];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Plot the full degeneracy spectrum, net meaures on L, diffusion map projections, and do slope analysis

% Full degeneracy spectrum
num_plot=200;
skip=10; [X,Y] = meshgrid(1:num_plot,hlist);
for ee=1:skip:length(eps_list)
    figure; surf(X,Y, squeeze(evals_tot_avg(ee,:,1:num_plot)))
    xlabel('n^{th} eval'); ylabel('h'); zlabel('Eigenvalues'); title(['Eigenvalue Decay \epsilon=' num2str(eps_list(ee))]);
    figname=[savedir '/3D_Eval_decay_' num2str(N) '_eps_' num2str(eps_list(ee))];
    saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])
end
% Not super helpful at face value/ with this representation. 


%% Plotting averages

% saving the data from the run for potential finite size scaling 

% if spacings
% figure; plot(hlist(4:end),level_spacing_avg(4:end).'/num_seeds,'k.-'); title('Level Spacings Ratio Avg')
% saveas(gcf,[savedir '/Level_Spacings_Ratio.png'],'png')
% end
% savedir=['plots/' ensemble '_ensemble/' lattice per '/' state '/N_' num2str(N) ]; mkdir(savedir)
figure; plot(hlist,nunique_avg.'/num_seeds,'k.-'); title('Number of Unique Samples Avg')
saveas(gcf,[savedir '/avg_unique_samples.png'],'png')

for ee =1:length(eps_list)
    epsilon=eps_list(ee);
    eps_epsstar=(sample_size/(2*pi))*epsilon;
    fig1=setfig(16,'h',0,'Degeneracy',0,['Avg\ Degeneracy\ for\ \epsilon=' num2str(epsilon)],'off'); 
    plot(hlist,ndegenerate_avg(ee,:)/num_seeds,'.-')
    figname=[savedir '/diffusionmap/degenerate_evals/Ndegen_Evals_epsepsstar_' num2str(eps_epsstar)];
    savefig(fig1,[figname '.fig']); saveas(fig1,[figname '.png'],'png')
    
%     fig2=setfig(16,'h',0,'Avg\ log(1-\lambda_2)',0,'Lambda\ measures','off');
%     subplot(2,1,1); plot(hlist,(log10(1-secondeval_avg(ee,:)/num_seeds)),'k-*'); 
%     title(['$ log10(1-\lambda_2)\ vs\ h\ for\ \epsilon=' num2str(epsilon) '$'],'interpreter','latex')
%     xlabel('h','interpreter','latex'); ylabel('log_{10}(1-\lambda_2)','interpreter','latex')
%     subplot(2,1,2); plot(hlist,eval_avg_avg(ee,:)/num_seeds,'k-*');
%     title(['$ mean(\lambda_i)\ vs\ h\ for\ \epsilon=' num2str(epsilon) '$'],'interpreter','latex')
%     xlabel('h','interpreter','latex'); ylabel('$ \Epsilon_i^N(\lambda_i)/N $','interpreter','latex')
%     figname=[savedir '/diffusionmap/evals2_avg/avg_and_2ndeval_epsepsstar_' num2str(eps_epsstar)];
%     savefig(fig2,[figname '.fig']); saveas(fig2,[figname '.png'],'png')
end

figure; plot(hlist,energy_avg.'/num_seeds,'k.-'); title('Energy avg'); saveas(gcf,[savedir '/energy.png'],'png')
figure; plot(hlist,E_gap_avg.'/num_seeds,'k.-'); title('Energy Gap'); saveas(gcf,[savedir '/E_gap.png'],'png')
figure; plot(hlist,energy_density.'/num_seeds,'k.-'); title('Avg Energy density'); saveas(gcf,[savedir '/' state '_energy_density.png'],'png')
figure; plot(hlist,Kpred_avg.'/num_seeds,'k.-'); title('K-means Prediction avg'); saveas(gcf,[savedir '/PCA_Kmeans.png'],'png')
figure; plot(hlist,PCA_pred_avg.'/num_seeds,'k.-'); title('PCA prediction avg'); saveas(gcf,[savedir '/PCApred.png'],'png')
if exist('init_state_overlap_avg','var')
figure; plot(hlist,init_state_overlap_avg/num_seeds,'k.-'); title('Initial State Overlap avg'); saveas(gcf,[savedir '/Init_state_overlap.png'],'png')
end
% Plotting overlayed numdegen lines from different epsilon regimes
if ~exist('deps','var')
deps=eps_list(2)-eps_list(1);
end
to_25=int8((0.0025+0*deps)/deps); 
setfig(16,'h',0,'Degeneracy',0,'Normalized\ Average\ Degeneracy\ for\ Different\ \epsilon','on');
plot(hlist, [ndegenerate_avg(1,:)./sum(ndegenerate_avg(1,:)); ndegenerate_avg(1+to_25,:)./sum(ndegenerate_avg(1+to_25,:))],'.-')...
%    ; ndegenerate_avg(1+2*to_25,:)./sum(ndegenerate_avg(1+2*to_25,:))],'.-')
% plot(hlist, ndegenerate_avg(1+to_25,:)./num_seeds,'.-')
% legend(['\epsilon=' num2str(eps_list(1+to_25))])
legend({['\epsilon=' num2str(eps_list(1))],['\epsilon=' num2str(eps_list(1+to_25))]});%,['\epsilon=' num2str(eps_list(1+2*to_25))]})
saveas(gcf,[savedir '/NdegenComb.png'],'png');savefig([savedir '/NdegenComb.fig'])

% interpolate and store the intercept points 
f1=ndegenerate_avg(1,:)./sum(ndegenerate_avg(1,:)); 
f2=ndegenerate_avg(1+to_25,:)./sum(ndegenerate_avg(1+to_25,:));
intercept_pt=find((f1-f2)<0); % should be intercept point
y1=f1(intercept_pt-1); y2=f2(intercept_pt-1); 
m1=f1(intercept_pt)-f1(intercept_pt-1); m2=f2(intercepot_pt)-f2(intercept_pt-1);
x_intercept=(y2-y1)/(m1-m2);

% if ~finite_size
%     pause
% end


%% Histogram of the distance of the samples to the initial state (and just distances present in general)

figure=setfig(16,'h',0,'States\ Equal\ Distance\ d',0,['Inter-Sample\ Distances;\ State\ ' state],'on'); 
for k=1:(N/2+1)
    plot(hlist,distance_avgs(k,:)/num_seeds,'-o','DisplayName',['d=' num2str((k-1)*2) ])
    hold on
end
legend; figname=[savedir '/distances_vs_h_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

figure=setfig(16,'Distance',0,'States\ Equal\ Distance\ d',0,['Inter-Sample\ Distances;\ State\ ' state],'on'); 
for kk=1:6:length(hlist)
    plot(0:2/N:1,distance_avgs(:,kk)/num_seeds,'-o','DisplayName',['h=' num2str(hlist(kk)) ])
    hold on
end
legend(); figname=[savedir '/Distance_distribution_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

% Standard deviation or Dispersion calculation
std_dev=zeros([length(hlist),1]);
for kk=1:length(hlist)
%     std_dev(kk)=std(distance_avgs(:,kk)/num_seeds);
    std_dev(kk)= std(N*Ys_avg(:,:,kk)/num_seeds,0,'all');
end
figure=setfig(16,'h',0,'Standard\ Deviation',0,'Distance\ Distributions\ Standard\ Deviation','on'); 
plot(hlist,std_avg/num_seeds,'.-')
% plot(hlist,std_dev,'.-')
figname=[savedir '/Distance_distribution_MBL_Heis_std_dev_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

if ~exist('dist_to_initial_state','var')
dist_to_initial_state = zeros([sample_size*length(hlist),1]);
for seed=1:num_seeds
%     for h=1:length(hlist)
    samples=Zred(samplepos_tot(seed,:),:);
%     end
    dist_to_initial_state=dist_to_initial_state+sum((mix-samples).^2,2);
end
end

figure=setfig(16,'Average\ Distance',0,'Counts',0,['Sample\ Distances\ to\ Initial\ State:\ ' state],'on'); 
for kk=1:1:length(hlist)
    histogram(dist_to_initial_state((kk-1)*sample_size+1:kk*sample_size)./20,'FaceAlpha',0.5,...
      'DisplayName',['h=' num2str(hlist(kk))]);
    hold on
end
legend(); figname=[savedir '/Hist_dist_to_init_state_' state '_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

%% Finite sized scaling processing
ensemble= 'long_time';  % ensemble to do processing for
dir= ['plots/' ensemble '_ensemble/'];
file_list=glob([dir '**.mat']);

% fig1=setfig(16,'h',0,'Degeneracy',0,'','on'); 
% fig2=setfig(16,'h',0,'Avg\ log(1-\lambda_2)',0,'Lambda\ measures','on');
fig3=setfig(16,'h',0,'Avg\ Samples',0,'Number\ of\ Unique\ Samples\ Avg','on');
N_list=[];
for ff = 1:length(file_list)
    file=file_list{ff}
    load(file); % load the data
    
    N_list(end+1)=N;
    if ff==1
        ndegenerate_avg_n=zeros([size(ndegenerate_avg),length(file_list)]);
        secondeval_avg_n=zeros([size(secondeval_avg),length(file_list)]);
    end
    
    ndegenerate_avg_n(:,:,ff)=ndegenerate_avg(:,:);
    secondeval_avg_n(:,:,ff)=secondeval_avg(:,:);
       
    fig3;   plot(hlist,nunique_avg.'/sum(nunique_avg),'.-','DisplayName',num2str(N)); hold on;  legend
end

mkdir([dir 'Finite_size/evals2_avg'])
mkdir([dir 'Finite_size/degenerate_evals'])
for ee =1:length(eps_list)
    fig1=setfig(16,'h',0,'Degeneracy',0,['Normalized\ Average\ Degeneracy\ for\ \epsilon=' num2str(epsilon)],'off'); 
    epsilon=eps_list(ee);
    eps_epsstar=(sample_size/(2*pi))*epsilon;
    for ff=1:length(file_list)
        plot(hlist,ndegenerate_avg_n(ee,:,ff)./sum(ndegenerate_avg_n(ee,:,ff)),...
            '.-','DisplayName',['L=' num2str(N_list(ff))] ); hold on
    end
    hold off;     legend('show')
    figname=[dir 'Finite_size/degenerate_evals/Ndegen_Evals_epsepsstar_' num2str(eps_epsstar)];
%     title(['$Avg\ Degeneracy\ for\ \epsilon=' num2str(epsilon) '$'],'interpreter','latex');
    savefig(fig1,[figname '.fig']); saveas(fig1,[figname '.png'],'png')

%     fig2=setfig(16,'h',0,'log_{10}(1-\lambda_2)',0,'Lambda\ measures','off');
%     for ff=1:length(file_list)
%         plot(hlist,(log10(1-secondeval_avg_n(ee,:,ff)./mean(secondeval_avg_n(ee,:,ff)))),...
%             '-*','DisplayName',num2str(N_list(ff))); hold on
%     end
%     title(['$ log10(1-\lambda_2)\ vs\ h\ for\ \epsilon=' num2str(epsilon) '$'],'interpreter','latex')
% %     xlabel('h','interpreter','latex'); ylabel('log_{10}(1-\lambda_2)','interpreter','latex')
%     legend;  figname=[dir 'Finite_size/evals2_avg/avg_and_2ndeval_epsepsstar_' num2str(eps_epsstar)];
%     savefig(fig2,[figname '.fig']); saveas(fig2,[figname '.png'],'png')
    
end

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

n_runs=20;  
eps_list=10^(-1.6); 
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
figure; plot(hlist,std(ndegen_tot,1),'.-')

fig1=setfig(16,'J_2',0,'Degeneracy',0,['Max\ Eigenvalue\ Degeneracy\ for\ \epsilon=' num2str(eps_list)],'on');
e=errorbar(hlist,mean(ndegen_tot,1),std(ndegen_tot,1),'.-');
% plot(hlist(2:end),mean(ndegen_tot(:,2:end),1))
%% DBSCAN over a range
eps_l=linspace(1.36,1.43,20); minpts=1;

parpool(4)
num_clusters_avg=zeros([length(hlist),length(eps_l)]);
parfor seed=1:num_seeds
    warning('off','MATLAB:MKDIR:DirectoryExists')
    fprintf(['run ' num2str(seed) ' started \n'])

    samples=Zred(samplepos_tot(seed,:),:);
    
    num_clusters=zeros([length(hlist),length(eps_l)]); % N_outlier=zeros( [length(hlist),length(eps_l)]);
    for ee=1:length(eps_l)
        idx=dbscan(samples,eps_l(ee),minpts);

        for h=1:length(Jlist)
            num_clusters(h,ee)=length(unique(idx((h-1)*sample_size+1:h*sample_size))); 
%             N_outlier(h,ee)=sum(idx((h-1)*sample_size+1:h*sample_size)==-1);
        end
    end
    
    num_clusters_avg=num_clusters_avg+num_clusters;
    fprintf(['run ' num2str(seed) ' completed \n'])
    
end
delete(gcp('nocreate'));

%% plot DBSCAN vs epsilon

setfig(16,'DBSCAN\ \epsilon',0,'Number\ of\ Clusters',0,['Number\ of\ Clusters\ vs\ radius,\ minpts=' num2str(minpts)],'on');
lines=[]; skip=2;
for kk=1:skip:length(hlist)
        lines(end+1)=plot((eps_l(:)), (num_clusters_avg(kk,:)/num_seeds),'.-','MarkerSize',10,'DisplayName',['h=' num2str(hlist(kk))]);
        hold on;
end
c_map=parula(length(lines)); %c_map(find(abs(v_list-0.5)<1e-4)/skip,:)=[1,0,0]; % makes J_2=0.5 red
set(lines,{'color'},num2cell(c_map,2));
legend('FontSize',8)
figname=[savedir '/DBSCAN_Nclusters_N_' num2str(N)];
saveas(gcf,[figname '.png'],'png') ; savefig([figname '.fig'])

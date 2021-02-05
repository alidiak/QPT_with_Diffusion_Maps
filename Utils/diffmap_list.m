function [ndegen_tot,secondeval,evals,Ys,Ls,diffmap_vecs]=diffmap_list(samples,vlist,eps_list,varargin)

% nevals=30; 
sample_size=size(samples,1)/numel(vlist);
degen_tol=1e-4; str=''; dist=true;  
varname='var'; alpha=1; plot_on=true; % default
LS=false;

if sample_size<10
    n_comp=sample_size;
else
    n_comp=10;
end

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'distance'
            dist=varargin{2};
        case 'save_dir'
            save_dir=varargin{2};
        case 'nevals' 
            nevals=varargin{2}; % initial config
        case 'degen_tol'
            degen_tol=varargin{2}; % for MPS
        case 'sample_type'
            str=varargin{2};
        case 'fig_name'
            fig_name=varargin{2};
        case 'fig_title'
            fig_title=varargin{2};
        case 'var_name'
            varname=varargin{2};
        case 'dist'
            dist_type=varargin{2}; 
        case 'alpha'
            alpha=varargin{2}; % if alpha=1, uses a modified kernel that is more robust to # samples
        case 'plot'
            plot_on=varargin{2};
        case 'ncomp'
            n_comp=varargin{2};
        case 'Ls'
            LS=varargin{2};
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

ndegen_tot=zeros([length(eps_list),length(vlist)]);
secondeval=zeros([length(eps_list),length(vlist)]);
evals=zeros([length(eps_list),length(vlist),sample_size]);
eval_avg=zeros([length(eps_list),length(vlist)]);
if dist
    Ys=zeros([sample_size,sample_size,length(vlist)]); 
else
    Ys=0;
end
% Y_avg=zeros([1,length(vlist)]);
% Y_unique=zeros([1,length(vlist)]);

if LS
Ls=zeros([sample_size,sample_size,length(vlist),length(eps_list)]);
else
    Ls=0;
end
diffmap_vecs=zeros([sample_size,n_comp,length(vlist),length(eps_list)]); 

for ee=1:length(eps_list)
    ndegenerate=zeros([length(vlist),1]); N_samples=[];
    for kk =1:max(size(vlist)) % applying diffusion maps to each Jx

        S=samples(((kk-1)*sample_size+1):((kk)*sample_size),:);
        % randomize samples and labels
        k=randperm(size(S,1));
        S=S(k,:);
        % labels=labels(k);
    
        epsilon=eps_list(ee);
        
%         epsilon=eps_list(ee)*(2*pi/sample_size);
        eps_epsstar=(sample_size/(2*pi))*epsilon;
        if exist('dist_type','var')
%             options={'dist',dist_type};
%             [vect2, vals2]= diffusionmaps(S',epsilon,alpha,n_comp,options{:});
                [vect2,vals2,N_samples(end+1),K,Y]=diffusionmaps_dist2(S.',epsilon,alpha,n_comp);
        else
            if dist
                if LS
                [diffmap_vecs(:,:,kk,ee), vals2,Ys(:,:,kk),Ls(:,:,kk,ee)]= diffusionmaps(S.',epsilon,alpha,n_comp);
                else
                [diffmap_vecs(:,:,kk,ee), vals2,Ys(:,:,kk),~]= diffusionmaps(S.',epsilon,alpha,n_comp);
                end
            else
                [vect2, vals2,~]= diffusionmaps(S.',epsilon,alpha,n_comp); 
            end
%                 [vals2,vect2]=dmaps(S,epsilon,sample_size); % sample_size is always the max possible number of degenerate evals
        end
        
        ndegenerate(kk)=sum(abs(ones(size(vals2))-vals2)<degen_tol);
        secondeval(ee,kk)=vals2(2);
        eval_avg(ee,kk)=mean(vals2);
        evals(ee,kk,:)=vals2;
%         Y_avg(kk)=mean(mean(Y));
%         Y_uniq(kk)=length(unique(Y));
        
    end
        
    if exist('fig_title','var')
        Title = [fig_title ' for\ \epsilon=' num2str(epsilon)]; % \frac{\epsilon}{\epsilon^*}=' num2str(eps_epsstar)];
    else
        Title= ['Diffusion\ Map\ Eigenvalue\ Degeneracy\ vs\ ' varname ...
            '\ for\ \epsilon=' num2str(epsilon)];
            %\frac{\epsilon}{\epsilon^*}=' num2str(eps_epsstar)];
    end
    
%     fig3=setfig(16,varname,0,'Avg\ Distance\ and\ unique\ distance\',0,'Distance');
%     plot(vlist, [Y_avg(ee,:) ; Y_uniq],'-*')% xlabel('Jx'),ylabel('# of evals=1')
%     legend({'# unique distances','Avg Distance'})
    
    if plot_on
        fig1=setfig(16,varname,0,'Degeneracy',0,Title,'off');
        plot(vlist, ndegenerate,'.-'); title(['$' Title '$'], 'interpreter','latex')
        ndegen_tot(ee,:)=ndegenerate;

        fig2=setfig(16,varname,0,'log_{10}(1-\lambda_2)',0,'Lambda measures','off');
        subplot(2,1,1); %
        plot(vlist,(log10(1-secondeval(ee,:))),'.-'); 
        title(['$ log_{10}(1-\lambda_2)\ vs\ ' varname '\ for\ \epsilon=' num2str(epsilon) '$'], 'interpreter','latex')
        xlabel(['$' varname '$'],'interpreter','latex'); 
        ylabel('$ log_{10}(1-\lambda_2) $', 'interpreter','latex')
        subplot(2,1,2);
        plot(vlist,eval_avg(ee,:),'.-');
        title(['$ mean(\lambda_i)\ vs\ ' varname '\ for\ \epsilon=' num2str(epsilon) '$'], 'interpreter','latex')
        xlabel(['$' varname '$'],'interpreter','latex'); 
        ylabel('$ \sum_i^N (\lambda_i) /N $', 'interpreter','latex')
    end
    
    if exist('save_dir','var')
        save_dir2=[save_dir '/diffusionmap/'];
        mkdir([save_dir2 '/degenerate_evals'])
        mkdir([save_dir2 '/evals2_avg'])

        figname=[save_dir2 'degenerate_evals/' fig_name '_epsepsstar_' num2str(eps_epsstar) str];
        savefig(fig1,[figname '.fig']); saveas(fig1,[figname '.png'],'png')
        
        figname=[save_dir2 'evals2_avg/avg_and_2ndeval_epsepsstar_' num2str(eps_epsstar) str];
        savefig(fig2,[figname '.fig']); saveas(fig2,[figname '.png'],'png')
%         
%         figname=[save_dir2 'evals2_avg/Distance_avg_epsepsstar_' num2str(eps_epsstar) str];
%         savefig(fig3,[figname '.fig']); saveas(fig3,[figname '.png'],'png')
    end
    
end
if ~isempty(N_samples)
    figure, plot(vlist,N_samples,'.-'), xlabel('var'), ylabel('# samples used')
end

end
function [Kpred,PCA_pred]=PCA_kmeans(samples,dpoints_per_v,K, vlist, adjmat, savedir, varargin)

% symm='rotational';
kmeans_dim=3; no_show=false; % default val

while ~isempty(varargin)
    switch lower(varargin{1})
        %case 'kmeansdata' 
         %   kmeansdata=varargin{2}; % symmetry to impose
        case 'symm'
            symm=varargin{2};
        case 'kmeans_dim'
            kmeans_dim=varargin{2};
        case 'no_plot'
            no_show=varargin{2};
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

[U,S] = svd((samples')*samples);

sample_size=max(size(samples))/numel(vlist);
N=size(S,1);

if ~no_show
figure, ax1 = subplot(3,1,1); ax2 = subplot(3,1,2); ax3=subplot(3,1,3);
plot(ax1,1:N,diag(S)/norm(diag(S)),'-o')
title(ax1,['PCA Evals'])
ylabel(ax1,'Ordered Evals')

plot(ax2,1:N,U(:,1),'-*'), 
title(ax2,['Largest Evec'])
xlabel(ax2,'Lattice Site'), ylabel(ax2,'Lattice Site Evec Val.')

plot(ax3,1:N,U(:,2),'-*'), 
title(ax3,['2nd Largest Evec'])
xlabel(ax3,'Lattice Site'), ylabel(ax3,'Lattice Site Evec Val.')
figname=[savedir '/PCA_Results'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

figure, plot(vlist,abs(mean(squeeze(reshape((samples*U(:,1)),...
    [sample_size,size(vlist)])))),'-o','DisplayName',['Evec' num2str(1)])
hold on
plot(vlist,abs(mean(squeeze(reshape((samples*U(:,2)),...
    [sample_size,size(vlist)])))),'-o','DisplayName',['Evec' num2str(2)])
hold off, legend('show')
title('Data projection of evecs as a function of h')
figname=[savedir '/PCA_data_Projection'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

% Plotting it on the non-periodic latttice regardless of if it is periodic
% or not
% if ~strcmp(lattice,'1dchain')
% nonper = importdata(['../Edgemats/' lattice '_edges_' num2str(N) '_vertices.txt'],' ' ,3);
% edges = nonper.data(:,:);
% adjmat = sparse([edges(:,1) edges(:,2)],[edges(:,2) edges(:,1)],[edges(:,3) edges(:,3)],N,N);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
G = graph(adjmat);
p=plot(G);
G.Nodes.NodeColors = real(U(:,1));
p.NodeCData = G.Nodes.NodeColors;
colorbar
title('Colormap for largest evec')
figname=[savedir '/Graph_1stEvec'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

figure
p2=plot(G);
G.Nodes.NodeColors = real(U(:,2));
p2.NodeCData = G.Nodes.NodeColors;
colorbar
title('Colormap for 2nd largest evec')
figname=[savedir '/Graph_2ndEvec'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')
figure; 
end 

% dpoints_per_v=sample_size; 
skip=sample_size/dpoints_per_v; % the number of points to be skipped in samples
ncolors=numel(vlist); % ncolors=size(samples,1)/sample_size;
C=repmat(1:ncolors,sample_size/skip,1); % for mapping each array to a color
y1=samples(1:skip:end,:)*U(:,1);y2=samples(1:skip:end,:)*U(:,2);
y3=samples(1:skip:end,:)*U(:,3); % projections onto PCA subspace
if ~no_show
s=repmat(5*ones(ncolors,1),sample_size/skip,1); scatter3(y1,y2,y3,s(:),C(:)) 
xlabel('Evec1 Projection'); ylabel('Evec2 Projection'), zlabel('Evec3 Projection')
title('Data projection on PCA for evecs 123'); colormap(jet(ncolors))
c=colorbar('Ticks',1:2:ncolors,'TickLabels',vlist(1:2:end)); c.Label.String='Transverse field Jx';
figname=[savedir '/PCA_Projection_evecs'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

% performing k-means on the PCA data (k-means doesn't work on complex vals)
figure;
end
if exist('symm','var')
    if strcmp(symm,'circular')
        if ~no_show
        scatter(y1.^2+y2.^2,y3,s(:),C(:)) % using circular symmetry.
        end
        [idx,centroids]=kmeans([y1.^2+y2.^2 y3],K, 'Replicates',10,'Start','plus');%,'Display','final');
        xlabel('Evec1.^2 +Evec2.^2'); ylabel('Evec3 Projection')
    elseif strcmp(symm,'spherical')
        [idx,centroids]=kmeans([y1.^2+y2.^2+y3.^2],K, 'Replicates',10,'Start','plus');%,'Display','final');
        if ~no_show
        scatter(y1.^2+y2.^2+y3.^2, samples(1:skip:end,:)*U(:,4),s(:),C(:))
        xlabel('Evec1.^2 +Evec2.^2+Evec3.^2'); ylabel('Evec4 Projection'), 
        end
    end
else
    if ~no_show
    scatter(y1, y2,s(:),C(:))
    end
    yk=samples(1:skip:end,:)*U(:,1:kmeans_dim);
    [idx,centroids]=kmeans(yk,K, 'Replicates',10,'Start','plus');%,'Display','final');
end
if ~no_show
xlabel('Evec1'); ylabel('Evec2 Projection'), 
colormap(jet(ncolors))
c=colorbar('Ticks',1:5:ncolors,'TickLabels',vlist(1:5:end));
c.Label.String='Transverse Field B';
figname=[savedir '/PCA_Projection_' symm 'symm'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

figure; gscatter(y1,y2,idx)
title('kmeans clustering of PCA projection')
figname=[savedir '/PCA_Projection_kmeansclustering'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

fig=setfig(12,'f',0,'avg\ cluster\ number',0,'K-means\ phase\ transition\ prediction','off');
plot(vlist,mean(reshape(idx,[dpoints_per_v,ncolors]),1),'-o')
% xlabel('Transverse Field B'); ylabel('avg of cluster number'); 
figname=[savedir '/kmeans_PT_prediction'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')
end

Kpred=mean(reshape(idx,[dpoints_per_v,ncolors]),1);
PCA_pred=(mean(squeeze(reshape((samples*U(:,1)),...
    [sample_size,size(vlist)]))));

end

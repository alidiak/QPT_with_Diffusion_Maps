function [Kpred,Pred,Pred2,autoenc]=Autoenc_kmeans(samples,dpoints_per_v,latent_dim,K, vlist, adjmat, savedir, varargin)

%kmeans_dim=3; % default val
no_show=false; symm='none';
sample_size=max(size(samples))/numel(vlist);
N=min(size(samples)); ncolors=numel(vlist); 
skip=sample_size/dpoints_per_v; % the number of points to be skipped in samples

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'symm'
            symm=varargin{2};
%         case 'kmeans_dim'
%             kmeans_dim=varargin{2};
        case 'apply_transferfunc'
            apply_transferfunc=varargin{2};
        case 'no_plot'
            no_show=varargin{2};
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

autoenc=trainAutoencoder(samples(1:skip:end,:).',latent_dim,'ShowProgressWindow', ~no_show);

if exist('apply_transferfunc','var')
    if strcmp(autoenc.EncoderTransferFunction,'logsig')
        X=logsig(autoenc.EncoderWeights);
    end
else
    X=autoenc.EncoderWeights;
end

C=repmat(1:ncolors,sample_size/skip,1);
s=repmat(5*ones(ncolors,1),sample_size/skip,1);
G = graph(adjmat);

if ~no_show
figure; p=plot(G); G.Nodes.NodeColors =X(1,:).';
p.NodeCData = G.Nodes.NodeColors; colormap(parula(5)); 
colorbar; title('Colormap for largest Autoenc evec')
figname=[savedir '/Graph__1stAutoencEvec'];
savefig(figname); saveas(gcf,[figname '.png'],'png')

if latent_dim>1
    figure; p=plot(G); G.Nodes.NodeColors =X(2,:).';
    p.NodeCData = G.Nodes.NodeColors; colormap(parula(5)); 
    colorbar; title('Colormap for 2nd largest Autoenc evec')
    figname=[savedir '/Graph_2ndAutoencEvec'];
    savefig(figname); saveas(gcf,[figname '.png'],'png')
end

subplot(3,1,1); plot(1:latent_dim,autoenc.EncoderBiases,'-o')
title(['Autoenc Biases']); ylabel('Encoder Biases')
subplot(3,1,2); plot(1:N,X(1,:).','-*'), 
title(['Largest Autoenc Evec'])
xlabel('Lattice Site'), ylabel('Lattice Site Evec Val.')
if latent_dim>1
    subplot(3,1,3); plot(1:N,X(2,:).','-*'), 
    title(['2nd Largest Autoenc Evec'])
    xlabel('Lattice Site'), ylabel('Lattice Site Evec Val.')
end
figname=[savedir '/Autoenc_Results'];
savefig(figname); saveas(gcf,[figname '.png'],'png')
end

y1=samples(1:skip:end,:)*(X(1,:).');
if latent_dim==2
    y2=samples(1:skip:end,:)*(X(2,:).');
    if ~no_show
        figure; scatter(y1,y2,s(:),C(:)); 
    end
elseif latent_dim>=3
     y2=samples(1:skip:end,:)*(X(2,:).'); y3=samples(1:skip:end,:)*(X(3,:).');
    if ~no_show
    figure; scatter3(y1,y2,y3,s(:),C(:)); zlabel('Evec3 Projection')
    end
elseif latent_dim==1 && ~no_show
    figure; plot(vlist,mean(reshape(y1,[dpoints_per_v,ncolors]),1),'.-')
end

if ~no_show
xlabel('Evec1 Projection'); ylabel('Evec2 Projection')
colormap(jet(ncolors));  c=colorbar('Ticks',1:5:ncolors,'TickLabels',vlist(1:5:end));
c.Label.String='Transverse Field B';
figname=[savedir '/Autoenc_LatentSpace_'  num2str(latent_dim) 'Data_Projection'];
savefig(figname); saveas(gcf,[figname '.png'],'png')

figure, plot(vlist,abs(mean(squeeze(reshape((samples*X(1,:).'),...
    [sample_size,size(vlist)])))),'-o','DisplayName',['Evec' num2str(1)])
title('Data projection of evecs as a function of var')
figname=[savedir '/Autoenc_data_Projection'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

end
 

if strcmp(symm,'circular')
    [idx,centroids]=kmeans([y1.^2+y2.^2],K, 'Replicates',10,'Start','plus');%,'Display','final');
    if latent_dim>=3 && ~no_show
        figure; scatter(y1.^2+y2.^2,y3,s(:),C(:)) % using circular symmetry.
        xlabel('Evec1.^2 +Evec2.^2'); ylabel('Evec3 Projection')
    end
elseif strcmp(symm,'spherical') && latent_dim>=3
    [idx,centroids]=kmeans([y1.^2+y2.^2+y3.^2],K, 'Replicates',10,'Start','plus');%,'Display','final');
    % scatter(y1.^2+y2.^2+y3.^2, samples(1:skip:end,:)*U(:,4),s(:),C(:))
else
    if ~no_show
        figure; scatter(y1, y2,s(:),C(:))
    end
    yk=samples*(X(:,:).');
    [idx,centroids]=kmeans(yk,K, 'Replicates',10,'Start','plus');%,'Display','final');
end
if ~no_show
xlabel('Evec1'); ylabel('Evec2 Projection'), colormap(jet(ncolors))
c=colorbar('Ticks',1:5:ncolors,'TickLabels',vlist(1:5:end));
c.Label.String='Transverse Field B';
figname=[savedir '/Autoenc_Projection_' symm 'symm'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')

figure; gscatter(y1,y2,idx); title('kmeans clustering of Autoenc projection');
figname=[savedir '/Autoenc_LatentSpace_'  num2str(latent_dim) '_kmeansclustering'];
savefig(figname); saveas(gcf,[figname '.png'],'png')

fig=setfig(12,'f',0,'avg\ cluster\ number',0,'K-means\ phase\ transition\ prediction','off');
plot(vlist,mean(reshape(idx,[dpoints_per_v,ncolors]),1),'-o')
% xlabel('Transverse Field B'); ylabel('avg of cluster number'); title('K-means phase transition prediction')
figname=[savedir '/Autoenc_LatentSpace_'  num2str(latent_dim) '_kmeans_PT_prediction'];
savefig(figname); saveas(gcf,[figname '.png'],'png')
end 
Kpred=mean(reshape(idx,[dpoints_per_v,ncolors]),1);
Pred=(mean(squeeze(reshape((samples*X(1,:).'),...
    [sample_size,size(vlist)]))));
Pred2=(mean(squeeze(reshape((samples*X(2,:).'),...
    [sample_size,size(vlist)]))));
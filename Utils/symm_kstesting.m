function [ks_tests,ks_stats] = symm_kstesting(samples1,samples2,hlist,S,symmstr,varstr,varargin)
% This function is for comparing two sample sets samples1 and samples2
% using the Kolmogorov-Smirnov test/statistic. Can either input samples1 
% and samples2 as the samples in linear index representation or as the
% samples themselves (then represented by the function getR).

plot_on=true; statistic='KS';
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'statistic'
            statistic=varargin{2};
        case 'plot'
            plot_on=varargin{2};
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

if size(samples1)~=size(samples2)
    error('sample sizes are different')
end

sample_size=max(size(samples1))/length(hlist);

if min(size(samples1))~=1 % means samples were entered
    indices=getR(samples1,S); % converts to a linear index
else
    indices=samples1;
end
if min(size(samples2))~=1
    symm_indices=getR(samples2,S); 
else
    symm_indices=samples2;
end

ks_tests=zeros([max(size(hlist)),1]); ks_hypothesis=zeros([max(size(hlist)),1]); 
ks_stats=zeros([max(size(hlist)),1]); % avg_mag=zeros([max(size(hlist)),1]); 
%kuip_stat=zeros([max(size(hlist)),1]);  cramer_stat=zeros([max(size(hlist)),1]); options={'cramer',true};
if exist('statistic','var')
    options={statistic,true};
else
    options={'ks',true};
end
for kk=1:max(size(hlist))
    samplepos_h=indices(((kk-1)*sample_size+1):(kk*sample_size));
    cmprto=symm_indices(((kk-1)*sample_size+1):(kk*sample_size));
    %kuip_stat(kk)=kuipertest2(samplepos_h,cmprto,500,false);
    %cramer_stat(kk)=kuipertest2(samplepos_h,cmprto,500,false,options{:});
    [ks_hypothesis(kk),ks_tests(kk),~]=kstest2(samplepos_h,cmprto);
    ks_stats(kk)=kuipertest2(samplepos_h,cmprto,sample_size,false,options{:});
    % avg_mag(kk)=mean(mean(Z(samplepos_h,:)));
end

%subplot(5,1,1);plot(hlist,kuip_stat,'-*');  title(['Kuiper Statistic comparing with symmetry: ' symm])
%subplot(5,1,2);plot(hlist,cramer_stat,'-*');  title(['Cramer Statistic comparing with symmetry: ' symm])
% subplot(2,1,1); plot(hlist,ks_tests,'-*'); title(['K-S test p value comparing with symmetry: ' symmstr])
if plot_on
title_str=[statistic '\ statistic\ comparing\ with\ symmetry:\ ' symmstr];
fig=setfig(16,varstr,0,[statistic '\ Statistic'],0,title_str,true);
plot(hlist,ks_stats,'k-*'); 
figname=['./' statistic 'plot'];
savefig([figname '.fig']); saveas(gcf,[figname '.png'],'png')
end
%subplot(3,1,3); plot(hlist,avg_mag,'-*'); title(['Average Magnetization'])

end
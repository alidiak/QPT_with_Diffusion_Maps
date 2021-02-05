function [stat]=kuipertest2(s1,s2,res, plt,varargin)
% This function compares input samples one and two based on the 
% Kuiper Test and Kuiper Statistic. It creates two empirical cdfs 
% based on the two samples and then returns the Kuiper statistic.
% Also has the option to plot the cdfs as well as the D- and D+ values 
% that determined the Kuiper statistic.
%
% Inputs:
% s1 - sample set 1
% s2 - sample set 2
% res - resolution for cdf
% plt - true or false for plotting
% 

kuiper=true; cramer=false; ks=false; %default vals
testname='Kuiper';

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'kuiper'
            kuiper=varargin{2}; testname='Kuiper';
        case 'cramer'
            cramer=varargin{2}; testname='Cramer'; kuiper=false;
        case 'ks'
            ks=varargin{2}; testname='KS'; kuiper=false;
        otherwise
            error(['Invalid Input: ' varargin{1}])        
    end
    varargin(1:2)=[];
end

n=length(s1(:));
m=length(s2(:));

% like ecdf function but with a fixed resolution (easier to compare
% ecdfs this way)
x=linspace(0,max(max(s1),max(s2)),res);
f1=zeros([res,1]);f2=zeros([res,1]);
for nn=1:res
    f1(nn)= sum(s1<=x(nn))/n;% this is the cumulative dist func for s1
    f2(nn)=sum(s2<=x(nn))/m;
end

if kuiper
    [dplus, dpind]=max([0 ; f1-f2]);
    [dminus,dmind]=max([0 ; f2-f1]);

    stat=dplus+dminus;
elseif cramer % the cramer-von-mises test
    % normalizing factor from https://www.microsoft.com/en-us/research/
    % wp-content/uploads/2016/02/williams2008specom.pdf
%     alpha=sqrt( (12*res)/(4*res^2 -1)); 
    alpha=1/(sqrt(res));
    stat=alpha*norm(f2-f1); % order doesn't matter, L2 norm
elseif ks
    stat=max(abs(f2-f1));
end
    
if plt
    figure;
    plot(x,[f1,f2],'.-'); hold on
    str={'sample1 ecdf','sample2 ecdf'}; 
    if kuiper
        plot([x(dpind-1),x(dpind-1)],[f1(dpind-1),f2(dpind-1)],'k-') % plotting vertical lines
        plot([x(dmind-1),x(dmind-1)],[f1(dmind-1),f2(dmind-1)],'k-') % plotting vertical lines
        str={'sample1 ecdf','sample2 ecdf','dplus','dminus'};
    end
    legend(str)
    legend('Location','northwest') 
    title([testname ' Test Results'])
end

end
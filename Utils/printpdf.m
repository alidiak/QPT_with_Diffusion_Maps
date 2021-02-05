function printpdf(fig,fn,Legend,LegLoc)
% Print the figure handle fig to pdf file with name fn
% if no fn specified, using the 'tempfig.pdf'
% Legend and legend location optional, but must be provided both at the same time
% Require ghostscript to be installed so we can use 'esptopdf' to efficiently converse eps to pdf

if nargin>2
    Loc={'northeast','northwest','southwest','southeast','best'}; %Using quardrant index 1,2,3,4,(5)
    legend(Legend,'interpreter','latex','location',Loc{LegLoc});
end
print(fig,'tempfig','-depsc');

PA=getenv('PATH');
LD=getenv('DYLD_LIBRARY_PATH');
setenv('PATH','/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin');
setenv('DYLD_LIBRARY_PATH','/opt/local/lib/:');
system('epstopdf tempfig.eps'); %epstopdf require the above environmental variables
system('rm tempfig.eps');
setenv('DYLD_LIBRARY_PATH',LD); %restore back the original matlab's environmental variables
setenv('PATH',PA);

if nargin==2 || nargin==4
    system(['mv tempfig.pdf ',fn,'.pdf']);
%     system(['qlmanage -p ',fn,'.pdf']);
    open([fn,'.pdf']);

else
    open('tempfig.pdf');
%     system('qlmanage -p tempfig.pdf >& /dev/null');
end


end

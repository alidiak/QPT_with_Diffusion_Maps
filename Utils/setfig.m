function fig=setfig(FontSize,XLabel,XLog,YLabel,YLog,Title, show)
% Universal script for setting up figure for exporting to pdf/eps
    set(0,'DefaultAxesFontSize',FontSize);
    set(0,'DefaultAxesFontName','Times'); 
    
    fig=figure; 
    set(fig,'PaperPosition', [0 0 8 6]);
    set(fig,'PaperSize',[8 6]); 
    set(fig,'Visible',show);  
    box on; axis tight; grid on; hold on;
    if nargin>0
        if ~isempty(XLabel)
            xlabel(['$',XLabel,'$'],'color','blue','interpreter','latex');
        end
        if ~isempty(YLabel)
            ylabel(['$',YLabel,'$'],'color','blue','interpreter','latex');
        end
        if XLog
            set(gca,'xscale','log');
        end
        if YLog
            set(gca,'yscale','log');
        end
        if ~isempty(Title)
            title(['$',Title,'$'],'color','black','interpreter','latex');
        end
    end
end

function map_plot(matrix,v_list,title_str,figname)
% just for plotting the spin clock phase map

fig=setfig(16,'\theta',0,'f',0,title_str,'on');
imagesc(matrix); colormap('jet');  % hot is closest to the inverse blackbody that Johannes mentioned
% use gray for grayscale, and jet is very clear
caxis([min(min(matrix)) max(max(matrix))+1e-8]);
colorbar; xticks([1,(size(matrix,2)/2),size(matrix,2)]);
xticklabels({'0','\pi/6','\pi/3'}); 
% xtk=get(gca, 'XTick'); ytk=get(gca, 'YTick');
%set(gca, 'XTick', [xtk(1),xtk(round((max(size(xtk))+1)/2)),xtk(end)],'XTickLabel', {'0','\pi/6','\pi/3'},...
% names={num2str(v_list(1)), num2str(v_list(round(max(size(v_list))/2))),num2str(v_list(end))};
names={'0', '0.5','1'};
yticks([1,(size(matrix,1)/2),size(matrix,1)]);
yticklabels(names)
set(gca,'Xcolor','blue','Ycolor','blue');
% set(gca, 'YTick', [ytk(1),ytk(round((max(size(ytk))+1)/2)),ytk(end)],'YTickLabel', names,...
%     'Ycolor','blue');
colorbar, % xlabel('theta'), ylabel('f (with J=1-f)') title(title_str)    
saveas(gca,[figname '.png'],'png'), savefig([figname '.fig'])

end
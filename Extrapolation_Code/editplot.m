function editplot
xlabel(gca,get(get(gca,'XLabel'),'String'),'Units','normalized','Position',[0.5 -0.07],'FontName','Times New Roman','FontSize',16)
ylabel(gca,get(get(gca,'YLabel'),'String'),'Units','normalized','Position',[-0.07 0.5],'FontName','Times New Roman','FontSize',16)
legend(gca,'boxoff')
title(gca,get(get(gca,'Title'),'String'),'FontName','Times New Roman','FontSize',11,...
    'Units','normalized','Position',[0.5 1])

set(gcf,'Position', [110   143   781   625],...
       'OuterPosition', [101   134   799   723],...
       'Units', 'pixels','WindowStyle','normal',...
       'PaperPosition', [1.9473    8.0080   17.0882   13.6749],...
       'PaperUnits', 'centimeters','PaperPositionMode','auto',...
       'InvertHardcopy','off','Color','white');

set(gca,'OuterPosition', [-0.0382   -0.0095    1.1218    1.0464],...
'Position', [0.1076    0.1056    0.8694    0.8528],...
'Units','normalized',... %'LineWidth',1, 
'FontName', 'Times New Roman', 'FontSize', 14)
axis(gca,'tight')
end
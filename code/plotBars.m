function hbar = plotBars(nmutPCR)

figure
hbar = bar(1:numel(nmutPCR.uniCase), nmutPCR.PCR, 0.5, 'stack');
colormap([1 0 0; 1 1 0; 0 0 0]) % red yellow and black

set(gca, 'XTick', 1:numel(nmutPCR.uniCase))
xtickangle(90)
ylabel('Number of Somatic Mutations','fontsize',20,'fontweight','bold')


set(gca,'ylim',[0,13200],'xticklabel',nmutPCR.uniCase',...
    'linewidth',2,'tickdir','out','fontsize',8,'fontweight','bold',...
    'xlim',[0.5 numel(nmutPCR.uniCase)+0.5],'ticklength',[0.001 0.002],...
    'ytick',[0:100:400,420,13200],'yticklabel',{'0','100','200','300','400',...
    '400+','13200'})
legend('Untreated Only', 'Common', 'Relapsed Only','Location','north' )

breakyaxis([420 13150]);
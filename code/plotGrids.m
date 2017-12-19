function hgrid = plotGrids(Drivers,saviCase)
% plot grids of PCR frequencies of driver genes in each sample.

[nDriver, nCase] = size(Drivers.gPCR);

figure
hgrid = imagesc(Drivers.gPCR);
set(gca, 'XTick',1:nCase ,'XTickLabel',saviCase.uniCase);
set(gca, 'YTick',1:nDriver, 'YTickLabel',Drivers.GeneID);
set(gca, 'TickDir','out', 'TickLength',[0.001 0.002])
xtickangle(90)

colormap([1 1 1; 1 0 0; 0 0 0; 1 1 0]);
% [none = white; primary = red; recurrent = black; common = yellow]

for i = 1.5:1:(nDriver-0.5)
    line([0.5 nCase+0.5], [i i],'Color','k');
end

for j = 1.5:1:(nCase-0.5)
    line([j j], [0.5 nDriver+0.5],'Color','k');
end
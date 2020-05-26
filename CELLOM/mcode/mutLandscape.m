function hland = mutLandscape(savi, kdstr)
% Build a landscape of mutations for longitudinal data
% Bar + Grid

ng = numel(kdstr.kdlist);
nc = numel(kdstr.unicase);
PCR = zeros(nc, 3);
vafcut = 5;
for i = 1:nc
    savi1 = savi(savi.caseidx == i,:);
    PCR(i,1) = nnz(savi1.Primary_freq >= vafcut & savi1.Recurrent_freq < vafcut);
    PCR(i,2) = nnz(savi1.Primary_freq >= vafcut & savi1.Recurrent_freq >= vafcut);
    PCR(i,3) = nnz(savi1.Primary_freq < vafcut & savi1.Recurrent_freq >= vafcut);
end

hland = figure('position',[0 0 2400 550]);
% Bars
subplot(10,1,1:4)
hold on
hb = bar(1:nc,PCR,0.5,'stacked');
hb(1).FaceColor = [1 0 0];
hb(2).FaceColor = [1 1 0];
hb(3).FaceColor = [0 0 0];

xlim([0.5, nc + 0.5])
xticks(1:nc)
xticklabels({})
ylabel('Somatic Mutations','fontsize',16)
ylim([0 1000])
yticks([0 100 200 300 400 500 1000])
yticklabels({'0','100','200','300','400','500','>1,000'})
set(gca,'tickdir','out','TickLength',[0.003 0.003],'fontsize',16,'box','off','linewidth',1.5)

breakyaxis([500 950]);
hold off

%% Grids

Gmat = zeros(nc,ng);
Gmat(kdstr.Pmat == 1) = 1;
Gmat(kdstr.Rmat == 1) = 3;
Gmat(kdstr.Cmat == 1) = 2;
Gmat = Gmat';
Gmat = flipud(Gmat);

subplot(10,1,5:10)
hold on
imagesc(Gmat)
xticks(1:nc)
xticklabels(kdstr.unicase)
xtickangle(90)
xlim([0.5, nc + 0.5])
ylim([0.5, ng + 0.5])
yticks(1:ng)
yticklabels(fliplr(kdstr.kdlist))
colormap([1 1 1; 1 0 0; 1 1 0; 0 0 0]);

for i = 1.5:1:(ng+0.5)
    line([0.5 nc+0.5], [i i],'Color','k');
end

for j = 1.5:1:(nc+0.5)
    line([j j], [0.5 ng+0.5],'Color','k');
end

set(gca,'tickdir','out','TickLength',[0.003 0.003],'fontsize',16,'box','off','linewidth',1.5)
hold off



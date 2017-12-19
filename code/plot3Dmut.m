function h3d = plot3Dmut(geneTable)

ng = size(geneTable,1);

n01 = zeros(ng,1);
n10 = zeros(ng,1);
n11 = zeros(ng,1);

for i = 1:ng
    n01(i) = nnz(geneTable.gPCR(i,:) == 2); % recurrent
    n10(i) = nnz(geneTable.gPCR(i,:) == 1); % primary
    n11(i) = nnz(geneTable.gPCR(i,:) == 3); % common
end

h3d = figure;
hold on
plot([0 15],[0 -15],'color','k','linewidth',2)
plot([0 -15],[0 -15],'color','r','linewidth',2)
plot([0 0],[0 45],'color','y','linewidth',2)

arrow([0 0],[15 -15],20,'BaseAngle',50,'EdgeColor','k','FaceColor','k')
arrow([0 0],[-15 -15],20,'BaseAngle',60,'EdgeColor','k','FaceColor','r')
arrow([0 0],[0 46],20,'BaseAngle',60,'EdgeColor','k','FaceColor','y')

x = n01-n10;
y = 2*n11-n10-n01;

colors = n01./(n01+n10+n11)*[0 0 0]+n10./(n01+n10+n11)*[1 0 0]+n11./(n01+n10+n11)*[1 1 0];

n00 = 90-n01-n10-n11;

s = 90-n00+10;

for i = 1:ng
    if ~isnan(colors(i,1))
        plot( x(i) , y(i) ,'o', 'MarkerSize' , s(i) , 'MarkerFaceColor' , colors(i,:) , 'MarkerEdgeColor', 'k')
        text( x(i)+1, y(i), geneTable.GeneID{i} ,'fontsize',16,'fontweight','bold','FontAngle','italic')
    end
end


set(gca,'xlim',[-30 30] ,'ylim',[-35 60] ,'box','on','tickdir','out','linewidth',2,'fontsize',16)
axis off
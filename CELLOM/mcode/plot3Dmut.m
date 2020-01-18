function h3d = plot3Dmut(kdstr)
% Plot mutations on 3D space: primary, common and recurrence.

genelist = kdstr.kdlist;
ng = numel(genelist);

GC = kdstr.Cmat';
GP = kdstr.Pmat';
GR = kdstr.Rmat';

n01 = zeros(ng,1);
n10 = zeros(ng,1);
n11 = zeros(ng,1);

for i = 1:ng
    n01(i) = nnz(GR(i,:)); % recurrent
    n10(i) = nnz(GP(i,:)); % primary
    n11(i) = nnz(GC(i,:)); % common
end

%%

h3d = figure('position',[1200 600 900 600]);
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

s = (n01 + n10 + n11)/3 + 10;

rng(1)
for i = 1:ng
    if ~isnan(colors(i,1))
        plot(x(i), y(i), 'o', 'MarkerSize' , s(i) , 'MarkerFaceColor' , colors(i,:) , 'MarkerEdgeColor', 'k')
        if x(i) < 2
            text(x(i)-1, y(i)+1, genelist{i} ,'fontsize',14, 'horizontalalignment','right')
        else
            text(x(i)+1, y(i)+1, genelist{i} ,'fontsize',14, 'horizontalalignment','left')
        end
    end
end


set(gca,'xlim',[-30 30] ,'ylim',[-35 60] ,'box','on','tickdir','out','linewidth',1.5,'fontsize',12)
axis off

hold off
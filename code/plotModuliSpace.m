function hmod = plotModuliSpace(caseTable)
% plot 3-D triangle for evolutionary moduli space

hmod = figure;
hold on

[x,y,z] = sphere(40); % generate coordinates of a sphere: 40-by-40 faces.
Q = (x >= 0 & y >= 0 & z >= 0); % 1st Quadrant
surf(x.*Q,y.*Q,z.*Q, 'FaceColor','white','FaceLighting','none')
view([1.2,1.2,1.2]);
xlim([-0.2 1.2])
ylim([-0.2 1.2])
zlim([-0.2 1.2])
axis off

ns = numel(caseTable.uniCase);

F = caseTable.PCR ./ repmat(sum(caseTable.PCR,2),1,3); % normalization
B = 1./(F(:,1).^2+F(:,2).^2+ F(:,3).^2)*ones(1,3);
S = F.*sqrt(B);

colors = {'r';'y';'k'};
idx = kmeans(S,3,'Replicates',1000);

r = 50*0.0005+0.006;
for i = 1:ns
    [xx,yy,zz] = meshgrid(S(i,1):(r/10):(2*r)+S(i,1),...	% xx = primary
        S(i,3):(r/10):(2*r)+S(i,3),...                      % yy = recurrent
        S(i,2):(r/10):(2*r)+S(i,2));                        % zz = common
    
    rr = r*sqrt((xx-S(i,1)-r).^2 + (yy-S(i,3)-r).^2 + (zz-S(i,2)-r).^2);
    
    p = patch(isosurface(xx,yy,zz,rr));
    
    set(p,'FaceColor',colors{idx(i)},'EdgeColor','none','FaceLighting','gouraud');
end
camlight
axis square

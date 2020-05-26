function hmod = mutTreeClustering(savi)
% Embed evolutionary trees into Moduli Space

[unicase, ~, cix] = unique(savi.CaseID);
numcases = numel(unicase);
PCR = zeros(numcases, 3);
vafcut = 5;
for i = 1:numcases
    savi1 = savi(cix == i,:);
    PCR(i,1) = nnz(savi1.Blood_freq <= 1 & savi1.Primary_freq >= vafcut & savi1.Recurrent_freq < vafcut);
    PCR(i,2) = nnz(savi1.Blood_freq <= 1 & savi1.Primary_freq >= vafcut & savi1.Recurrent_freq >= vafcut);
    PCR(i,3) = nnz(savi1.Blood_freq <= 1 & savi1.Primary_freq < vafcut & savi1.Recurrent_freq >= vafcut);
end

hmod = figure('position',[600 0 600 600]);
hold on

[x,y,z] = sphere(40); % generate coordinates of a sphere: 40-by-40 faces.
Q = (x >= 0 & y >= 0 & z >= 0); % 1st Quadrant
surf(x.*Q,y.*Q,z.*Q, 'FaceColor','white','FaceLighting','none')
view([1.2,1.2,1.2]);
xlim([-0.2 1.2])
ylim([-0.2 1.2])
zlim([-0.2 1.2])
axis off

%%
F = PCR./repmat(sum(PCR,2), 1, 3); % normalization
B = 1./(F(:,1).^2+F(:,2).^2+ F(:,3).^2)*ones(1,3);
S = F.*sqrt(B);

rng(1)
idx = kmeans(S,3,'Replicates',5);
colors = cell(3,1);
[~,i1] = max(S(:,1));
colors{idx(i1)} = 'r';
[~,i2] = max(S(:,2));
colors{idx(i2)} = 'y';
[~,i3] = max(S(:,3));
colors{idx(i3)} = 'k';

r = 50*0.0005+0.006;
for i = 1:numcases
    [xx,yy,zz] = meshgrid(S(i,1):(r/10):(2*r)+S(i,1),...	% xx = primary
        S(i,3):(r/10):(2*r)+S(i,3),...                      % yy = recurrent
        S(i,2):(r/10):(2*r)+S(i,2));                        % zz = common
    
    rr = r*sqrt((xx-S(i,1)-r).^2 + (yy-S(i,3)-r).^2 + (zz-S(i,2)-r).^2);
    
    p = patch(isosurface(xx,yy,zz,rr));
    
    set(p,'FaceColor',colors{idx(i)},'EdgeColor','none','FaceLighting','gouraud');
end
camlight
axis square
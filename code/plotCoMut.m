function [hcom,sigNum] = plotCoMut(geneTable,doPlot)
% plot mut-mut relationships: co-occurrence vs. mutually exclusive in
% primary and recurrent tumor samples, respectively.

ng = size(geneTable,1);

M = geneTable.gPCR;
MP = zeros(size(M));
MR = zeros(size(M));
MP(M == 1 | M == 3) = 1;
MR(M == 2 | M == 3) = 1;

dotx = zeros(ng*ng - ng, 1);
doty = zeros(ng*ng - ng, 1);
sz = zeros(ng*ng - ng, 1);
cr = zeros(ng*ng - ng, 3);
% cr = color
% red = [1 0 0]
% grey = [0.75 0.75, 0.75]
% black = [0 0 0]
% green = [0 1 0]

sigNum = zeros(ng,1);

t = 0;
for i = 1:ng
    for j = 1:ng
        Z = zeros(2,2);
        if i < j % primary
            t = t + 1;
            dotx(t) = i;
            doty(t) = j;
            Z(1,1) = nnz(MP(i,:) == 1 & MP(j,:) == 1);
            Z(1,2) = nnz(MP(i,:) == 1 & MP(j,:) == 0);
            Z(2,1) = nnz(MP(i,:) == 0 & MP(j,:) == 1);
            Z(2,2) = nnz(MP(i,:) == 0 & MP(j,:) == 0);
            Z2 = Z + ones(2,2);
            oz = Z2(1,1)*Z2(2,2)/(Z2(1,2)*Z2(2,1));
            [~,p,~] = fishertest(Z);
            if p < 0.1
                sigNum(i) = sigNum(i) + 1;
                sz(t) = -log10(p);
                if oz > 1
                    cr(t,:) = [1 0 0]; % red
                else
                    cr(t,:) = [0 1 0]; % green
                end
            else
                sz(t) = 1;
                cr(t,:) = [0.75 0.75 0.75]; % grey
            end

        elseif i > j % relapsed
            t = t + 1;
            dotx(t) = i;
            doty(t) = j;
            Z(1,1) = nnz(MR(i,:) == 1 & MR(j,:) == 1);
            Z(1,2) = nnz(MR(i,:) == 1 & MR(j,:) == 0);
            Z(2,1) = nnz(MR(i,:) == 0 & MR(j,:) == 1);
            Z(2,2) = nnz(MR(i,:) == 0 & MR(j,:) == 0);
            Z2 = Z + ones(2,2);
            oz = Z2(1,1)*Z2(2,2)/(Z2(1,2)*Z2(2,1));
            [~,p,~] = fishertest(Z);
            if p < 0.1
                sigNum(i) = sigNum(i) + 1;
                sz(t) = -log10(p);
                if oz > 1
                    cr(t,:) = [0 0 0]; % black
                else
                    cr(t,:) = [0 1 0]; % green
                end
            else
                sz(t) = 1;
                cr(t,:) = [0.75 0.75 0.75]; % grey
            end
        end
    end
end

if strcmp(doPlot, 'plot')
    figure
    hcom = scatter(dotx, doty, sz*100, cr, 'filled');
    hold on
    axis square
    xticks(1:ng)
    xticklabels(geneTable.GeneID)
    xtickangle(90)
    yticks(1:ng)
    yticklabels(geneTable.GeneID)
    xlim([0,ng+1])
    ylim([0,ng+1])
    set(gca, 'TickDir','out')
    box on
    hold off
else
    hcom = 0;
end

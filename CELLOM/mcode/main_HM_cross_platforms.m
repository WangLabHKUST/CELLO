% CELLO extension of hypermutation detection to WES, targeted and RNA-seq

close all
clear
clc

% Reading ...
W = readtable('input.wes.savi.txt');
T = readtable('input.targeted.savi.txt');
R = readtable('input.rna.savi.txt');

% Calculating ...
HW = calcHMmultiseq(W);
HT = calcHMmultiseq(T);
HR = calcHMmultiseq(R);

% Plotting ...
hmcut = 1.3;
fz = 24;
fa = 0.5;

figure('position',[0 0 1800 600])
subplot(1,3,1)
hold on
numcases = size(HW,1);
tmbcut = 350;
scatter(HW.mutload, HW.tmz, ones(numcases,1)*100, repmat(hex2rgb('#7F152E'),numcases,1),'o', 'LineWidth',1.5)
numhm = nnz(HW.tmz > hmcut);
disp(numhm)
scatter(HW.mutload(HW.tmz > hmcut), HW.tmz(HW.tmz > hmcut), ...
    ones(numhm,1)*100, repmat(hex2rgb('#7F152E'),...
    numhm,1),'filled', 'LineWidth',1.5,'MarkerFaceAlpha',fa)
line([1 10^4],[hmcut hmcut],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
line([tmbcut tmbcut],[eps 2],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
xlabel('Tumor Mutation Load')
ylabel('HM Score')
yticks([0 1 1.3 2])
yticklabels({'0','1','1.3','2'})
axis square
title('Whole-Exome Sequencing','FontWeight','normal')
set(gca,'tickdir','out','TickLength',[0.01 0.01],'fontsize',fz,'box','off','linewidth',1.5,'XScale','log')
hold off

%==========================================================================

subplot(1,3,2)
hold on
numcases = size(HT,1);
tmbcut = 60;
scatter(HT.mutload, HT.tmz, ones(numcases,1)*100, repmat(hex2rgb('#DB9501'),numcases,1),'o', 'LineWidth',1.5)
numhm = nnz(HT.tmz > hmcut);
disp(numhm)
scatter(HT.mutload(HT.tmz > hmcut), HT.tmz(HT.tmz > hmcut), ...
    ones(numhm,1)*100, repmat(hex2rgb('#DB9501'),...
    numhm,1),'filled', 'LineWidth',1.5,'MarkerFaceAlpha',fa)
line([1 10^4],[hmcut hmcut],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
line([tmbcut tmbcut],[eps 2],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
xlabel('Tumor Mutation Load')
ylabel('HM Score')
yticks([0 1 1.3 2])
yticklabels({'0','1','1.3','2'})
axis square
title('Targeted-DNA Sequencing','FontWeight','normal')
set(gca,'tickdir','out','TickLength',[0.01 0.01],'fontsize',fz,'box','off','linewidth',1.5,'XScale','log')
hold off

%==========================================================================

subplot(1,3,3)
hold on
numcases = size(HR,1);
tmbcut = 60;
scatter(HR.mutload, HR.tmz, ones(numcases,1)*100, repmat(hex2rgb('#258039'),numcases,1),'o', 'LineWidth',1.5)
numhm = nnz(HR.tmz > hmcut);
disp(numhm)
scatter(HR.mutload(HR.tmz > hmcut), HR.tmz(HR.tmz > hmcut), ...
    ones(numhm,1)*100, repmat(hex2rgb('#258039'),...
    numhm,1),'filled', 'LineWidth',1.5,'MarkerFaceAlpha',fa)
line([1 10^4],[hmcut hmcut],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
line([tmbcut tmbcut],[eps 2],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
xlabel('Tumor Mutation Load')
ylabel('HM Score')
yticks([0 1 1.3 2])
yticklabels({'0','1','1.3','2'})
axis square
title('RNA Sequencing','FontWeight','normal')
set(gca,'tickdir','out','TickLength',[0.01 0.01],'fontsize',fz,'box','off','linewidth',1.5,'XScale','log')
hold off
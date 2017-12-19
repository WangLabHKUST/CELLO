function [predDrivers,predStats,missGene,AUC,F1] = Ensemble(saviGene,doPlot)
% ensemble = RF + SVM

candidateGene = saviGene(saviGene.nNonHMsamples >= 3,:);
% not necessary feature:
%[~,sigmut] = plotCoMut(candidateGene,'noPlot');
%candidateGene.sigmut = sigmut;

ng = size(candidateGene,1);
ntest = 100;
ntrees = 9;

F1 = zeros(ntest*2,1);
AUC = zeros(ntest*2,1);
Ytrue = candidateGene.isKnownDriver;
Yp = zeros(ng, ntest*2);
c = 0;

for i = 1:ntest
    c = c + 1;
    [Yl,Ys] = computeDriver(candidateGene,'RandomForest',ntrees,'CV','noWeightPlot');
    [AUC(c),~,~] = calcauroc(Ytrue,Ys,'noplot','Random Forest');
    Yp(:,c) = Yl;
    f1 = 2*nnz(Yl + Ytrue == 2)/(2*nnz(Yl + Ytrue == 2) + nnz(Yl + Ytrue == 1));
    F1(c) = f1;
    disp(['RF test ',num2str(i),' out of ',num2str(ntest),': AUC = ',num2str(AUC(c)),' and F1 = ',num2str(F1(c))])
end

for i = 1:ntest
    c = c + 1;
    [Yl,Ys] = computeDriver(candidateGene,'SVM','linear','CV','noWeightPlot');
    [AUC(c),~,~] = calcauroc(Ytrue,Ys,'noplot','SVM');
    Yp(:,c) = Yl;
    f1 = 2*nnz(Yl + Ytrue == 2)/(2*nnz(Yl + Ytrue == 2) + nnz(Yl + Ytrue == 1));
    F1(c) = f1;
	disp(['SVM test ',num2str(i),' out of ',num2str(ntest),': AUC = ',num2str(AUC(c)),' and F1 = ',num2str(F1(c))])
end

hitscutoff = 2*ntest/100;
allhits = sum(Yp,2);
[sorthits, iy] = sort(allhits, 'descend');
selectedgeneix = iy(1:nnz(sorthits >= hitscutoff));
predDrivers = candidateGene(selectedgeneix,:);

GeneName = candidateGene.GeneID(selectedgeneix);
numHits = allhits(selectedgeneix);
knownDriverTag = Ytrue(selectedgeneix);
for i = 1:numel(selectedgeneix)
    if knownDriverTag(i) == 0
        GeneName{i} = [GeneName{i},' (new)'];
    end
end
rfHits = sum(Yp(:,1:ntest),2);
rfHits = rfHits(selectedgeneix);
svmHits = sum(Yp(:,(ntest+1):end),2);
svmHits = svmHits(selectedgeneix);
predStats = table(GeneName,numHits,knownDriverTag,rfHits,svmHits);
missGene = candidateGene(allhits < hitscutoff & Ytrue == 1,:);

if strcmp(doPlot,'plot')
    figure;
    bar([rfHits, svmHits],'stacked')
    legend({'RF Hits','SVM Hits'})
    npredgene = size(predStats,1);
    xlim([0 npredgene + 1])
    xticks(1:npredgene)
    xticklabels(predStats.GeneName)
    xtickangle(45)
    ylabel('Number of Hits')
end

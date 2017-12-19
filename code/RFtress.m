% test different number of trees used in RF
candidateGene = saviGene(saviGene.nNonHMsamples >= 3,:);

%[~,sigmut] = plotCoMut(candidateGene,'noPlot');
%candidateGene.sigmut = sigmut;

ntest = 3;
ntrees = [5,10,50,100,500,1000];
F = zeros(ntest,numel(ntrees));

Ytrue = candidateGene.isKnownDriver;
Yp = zeros(size(candidateGene,1), numel(ntrees)*ntest);
c = 0;

for t = 1:numel(ntrees)
    for i = 1:ntest
        c = c + 1;
        [Ypred,Yscore] = computeDriver(candidateGene,'RandomForest',ntrees(t),'CV','noWeightPlot');
        [auc,fpr,tpr] = calcauroc(Ytrue,Yscore,'noplot','Random Forest');
        Yp(:,c) = Ypred;
        F1 = 2*nnz(Ypred + Ytrue == 2)/(2*nnz(Ypred + Ytrue == 2) + nnz(Ypred + Ytrue == 1));
        disp(['Cross validation F1 = ',num2str(F1)]);
        F(i,t) = F1;
    end
end

figure;
boxplot(F);
set(gca, 'TickDir','out')
xticklabels(ntrees)
xlabel('Number of Trees')
ylabel('F1 Score')

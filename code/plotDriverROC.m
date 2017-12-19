function hroc = plotDriverROC(saviGene)

candidateGene = saviGene(saviGene.nNonHMsamples >= 3,:);
disp(['Number of known drivers = ',num2str(nnz(candidateGene.isKnownDriver))])

% test empirical
Ytrue = candidateGene.isKnownDriver;
[Ylabel,~] = computeDriver(candidateGene,'EmpiricalFilter','none','noCV','noWeightPlot');
tp = nnz(Ylabel == 1 & Ytrue == 1);
fp = nnz(Ylabel == 1 & Ytrue == 0);
fn = nnz(Ylabel == 0 & Ytrue == 1);
tn = nnz(Ylabel == 0 & Ytrue == 0);
tpr0 = tp/(tp + fn);
fpr0 = fp/(fp + tn);

% test RF
ntrees = [9,49,99];
AUCrf = zeros(numel(ntrees),1);
hroc = figure;
hold on
for i = 1:numel(ntrees)
    [~,Yscore] = computeDriver(candidateGene,'RandomForest',ntrees(i),'CV','noWeightPlot');
    [AUCrf(i),~,~] = calcauroc(candidateGene.isKnownDriver,Yscore,'plot1','Random Forest');
end

% test SVM
kernels = {'linear', 'gaussian', 'polynomial'};
AUCsvm = zeros(numel(kernels),1);
for i = 1:numel(kernels)
    [~,Yscore] = computeDriver(candidateGene,'SVM',kernels{i},'CV','noWeightPlot');
    [AUCsvm(i),~,~] = calcauroc(candidateGene.isKnownDriver,Yscore,'plot1',['SVM-',kernels{i}]);
end

plot(fpr0, tpr0, 'ko')
plot([0 1], [0 1], 'k--', 'LineWidth', 2)

legend({['RF: 9 trees, AUC = ',num2str(AUCrf(1),2)],['RF: 49 trees, AUC = ',num2str(AUCrf(2),2)],['RF: 99 trees, AUC = ',num2str(AUCrf(3),2)],...
    ['SVM: linear, AUC = ',num2str(AUCsvm(1),2)],['SVM: gaussian, AUC = ',num2str(AUCsvm(2),2)],['SVM: polynomial, AUC = ',num2str(AUCsvm(3),2)],...
    'Empirical Filtering', 'Random, AUC = 0.5'},'Location','southeast')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
axis square
box on
hold off

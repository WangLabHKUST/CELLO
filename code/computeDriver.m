function [Ylabel,Yscore] = computeDriver(geneTable,Model,param,doCV,plotWeights)
% Identify driver genes using direct selection or random forest.

ng = size(geneTable,1);
%nT = 20; % number of trees
%doCV = 'CV'; % whether do cross validation or not.
%plotWeights = 'WeightPlot';
cvFolds = 5; % ?-fold cross validation

if strcmp(Model,'EmpiricalFilter')
    driverFilter = (geneTable.nNonHMsamples >= 5) & ...	% var must occur in at least 5 non-HM sample
        (geneTable.nLowMut <= 3) & (geneTable.maxSgt1maxFreq >= 20) & ...	% Low-impact vars are less than 3; max freq in tumors is greater than 20
        ((geneTable.isKnownDriver == 1) | (geneTable.pvPRpropDiff < 0.05) | (geneTable.isSwitch > 0));
    Ylabel = zeros(ng,1);
    Ylabel(driverFilter) = 1;
    Yscore = Ylabel; % dummy
    
elseif strcmp(Model,'RandomForest')
    Y = geneTable.isKnownDriver; % column vector
	tmpTable = geneTable;
	tmpTable(:,'isKnownDriver') = [];
    tmpTable(:,'GeneID') = [];
    tmpTable(:,'SwitchCases') = [];
    tmpTable(:,'gPCR') = [];
    X = tmpTable{:,:};
    
    % B = TreeBagger(NumTrees,X,Y,'OOBPredictorImportance','on')
    if strcmp(doCV,'CV')
        [CX,CY,CI] = splitGenes(X,Y,cvFolds);
        Ylabel = zeros(ng,1);
        Yscore = zeros(ng,1);
        
        for i = 1:cvFolds
            tmpX = CX;
            testX = tmpX{i};
            tmpX{i} = [];
            trainX = cell2mat(tmpX);
            
            tmpY = CY;
            tmpY{i} = [];
            trainY = cell2mat(tmpY);
            
            testI = CI{i};
            
            B = TreeBagger(param,trainX,trainY,'Method','classification');
            [yp,ys] = predict(B,testX);
            Ylabel(testI) = cell2num(yp);
            Yscore(testI) = ys(:,2); % ys(:,2) = Prob(label = 1)
        end
        
    elseif strcmp(doCV,'noCV')
        ixpos = find(Y == 1);
        npos = numel(ixpos);
        mutepos = ceil(npos/cvFolds);
        ixmute = randperm(npos,mutepos);
        Y(ixpos(ixmute)) = 0;
        B = TreeBagger(param,X,Y,'Method','classification','OOBPredictorImportance','on');
        if strcmp(plotWeights,'WeightPlot')
            figure;
            nf = numel(B.OOBPermutedVarDeltaError);
            bar(B.OOBPermutedVarDeltaError)
            xlim([0 nf+1])
            xticklabels(tmpTable.Properties.VariableNames)
            xtickangle(45)
            ylabel('Predictor Importance')
        end
        [Yl,Ys] = predict(B,X);
        Yscore = Ys(:,2);
        Ylabel = cell2num(Yl);
    else
        error('unknown CV option!')
    end
    
elseif strcmp(Model,'SVM')
    
    Y = geneTable.isKnownDriver; % column vector
	tmpTable = geneTable;
	tmpTable(:,'isKnownDriver') = [];
    tmpTable(:,'GeneID') = [];
    tmpTable(:,'SwitchCases') = [];
    tmpTable(:,'gPCR') = [];
    X = tmpTable{:,:};
    
    % Mdl = fitcsvm(Tbl,Y)
    if strcmp(doCV,'CV')
        [CX,CY,CI] = splitGenes(X,Y,cvFolds);
        Ylabel = zeros(ng,1);
        Yscore = zeros(ng,1);
        
        for i = 1:cvFolds
            tmpX = CX;
            testX = tmpX{i};
            tmpX{i} = [];
            trainX = cell2mat(tmpX);
            
            tmpY = CY;
            tmpY{i} = [];
            trainY = cell2mat(tmpY);
            
            testI = CI{i};
            
            SVM = fitcsvm(trainX,trainY,'Standardize',true,'KernelFunction',param);
            [yp,ys] = predict(SVM,testX);
            Ylabel(testI) = yp;
            Yscore(testI) = ys(:,2); % ys(:,2) = Prob(label = 1)
        end
        
    elseif strcmp(doCV,'noCV')
        ixpos = find(Y == 1);
        npos = numel(ixpos);
        mutepos = ceil(npos/cvFolds);
        ixmute = randperm(npos,mutepos);
        Y(ixpos(ixmute)) = 0;
        SVM = fitcsvm(X,Y,'Standardize',true,'KernelFunction',param);
        [Ylabel,Ys] = predict(SVM,X);
        Yscore = Ys(:,2);
    else
        error('unknown CV option!')
    end
    
    
else
    error('unknown model!')
end

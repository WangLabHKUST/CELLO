function [HMplot, savi] = getHyperMut(savi)
% Identify hyper-mutated samples

varPresentcut = 15;
cutMutLoad = 400;
cutHMscore = 1.0;

nvar = size(savi,1);
uniPatient = unique(savi.CaseID);
npatient = numel(uniPatient);

label = [uniPatient; uniPatient]; % for primary and recurrent, respectively
mutload = zeros(npatient*2,1);
hmscore = zeros(npatient*2,1);

for i = 1:npatient
    % for each patient
    saviPatient = savi(strcmp(savi.CaseID, uniPatient{i}),:);
    
    % for primary
    label{i} = [label{i},'P'];
    
    Pvar = saviPatient.Blood_freq == 0 & saviPatient.Primary_freq > varPresentcut;
    nPvar = nnz(Pvar);
    mutload(i) = nPvar;
    saviPatientP = saviPatient(Pvar,:);
    
    nPctga = nnz((strcmp(saviPatientP.ref, 'C') & strcmp(saviPatientP.alt, 'T')) | ...
        (strcmp(saviPatientP.ref, 'G') & strcmp(saviPatientP.alt, 'A')));
    
    nPccgg = nnz((strcmp(saviPatientP.ref, 'C') & strcmp(saviPatientP.varSurffix, 'C')) | ...
        (strcmp(saviPatientP.ref, 'G') & strcmp(saviPatientP.varPrefix, 'G')));
    
    nPcggc = nnz((strcmp(saviPatientP.ref, 'C') & strcmp(saviPatientP.varSurffix, 'G')) | ...
        (strcmp(saviPatientP.ref, 'G') & strcmp(saviPatientP.varPrefix, 'C')));
    
    hmscore(i) = (nPccgg - nPcggc)/(nPctga + 1) + nPctga/(nPvar + 1);
    
    % for recurrent
    ii = i + npatient;
    label{ii} = [label{ii},'R'];
    
	Rvar = saviPatient.Blood_freq == 0 & saviPatient.Recurrent_freq > varPresentcut;
    nRvar = nnz(Rvar);
    mutload(ii) = nRvar;
    saviPatientR = saviPatient(Rvar,:);
    
    nRctga = nnz((strcmp(saviPatientR.ref, 'C') & strcmp(saviPatientR.alt, 'T')) | ...
        (strcmp(saviPatientR.ref, 'G') & strcmp(saviPatientR.alt, 'A')));
    
    nRccgg = nnz((strcmp(saviPatientR.ref, 'C') & strcmp(saviPatientR.varSurffix, 'C')) | ...
        (strcmp(saviPatientR.ref, 'G') & strcmp(saviPatientR.varPrefix, 'G')));
    
    nRcggc = nnz((strcmp(saviPatientR.ref, 'C') & strcmp(saviPatientR.varSurffix, 'G')) | ...
        (strcmp(saviPatientR.ref, 'G') & strcmp(saviPatientR.varPrefix, 'C')));
    
    hmscore(ii) = (nRccgg - nRcggc)/(nRctga + 1) + nRctga/(nRvar + 1);
    
end


savi.isHM = zeros(nvar,1);
hmmark = zeros(npatient*2,1);

for i = 1:npatient
    if (mutload(i) > cutMutLoad && hmscore(i) > cutHMscore)
        savi.isHM(strcmp(savi.CaseID, uniPatient{i})) = 1;
        hmmark(i) = 1;
    end
    
    ii = i + npatient;
	if mutload(ii) > cutMutLoad && hmscore(ii) > cutHMscore
        savi.isHM(strcmp(savi.CaseID, uniPatient{i})) = 1;
        hmmark(ii) = 1;
	end
end

% Stacked Bars for Mutation Types
MutTypes = zeros(3,6); % 3 bars with each stacked 6 mut types.
MutTypes(1,:) = countMutTypes(savi(((savi.Primary_freq > varPresentcut) & (savi.Recurrent_freq <= varPresentcut)),:));
MutTypes(2,:) = countMutTypes(savi(((savi.Primary_freq <= varPresentcut) & (savi.Recurrent_freq > varPresentcut) & (savi.isHM == 0)),:));
MutTypes(3,:) = countMutTypes(savi(((savi.Primary_freq <= varPresentcut) & (savi.Recurrent_freq > varPresentcut) & (savi.isHM == 1)),:));
MutTypes = MutTypes ./ (sum(MutTypes,2)*ones(1,6));
figure;
bar(MutTypes, 0.4, 'stacked') % width = 0.4
hold on
xticks([1,2,3])
xticklabels({'Primary','non-HM Relapse','HM Relapse'})
ylabel('Mutation Fraction (%)')
yticks(0:0.1:1)
yticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
set(gca, 'tickdir','out')
legend({'A>C / T>G','A>G / T>C','A>T / T>A','C>A / G>T','C>G / G>C','C>T / G>A'}, 'Location','eastoutside')
hold off

% Mutation Load vs. HM Score
figure;
scatter(mutload(1:npatient), hmscore(1:npatient), [], [1 0 0], 'filled')
hold on
scatter(mutload(npatient+1:end), hmscore(npatient+1:end), [], [0 0 0], 'filled')
scatter(mutload(hmmark == 1), hmscore(hmmark == 1), ones(sum(hmmark),1)*100, [0.75 0.75 0.75],'LineWidth',3)
legend('Primary','Relapsed','Hypermutated','Location','east');

set(gca, 'XScale', 'log')
xlabel('Mutation Load')
ylabel('HM Score')
box on
xlim([1 20000])
ylim([-0.6 1.6])
line([1 2e4],[cutHMscore cutHMscore],'Color','blue','LineStyle','--')
line([cutMutLoad cutMutLoad],[-0.6 1.6],'Color','blue','LineStyle','--')

hold off

HMplot = table(label, mutload, hmscore, hmmark);
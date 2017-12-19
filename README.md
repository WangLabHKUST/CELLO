# CELLO - Cancer EvoLution for LOngitudinal data

## Introduction
Cancer EvoLution for LOngitudinal data (CELLO) is a MATLAB toolbox for comprehensive analysis of longitudinal genomic sequencing data in cancer. It was originally developed by Jiguang Wang [1,2], and was then packed up by [Biaobin Jiang](https://github.com/bbjiang). The code was written in **MATLAB/R2016b**.

## Datasets

The input SAVI report (input.savi.report.txt) consists of a list of genetic variants from 90 glioblastoma patients.

## 1-Click to Run CELLO pipeline: play_cello.m
```matlab
% Main function
% CELLO: Cancer EvoLution for LOngitudinal data
close all
addpath('./code/')
% delete outMutContours directory to avoid warnings ...
tic

% Read savi report, the required 16 columns are:
% 'chr','pos','ref','alt','id_cosmic','Effect_Impact','Amino_Acid_Change','Amino_Acid_length','Gene_Name',
% 'Sgt1_max_frequency','Blood_freq','Primary_freq','Recurrent_freq','varPrefix','varSurfix','CaseID'

savi = readtable('input.savi.report.txt');

kdlist = {'ARID1A','ARID2','ATM','ATR','ATRX','BCOR','CDKN2A','CIC',...
    'DICER1','EGFR','FAT1','FAT2','FAT3','FBXW7','FUBP1','IDH1','LTBP4',...
    'MAX','MED12','MET','MTOR','NF1','NF2','NOTCH1','NOTCH2','PDGFRA',...
    'PIK3CA','PIK3CG','PIK3R1','PTEN','PTPN11','RB1','SMARCA4','TCF12','TP53'};
savi = markKnownDriver(kdlist,savi,'fast');

% Alternatively:
%load inputSavi.mat

filter = (savi.altdepth_Blood == 0 | (savi.altdepth_Blood == 1 & savi.refdepth_Blood >= 25)) & ...
    (savi.Sgt1_max_frequency >= 5);
savi = savi(filter,:);

[HMplot, savi] = getHyperMut(savi);

saviCase = getCaseInfo(savi);
[saviMut,savi] = getMutInfo(savi);
saviGene = getGeneInfo(savi);

nmutPCR = plotMutContours(savi.CaseID, savi.Primary_freq, savi.Recurrent_freq,5,'plot');
hbar = plotBars(nmutPCR);
hmod = plotModuliSpace(nmutPCR);

[predDrivers,predStats,missGene] = Ensemble(saviGene,'plot');
hroc = plotDriverROC(saviGene);

hgrid = plotGrids(predDrivers,saviCase);
hsw = plotSwitch(predDrivers,savi);
[hcom,~] = plotCoMut(predDrivers,'plot');
h3d = plot3Dmut(predDrivers);

evoGenes = {'EGFR','IDH1','RB1','MSH6','LRP1B','PIK3CA','PIK3R1','PIK3CG',...
    'PTEN','NF1','TP53','ATRX','LTBP4','PDGFRA'};
genTEDG(evoGenes, saviGene, saviCase.uniCase);

writetable(saviCase, 'saviCase.txt','Delimiter','\t')
writetable(saviGene, 'saviGene.txt','Delimiter','\t')
writetable(saviMut, 'saviMut.txt','Delimiter','\t')
writetable(predDrivers, 'predDrivers.txt','Delimiter','\t')

disp(['Total elapsed time: ',num2str(toc)])
```

## Reference

[1] Wang, J., Cazzato, E., Ladewig, E., Frattini, V., Rosenbloom, D. I., Zairis, S., ... & Lee, J. K. (2016). Clonal evolution of glioblastoma under therapy. **Nature Genetics**, 48(7), 768-776.

[2] Wang, J., Khiabanian, H., Rossi, D., Fabbri, G., Gattei, V., Forconi, F., ... & Pasqualucci, L. (2014). Tumor evolutionary directed graphs and the history of chronic lymphocytic leukemia. **Elife**, 3, e02869.

## Citation

The abstract of this project has been submitted to [AsianEvo 2018](http://www.asianevo.org/), and is now under review.

## Contact

For technical questions, please contact Biaobin via email: biaobinjiang@gmail.com

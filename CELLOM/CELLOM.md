## 1-Click to Run CELLO pipeline: play_cello.m
```matlab
% Main Function of CELLO.M
% CELLO.M: Cancer EvoLution for LOngitudinal data, a Matlab toolbox.

close all
clear
clc

tic

savi = readtable('input.savi.txt');

% Alternatively:
% load inputSavi.mat

%% Preprocessing: Filtering somatic mutations

somaticfilter = (savi.altdepth_Blood == 0 | (savi.altdepth_Blood == 1 & savi.refdepth_Blood >= 25)) & ...
    (savi.Sgt1_max_frequency >= 5);
savi = savi(somaticfilter,:);

%% Marking key driver genes ...

kdlist = {'TP53','ATRX','IDH1','EGFR','PTEN','PIK3CA','PIK3R1','PIK3CG','PDGFRA','RB1','NF1','PTPN11','LTBP4'};
[savi,kdstr] = markKnownDriver(kdlist,savi);

%% Mutational Landscape of Longitudinal Data

hland = getLandscape(savi, kdstr);

%% co-occurrence of mutations

hcom = corrMut(kdstr);

%% 3D mutation plot

h3d = plot3Dmut(kdstr);

%% Mutational Signature Analysis

hsig = getMutsigs(savi);

%% Moduli Space

hmod = getModuli(savi);

%% TEDG with deconv?

G = getTEDG(savi, true);

%% Clonal Switching

hsw = plotSwitch(savi,'PDGFRA');

%% Output

% save('outCELLO.mat')

disp(['Total elapsed time: ',num2str(toc)])

```

## Reference

[1] Wang, J., Cazzato, E., Ladewig, E., Frattini, V., Rosenbloom, D. I., Zairis, S., ... & Lee, J. K. (2016). Clonal evolution of glioblastoma under therapy. **Nature Genetics**, 48(7), 768-776.

[2] Wang, J., Khiabanian, H., Rossi, D., Fabbri, G., Gattei, V., Forconi, F., ... & Pasqualucci, L. (2014). Tumor evolutionary directed graphs and the history of chronic lymphocytic leukemia. **Elife**, 3, e02869.

## Citation

The abstract of this project has been accepted by [AsiaEvo 2018](http://www.asianevo.org/) for oral presentation.

## Contact

For technical questions, please contact Biaobin via email: biaobinjiang@gmail.com

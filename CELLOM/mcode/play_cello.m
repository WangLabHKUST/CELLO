% CELLO.M: Cancer EvoLution for LOngitudinal data, a Matlab toolbox.
% Main Function: play_cello.m

close all
clear
clc

tic

savi = readtable('../../input.savi.txt');

%% Preprocessing: Filtering somatic mutations

somaticfilter = savi.refdepth_Blood >= 20 & savi.altdepth_Blood <= 1 & savi.Sgt1_max_frequency >= 5;
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

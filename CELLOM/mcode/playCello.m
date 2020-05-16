% CELLO.M: Cancer EvoLution for LOngitudinal data, a Matlab toolbox.
% Main Function: playCello.m

close all
clear
clc

tic

saviTable = readtable('../../input.savi.txt');

%% Preprocessing: Filtering somatic mutations

somaticfilter = saviTable.refdepth_Blood >= 20 & saviTable.altdepth_Blood <= 1 & saviTable.Sgt1_max_frequency >= 5;
saviTable = saviTable(somaticfilter,:);

%% Marking key driver genes ...

knownDrivers = {'TP53','ATRX','IDH1','EGFR','PTEN','PIK3CA','PIK3R1','PIK3CG','PDGFRA','RB1','NF1','PTPN11','LTBP4'};
[saviTable,mutGeneTable] = mutStats(knownDrivers,saviTable);

%% Mutational Landscape of Longitudinal Data

hland = mutLandscape(saviTable, mutGeneTable);

%% co-occurrence of mutations

hcom = mutCorrelation(mutGeneTable);

%% 3D mutation plot

h3d = mutFrequency(mutGeneTable);

%% Mutational Signature Analysis

hsig = mutSignature(saviTable);

%% Moduli Space

hmod = mutTreeClustering(saviTable);

%% TEDG with deconv?

G = mutDirectedGraph(saviTable, true);

%% Clonal Switching

hsw = mutSwitch(saviTable,'PDGFRA');

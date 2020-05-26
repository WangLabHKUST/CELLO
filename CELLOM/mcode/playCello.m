% CELLO.M: Cancer EvoLution for LOngitudinal data, a Matlab toolbox.
% Main Function: playCello.m

close all
clear
clc

tic

%% Reading and filtering somatic mutations

saviTable = mutRead('../../input.savi.txt');

%% Calculating stats of known driver genes

[saviTable,mutGeneTable] = mutStats({'TP53','ATRX','IDH1','EGFR','PTEN','PIK3CA','PIK3R1','PIK3CG','PDGFRA','RB1','NF1','PTPN11','LTBP4'},saviTable);

%% Generating longitudinal mutational landscape

hland = mutLandscape(saviTable, mutGeneTable);

%% Analyzing mutation correlation

hcom = mutCorrelation(mutGeneTable);

%% Comparing mutation frequency

h3d = mutFrequency(mutGeneTable);

%% Extracting mutational signature

hsig = mutSignature(saviTable);

%% Projecting phylogenetic trees onto Moduli space

hmod = mutTreeClustering(saviTable);

%% Constructing tumor evolutionary directed graph

G = mutDirectedGraph(saviTable, true);

%% Identifying clonal switching events

hsw = mutSwitch(saviTable,'PDGFRA');

toc
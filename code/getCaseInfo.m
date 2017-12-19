function saviCase = getCaseInfo(savi)
% get information of each case and output them as Table saviCase

[uniCase, ~, rixcase] = unique(savi.CaseID); % rixcase is the case index for each variant in savi.
ncase = length(uniCase);

% Initialization
numvar = zeros(ncase,1);        % Number of variants per case
numPvar = zeros(ncase,1);       % Number of variants per primary case
numRvar = zeros(ncase,1);       % Number of variants per recurrent case
numCvar = zeros(ncase,1);       % Number of variants in common per case
isHM = zeros(ncase,1);          % Is the case hypermutated?
numGene = zeros(ncase,1);       % Number of mutated genes per case
numPcosmic = zeros(ncase,1);    % Number of cosmic mutations per primary case
numRcosmic = zeros(ncase,1);    % Number of cosmic mutations per recurrent case
numPkdriver = zeros(ncase,1);   % Number of mutated drivers per primary case
numRkdriver = zeros(ncase,1);   % Number of mutated drivers per recurrent case
maxPvaf = zeros(ncase,1);       % Max variant allele frequency per primary case
maxRvaf = zeros(ncase,1);       % Max variant allele frequency per recurrent case

for i = 1:ncase
    numvar(i) = nnz(rixcase == i);
    numPvar(i) = nnz(rixcase == i & savi.Primary_freq >= 5);
    numRvar(i) = nnz(rixcase == i & savi.Recurrent_freq >= 5);
    numCvar(i) = nnz(rixcase == i & savi.Primary_freq >= 5 & savi.Recurrent_freq >= 5);
    isHM(i) = unique(savi.isHM(rixcase == i));
    numGene(i) = numel(unique(savi.Gene_Name(rixcase == i)));
    numPcosmic(i) = nnz(~strcmp(savi.id_cosmic(rixcase == i & savi.Primary_freq >= 5), '-'));
    numRcosmic(i) = nnz(~strcmp(savi.id_cosmic(rixcase == i & savi.Recurrent_freq >= 5), '-'));
    numPkdriver(i) = nnz(savi.isKnownDriver(rixcase == i & savi.Primary_freq >= 5));
    numRkdriver(i) = nnz(savi.isKnownDriver(rixcase == i & savi.Recurrent_freq >= 5));
    maxPvaf(i) = max(savi.Primary_freq(rixcase == i));
    maxRvaf(i) = max(savi.Recurrent_freq(rixcase == i));
end

saviCase = table(uniCase, numvar, numPvar, numRvar, numCvar, isHM, numGene, ...
    numPcosmic, numRcosmic, numPkdriver, numRkdriver, maxPvaf, maxRvaf);

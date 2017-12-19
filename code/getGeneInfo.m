function saviGene = getGeneInfo(savi)
% get information of each gene

[GeneID,~,gix] = unique(savi.Gene_Name);
nGene = numel(GeneID);
if nGene ~= max(gix)
    error('Max GeneID is not equal to number of genes!')
end
uniCase = unique(savi.CaseID);
nCase = numel(uniCase);

nPatients = zeros(nGene,1);
nSamples = zeros(nGene,1);
nMutations = zeros(nGene,1);
nHMsamples = zeros(nGene,1);
nNonHMsamples = zeros(nGene,1);
nPsamples = zeros(nGene,1);
nRsamples = zeros(nGene,1);
nCsamples = zeros(nGene,1);
gPCR = zeros(nGene,nCase); % P=1;C=3;R=2
nLowMut = zeros(nGene,1);
nHighMut = zeros(nGene,1);
isSwitch = zeros(nGene,1);
SwitchCases = cell(nGene,1);
isKnownDriver = zeros(nGene,1);
ProteinLen = zeros(nGene,1);
maxSgt1maxFreq = zeros(nGene,1);
pvPRpropDiff = ones(nGene,1); % pv = 1 means insignificance
isCbioRecord = zeros(nGene,1);

for i = 1:nGene
    gi = find(gix == i);
    
    nPatients(i) = numel(unique(savi.CaseID(gi)));
    
    nMutations(i) = numel(gi);

    nHMsamples(i) = numel(unique(savi.CaseID(gix == i & savi.isHM == 1)));
    
    Psam = unique(savi.CaseID(gix == i & savi.Primary_freq > 5));   % Frequency cutoff for presence of one variant = 5?
    Rsam = unique(savi.CaseID(gix == i & savi.Recurrent_freq > 5));
    
	nPsamples(i) = numel(Psam);
    nRsamples(i) = numel(Rsam);
    nSamples(i) = nPsamples(i) + nRsamples(i);
    nNonHMsamples(i) = nRsamples(i) - nHMsamples(i);
    
    for u = 1:nPsamples(i)
        px = strcmp(uniCase, Psam{u});
        gPCR(i, px) = gPCR(i, px) + 1;
    end
	for v = 1:nRsamples(i)
        rx = strcmp(uniCase, Rsam{v});
        gPCR(i, rx) = gPCR(i, rx) + 2;
	end
    nCsamples(i) = nnz(gPCR(i,:) == 3);

    % Proportion test: [h,p,chi2stat,df] = prop_test(X,N,correct)
    [~,pvPRpropDiff(i),~,~] = prop_test([nPsamples(i),nRsamples(i)], [nCase,nCase], true);
    
    nLowMut(i) = nnz(savi.isLowImpact(gi) == 1);
    nHighMut(i) = nnz(savi.isHighImpact(gi) == 1);
    
    nonlowcase = unique(savi.CaseID(gix == i & savi.isLowImpact == 0));
    for u = 1:numel(nonlowcase)
        gsix = find(gix == i & savi.isLowImpact == 0 & strcmp(savi.CaseID, nonlowcase{u}));
        if any(savi.Primary_freq(gsix) <= 5 & savi.Recurrent_freq(gsix) >= 20)
            if any(savi.Primary_freq(gsix) >= 20 & savi.Recurrent_freq(gsix) <= 5)
                isSwitch(i) = isSwitch(i) + 1;
                SwitchCases{i,isSwitch(i)} = nonlowcase{u};
            end
        end
    end

    if any(savi.isKnownDriver(gi))
        isKnownDriver(i) = 1;
    end
    
    ProteinLen(i) = max(cell2num(strsplit(savi.Amino_Acid_length{gi(1)}, ', ')));
    
    maxprolen = max(savi.Sgt1_max_frequency(gi));
    if ~isnan(maxprolen)
        maxSgt1maxFreq(i) = maxprolen;
    end
    
    maxcbio = max(cell2num(savi.cbio_somatic_verified_mutation_count(gi)));
    if ~isnan(maxcbio)
        isCbioRecord(i) = maxcbio;
    end

end

saviGene = table(GeneID, nPatients, nSamples, nMutations, nHMsamples, ...
    nNonHMsamples, nPsamples, nRsamples, nCsamples, pvPRpropDiff, ...
    nLowMut, nHighMut, isSwitch, SwitchCases, isKnownDriver, ...
    ProteinLen, maxSgt1maxFreq, isCbioRecord, gPCR);

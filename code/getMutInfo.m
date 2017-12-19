function [saviMut,savi] = getMutInfo(savi)
% get information of each variant/mutation

bindMutID = strcat(savi.chr,':', num2str(savi.pos), savi.ref, savi.alt);
[MutID,~,mutix] = unique(bindMutID);
nMut = numel(MutID);
if nMut ~= max(mutix)
    error('Max MutID is not equal to number of mutations!')
end

nSamples = zeros(nMut,1);
chr = cell(nMut,1);
pos = zeros(nMut,1);
geneName = cell(nMut,1);
isLowImpact = zeros(nMut,1);
isHighImpact = zeros(nMut,1);
isEarly = zeros(nMut,1);
isOnlyPrimary = zeros(nMut,1);
isOnlyRecurrent = zeros(nMut,1);

savi.isLowImpact = zeros(size(savi,1),1);
savi.isHighImpact = zeros(size(savi,1),1);

for i = 1:nMut
    nSamples(i) = nnz(mutix == i);
    mx = find(mutix == i);
    
    chr{i} = savi.chr{mx(1)};
    pos(i) = savi.pos(mx(1));
    geneName{i} = savi.Gene_Name{mx(1)};
    
    if isempty(strfind(savi.Effect_Impact{mx(1)}, 'MODERATE'))
        if isempty(strfind(savi.Effect_Impact{mx(1)}, 'HIGH'))
            isLowImpact(i) = 1;
            savi.isLowImpact(mutix == i) = 1;
        elseif isempty(strfind(savi.Effect_Impact{mx(1)}, 'LOW'))
            isHighImpact(i) = 1;
            savi.isHighImpact(mutix == i) = 1;
        end
    end
    
    if savi.Primary_freq(mx(1)) >= 15 && savi.Recurrent_freq(mx(1)) >= 15
        isEarly(i) = 1;
    elseif savi.Primary_freq(mx(1)) < 5 && savi.Recurrent_freq(mx(1)) >= 15
        isOnlyPrimary(i) = 1;
    elseif savi.Primary_freq(mx(1)) >= 15 && savi.Recurrent_freq(mx(1)) < 5
        isOnlyRecurrent(i) = 1;
    end
end

saviMut = table(MutID, chr, pos, geneName, nSamples, isLowImpact, isHighImpact, isEarly, isOnlyPrimary, isOnlyRecurrent);

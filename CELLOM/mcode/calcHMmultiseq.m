function Y = calcHMmultiseq(savi)

nvar = size(savi,1);

[unicase, ~, cix] = unique(savi.CaseID);
numcases = numel(unicase);

mutload = zeros(numcases,1); 
tmz = zeros(numcases,1);

savi.isC2T = false(nvar,1);
savi.isC2T((strcmp(savi.ref,'C') & strcmp(savi.alt,'T')) ...
    | (strcmp(savi.ref,'G') & strcmp(savi.alt,'A'))) = true;

savi.isCC2TC = false(nvar,1); % tmz
savi.isCC2TC((strcmp(savi.ref,'C') & strcmp(savi.alt,'T') & strcmp(savi.varSurffix,'C')) ...
    | (strcmp(savi.ref,'G') & strcmp(savi.alt,'A') & strcmp(savi.varPrefix,'G'))) = true;

savi.isCT2TT = false(nvar,1); % age
savi.isCT2TT((strcmp(savi.ref,'C') & strcmp(savi.alt,'T') & strcmp(savi.varSurffix,'T')) ...
    | (strcmp(savi.ref,'G') & strcmp(savi.alt,'A') & strcmp(savi.varPrefix,'A'))) = true;

psedocount = 1;

for i = 1:numcases % for each sample
    psavi1 = savi(cix == i,:);
    mutload(i) = size(psavi1,1) + psedocount;
    tmz(i) = nnz(psavi1.isC2T)/mutload(i,1) ...
        + nnz(psavi1.isCC2TC)/(nnz(psavi1.isC2T) + psedocount) ...
        + sign(nnz(psavi1.isCC2TC) - nnz(psavi1.isCT2TT))*nnz(psavi1.isCT2TT)/(nnz(psavi1.isC2T) + psedocount);
end
tmz(tmz < 0) = 0;

Y = table(unicase, mutload, tmz);
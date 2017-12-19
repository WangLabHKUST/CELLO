function [CX,CY,CI] = splitGenes(X,Y,F)
% split X and Y into F subsets at random.
% the 0-1 labels in Y must be split uniformly to each subset.

ng = size(X,1);
if numel(Y) ~= ng
    error('X and Y must have the same rows, i.e., number of genes!');
end

LX = cell(2,F); % 1st row is pos; 2nd row is neg.
LY = cell(2,F);
LI = cell(2,F);

ix = find(Y == 1); % positive samples
shix = randperm(numel(ix));
bz = ceil(numel(ix)/F);
for i = 1:F
    start = 1 + (i-1)*bz;
    subset = start:min(start+bz-1,numel(ix));
	LI{1,i} = ix(shix(subset));
    LX{1,i} = X(ix(shix(subset)),:);
    LY{1,i} = Y(ix(shix(subset)));
end

ix = find(Y == 0); % negative samples
% uncomment below to select the same amount of negative samples with positive samples:
%ix = ix(randperm(numel(ix),nnz(Y))); 
shix = randperm(numel(ix));
bz = ceil(numel(ix)/F);
for i = 1:F
    start = 1 + (i-1)*bz;
    subset = start:min(start+bz-1,numel(ix));
	LI{2,i} = ix(shix(subset));
    LX{2,i} = X(ix(shix(subset)),:);
    LY{2,i} = Y(ix(shix(subset)));
end

CX = cell(F,1);
CY = cell(F,1);
CI = cell(F,1);

for i = 1:F
    CX{i} = [LX{1,i};LX{2,i}];
    CY{i} = [LY{1,i};LY{2,i}];
    CI{i} = [LI{1,i};LI{2,i}];
end

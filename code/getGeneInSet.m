function six = getGeneInSet(g, S)

C = strsplit(S, ', '); % S must be a string, like savi.Gene_Name{i}

n = length(C);

six = zeros(1,n);

for i = 1:n
    if strcmp(C{i}, g)
        six(i) = 1;
    end
end
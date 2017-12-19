function MT = countMutTypes(mutTable)
% count 6 mut types for the input mutation table.
% 6 mutational types:
% 1: A>C / T>G
% 2: A>G / T>C
% 3: A>T / T>A
% 4: C>A / G>T
% 5: C>G / G>C
% 6: C>T / G>A

MT = zeros(1,6);

MT(1) = nnz((strcmp(mutTable.ref, 'A') & strcmp(mutTable.alt, 'C')) | (strcmp(mutTable.ref, 'T') & strcmp(mutTable.alt, 'G')));
MT(2) = nnz((strcmp(mutTable.ref, 'A') & strcmp(mutTable.alt, 'G')) | (strcmp(mutTable.ref, 'T') & strcmp(mutTable.alt, 'C')));
MT(3) = nnz((strcmp(mutTable.ref, 'A') & strcmp(mutTable.alt, 'T')) | (strcmp(mutTable.ref, 'T') & strcmp(mutTable.alt, 'A')));
MT(4) = nnz((strcmp(mutTable.ref, 'C') & strcmp(mutTable.alt, 'A')) | (strcmp(mutTable.ref, 'G') & strcmp(mutTable.alt, 'T')));
MT(5) = nnz((strcmp(mutTable.ref, 'C') & strcmp(mutTable.alt, 'G')) | (strcmp(mutTable.ref, 'G') & strcmp(mutTable.alt, 'C')));
MT(6) = nnz((strcmp(mutTable.ref, 'C') & strcmp(mutTable.alt, 'T')) | (strcmp(mutTable.ref, 'G') & strcmp(mutTable.alt, 'A')));
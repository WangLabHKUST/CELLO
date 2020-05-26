function saviTable = mutRead(fileid)

saviTable = readtable(fileid);
somaticfilter = saviTable.refdepth_Blood >= 20 ...
    & saviTable.altdepth_Blood <= 1 ...
    & saviTable.Sgt1_max_frequency >= 5;
saviTable = saviTable(somaticfilter,:);
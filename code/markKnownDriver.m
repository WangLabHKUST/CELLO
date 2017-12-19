function savi = markKnownDriver(kdlist,savi,markmode)
% Add one column to savi table for known drivers

nkd = length(kdlist);
nvar = size(savi,1);

savi.isKnownDriver = zeros(nvar,1);

if markmode == 'slow'
    for i = 1:nkd
        disp(['Marking driver ',kdlist{i}])
        for k = 1:nvar
            if any(getGeneInSet(kdlist{i}, savi.Gene_Name{k}))
                savi.isKnownDriver(k) = 1;
            end
        end
    end
elseif markmode == 'fast'
    for i = 1:nkd
        disp(['Marking driver ',kdlist{i}])
        savi.isKnownDriver(strcmp(savi.Gene_Name, kdlist{i})) = 1;
    end
else
    error('unknown mark mode!')
end
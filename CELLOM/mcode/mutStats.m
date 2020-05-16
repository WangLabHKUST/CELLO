function [savi,kdtables] = mutStats(kdlist,savi)
kdtables.kdlist = kdlist;
ng = numel(kdlist);
nvar = size(savi,1);
savi.isKnownDriver = zeros(nvar,1);
savi.DriverName = cell(nvar,1);
for i = 1:ng
    savi.isKnownDriver(ismember(savi.Gene_Name, kdlist{i})) = i;
    savi.DriverName(ismember(savi.Gene_Name, kdlist{i})) = {kdlist{i}};
%     savi.isKnownDriver(~cellfun(@isempty,regexp(savi.Gene_Name, ['\<',kdlist{i},'\>']))) = i;
%     savi.DriverName(~cellfun(@isempty,regexp(savi.Gene_Name, ['\<',kdlist{i},'\>']))) = {kdlist{i}};
end

%%
[kdtables.unicase, ~, savi.caseidx] = unique(savi.CaseID);
nc = numel(kdtables.unicase); % number of patients

G = savi(savi.isKnownDriver > 0 & ~contains(savi.Effect_Impact,'LOW'),:); % Do not consider low-impact mutations.

kdtables.Pmat = zeros(nc, ng); % Mark mutation presence in primary
kdtables.Rmat = zeros(nc, ng); % Mark mutation presence in recurrence
kdtables.Cmat = zeros(nc, ng); % Mark mutation presence in common

vafcut = 5;

for i = 1:nc
    for k = 1:ng
        gsavi = G(G.caseidx == i & G.isKnownDriver == k,:);
        if ~isempty(gsavi)
            if nnz(gsavi.Primary_freq >= vafcut & gsavi.Recurrent_freq < vafcut) > 0
                kdtables.Pmat(i,k) = 1;
            end
            
            if nnz(gsavi.Primary_freq < vafcut & gsavi.Recurrent_freq >= vafcut) > 0
                kdtables.Rmat(i,k) = 1;
            end
            
            if kdtables.Pmat(i,k) == 1 && kdtables.Rmat(i,k) == 1
                kdtables.Cmat(i,k) = 1;
            end
            
            if nnz(gsavi.Primary_freq >= vafcut & gsavi.Recurrent_freq >= vafcut) > 0
                kdtables.Cmat(i,k) = 1;
            end
        end
    end
end



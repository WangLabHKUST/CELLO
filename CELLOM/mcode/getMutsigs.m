function hsig = getMutsigs(savi)
% Mutational Signature Analysis
% Temozolomide-induced signature

[unicase, ~, cix] = unique(savi.CaseID);
numcases = numel(unicase);

mutload = zeros(numcases,2); % col1: Primary; col2: Recurrence
tmz = zeros(numcases,2); % col1: Primary; col2: Recurrence

vafcut = 5;

savi = savi(ismember(savi.ref,{'A','T','C','G'}) | ismember(savi.alt,{'A','T','C','G'}),:);
nvar = size(savi,1);

savi.isC2T = false(nvar,1);
savi.isC2T((strcmp(savi.ref,'C') & strcmp(savi.alt,'T')) ...
    | (strcmp(savi.ref,'G') & strcmp(savi.alt,'A'))) = true;

savi.isCC2TC = false(nvar,1);
savi.isCC2TC((strcmp(savi.ref,'C') & strcmp(savi.alt,'T') & strcmp(savi.varSurffix,'C')) ...
    | (strcmp(savi.ref,'G') & strcmp(savi.alt,'A') & strcmp(savi.varPrefix,'G'))) = true;

savi.isCT2TT = false(nvar,1);
savi.isCT2TT((strcmp(savi.ref,'C') & strcmp(savi.alt,'T') & strcmp(savi.varSurffix,'T')) ...
    | (strcmp(savi.ref,'G') & strcmp(savi.alt,'A') & strcmp(savi.varPrefix,'A'))) = true;

psedocount = 1;

for i = 1:numcases
    % Primary
    psavi1 = savi(cix == i & savi.Blood_freq <= 1 & savi.Primary_freq >= vafcut,:);
    mutload(i,1) = size(psavi1,1) + psedocount;
    tmz(i,1) = (nnz(psavi1.isC2T) + psedocount)/mutload(i,1) + 2/3*nnz(psavi1.isCC2TC)/(nnz(psavi1.isC2T) + psedocount) + 1/3*nnz(psavi1.isCT2TT)/(nnz(psavi1.isC2T) + psedocount);
    
    % Recurrence
    rsavi1 = savi(cix == i & savi.Blood_freq <= 1 & savi.Recurrent_freq >= vafcut,:);
    mutload(i,2) = size(rsavi1,1) + psedocount;
    tmz(i,2) = (nnz(rsavi1.isC2T) + psedocount)/mutload(i,2) + 2/3*nnz(rsavi1.isCC2TC)/(nnz(rsavi1.isC2T) + psedocount) + 1/3*nnz(rsavi1.isCT2TT)/(nnz(rsavi1.isC2T) + psedocount);
end

%% tmz

hsig = figure('position',[0 0 600 600]);
hold on

scatter(log10(mutload(:,1)), tmz(:,1), ones(numcases,1)*100, repmat([1 0 0],numcases,1),'s', 'LineWidth',1.5)
scatter(log10(mutload(:,2)), tmz(:,2), ones(numcases,1)*100, repmat([0 0 0],numcases,1), 'LineWidth',1.5)

legend({'Initial','Recurrence'},'Location','southeast','Box','on','FontSize',16)

line([0 5],[1.2 1.2],'color',[.5 .5 .5],'linestyle','--','linewidth',2)
line([log10(500) log10(500)],[0.2 1.6],'color',[.5 .5 .5],'linestyle','--','linewidth',2)

xlabel('Number of Mutations (Log10 Scale)')
ylabel('Temozolomide-induced Signature')
axis square
set(gca,'tickdir','out','TickLength',[0.0075 0.0075],'fontsize',16,'box','off','linewidth',1.5)
hold off

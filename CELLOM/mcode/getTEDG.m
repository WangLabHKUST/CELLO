function G = getTEDG(savi,isdeconv)
% TEDG: Tumor Evolutionary Directed Graphs

evogenes = {'IDH1','PIK3CA','ATRX','PIK3R1','PTEN','TP53','PIK3CG','NF1','EGFR','MSH6','LTBP4'};

ng = numel(evogenes);
[~,kdstr] = markKnownDriver(evogenes,savi);

ns = numel(kdstr.unicase);
M = zeros(ns,ng);

M(kdstr.Pmat == 1) = 1;
M(kdstr.Rmat == 1) = 2;
M(kdstr.Cmat == 1) = 3;

Mtable = array2table(M);
Mtable.Properties.VariableNames = evogenes;
Mtable.Properties.RowNames = kdstr.unicase;

%% Construct the directed graph

A = zeros(ng,ng); % Adjacency Matrix
edgelabels = cell(ng,ng);

for i = 1:ns
    early = find(M(i,:) == 3);
    late = find(M(i,:) == 1 | M(i,:) == 2);
    
    if ~isempty(early) && ~isempty(late)
        A(early,late) = A(early,late) + 1;
    end
    
    for u = 1:length(early)
        uu = early(u);
        for v = 1:length(late)
            vv = late(v);
            if isempty(edgelabels{uu,vv})
                edgelabels{uu,vv} = kdstr.unicase{i};
            else
                elb = [edgelabels{uu,vv},';',kdstr.unicase{i}];
                edgelabels{uu,vv} = elb;
            end
        end
    end
    
end

%% unidirection: if A ->(3) B and B ->(1) A, then A ->(3) B

for i = 1:(ng-1)
    for j = (i+1):ng
        if A(i,j) >= A(j,i)
            A(j,i) = 0;
            edgelabels{j,i} = [];
        else
            A(i,j) = 0;
            edgelabels{i,j} = [];
        end
    end
end

%% Deconvolution

if isdeconv
    B = A + A';
    %[~,~,bv] = find(B);
    %T0 = full(mst(spones(B),'algname','kruskal','root',1,'edge_weight',-bv));
    T1 = full(adjacency(minspantree(graph(0 - B, evogenes),'Root',1)));
end

%% Start printing for Cytoscape ...

% Edge Table
fedge = fopen('../output/cyto.TEDG.edgetable.txt','w');
fprintf(fedge, 'EarlyIdx\tLateIdx\tEarlyGene\tLateGene\tNumOccurrences\tSampleIDs\n');

[u,v,ew] = find(A);
for i = 1:length(u)
    if isdeconv
        if ew(i) >= 3 || T1(u(i),v(i)) ~= 0
            fprintf(fedge,'%d\t%d\t%s\t%s\t%s\t%s\n', ...
                u(i), ...
                v(i), ...
                evogenes{u(i)}, ...
                evogenes{v(i)}, ...
                num2str(ew(i)), ...
                edgelabels{u(i),v(i)});
        end
    else
        fprintf(fedge,'%d\t%d\t%s\t%s\t%s\t%s\n', ...
            u(i), ...
            v(i), ...
            evogenes{u(i)}, ...
            evogenes{v(i)}, ...
            num2str(ew(i)), ...
            edgelabels{u(i),v(i)});
    end
end
fclose(fedge);

% Node Table

ins = sum(A,1)';    % in-degree
outs = sum(A,2);    % out-degree

pcdf = ones(ng,1); % for each gene in the network
Occurrence = zeros(ng,1); % frequently appear
for i = 1:ng
    if ins(i) < (ins(i)+outs(i))/2
        % y = binocdf(x,N,p) computes (x,y) that follow binomial dist (N,p)
        pcdf(i)  = binocdf(ins(i), ins(i)+outs(i), 0.5);
    else
        pcdf(i)  = binocdf(outs(i), ins(i)+outs(i), 0.5);
    end
    
    Occurrence(i) = nnz(M(:,i) == 1) + nnz(M(:,i) == 2) + nnz(M(:,i) == 3)*2;
end
Log2FC = log2((outs+1)./(ins+1)); % positive = early; negative = late.

fnode = fopen('../output/cyto.TEDG.nodetable.txt','w');
fprintf(fnode, 'GeneName\tP_CDF\tLog2FC\tOccurrence\n');
for i = 1:ng
    fprintf(fnode, '%s\t%f\t%f\t%f\n', evogenes{i}, pcdf(i), Log2FC(i), Occurrence(i));
end
fclose(fnode);

%% Network glance in Matlab
% NOTE: The following code works only for the default gene list with
% predetermined layout.

Edgetable = readtable('../output/cyto.TEDG.edgetable.txt');

if ~isdeconv
    Edgetable = Edgetable(Edgetable.NumOccurrences > 1 | strcmp(Edgetable.LateGene,'PIK3CG') | strcmp(Edgetable.LateGene,'PIK3R1'),:);
end

Nodetable = readtable('../output/cyto.TEDG.nodetable.txt');
G = digraph(Edgetable.EarlyIdx,Edgetable.LateIdx,Edgetable,Nodetable);

Nodetable.X = [0;-1.2;1;-1.8;-2.3;0.5;-2.3;-0.5;-3;2;3];
Nodetable.Y = [0;-1;-1;-2;-3;-3;-4;-4;-5;-5;-5];

figure('position',[1200 0 600 600])
H = plot(G,'XData',Nodetable.X,'YData',Nodetable.Y);
labelnode(H, 1:ng, Nodetable.GeneName)

fc = -5:0.01:5;
nc = length(fc);
colorscale = mycolorgrad3(nc,[0 0 0],[1 1 0],[1 0 0]);
nodecolors = zeros(G.numnodes,3);
for i = 1:G.numnodes
    [~,ix] = min(abs(fc - Nodetable.Log2FC(i)));
    nodecolors(i,:) = colorscale(ix,:);
end

H.NodeColor = nodecolors;
H.MarkerSize = Nodetable.Occurrence/3 + 5;
H.EdgeColor = 'k';

box('off')
set(gca,'visible','off')


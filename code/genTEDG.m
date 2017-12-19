function genTEDG(evoGenes, saviGene, samples)
% output:
%    cytoscape files
%      1. cyto.net
%      2. cyto.net.edge.atribute
%      3. cyto.net.node.atribute
% plot:
%    use cytoscape to load the network and all atributes

ng = numel(evoGenes);
ns = numel(samples);
gix = zeros(ng,1);
for i = 1:ng
    gix(i) = find(strcmp(saviGene.GeneID, evoGenes{i}));
end
evoSavi = saviGene(gix,:);
MP = zeros(ng,ns);
MC = zeros(ng,ns);
MR = zeros(ng,ns);
MP(evoSavi.gPCR == 1) = 1;
MC(evoSavi.gPCR == 3) = 1;
MR(evoSavi.gPCR == 2) = 1;
MP = MP';
MC = MC';
MR = MR';

%=========================================================================%

A = zeros(ng,ng); % Adjacency matrix of genes
edgeLabels = cell(ng,ng);

for i = 1:ns % for each sample
    prix = find(~~MP(i,:) | ~~MR(i,:)); % P-only or R-only gene indices
    for j = 1:numel(prix)
        cix = find(MC(i,:));
        A(prix(j), cix) = A(prix(j), ~~MC(i,:)) + 1;
        % A(j,i) += 1 if gene i is common and gene j is PR only.
        for k = 1:numel(cix) 
            % for each gene-gene interaction, add sample labels:
            if isempty(edgeLabels{prix(j),cix(k)})
                edgeLabels{prix(j),cix(k)} = samples{i};
            else
                edgeLabels{prix(j),cix(k)} = [edgeLabels{prix(j),cix(k)},';',samples{i}];
            end
        end
    end
    
end

A = A';
A = A - diag(diag(A)); % remove self-loops
edgeLabels = edgeLabels';

% Start printing ...
fedge = fopen('cyto.TEDG.edges.txt','w');
fprintf(fedge, 'GeneA\tGeneB\tWidth\tLabel\n');
[u,v] = find(A);
for i = 1:length(u)
	fprintf(fedge,'%s\t%s\t%s\t%s\n', ...
        evoGenes{u(i)}, ...
        evoGenes{v(i)}, ...
        num2str(A(u(i),v(i))), ...
        edgeLabels{u(i),v(i)});   
end
fclose(fedge);

ins = sum(A,1)';    % in-degree
outs = sum(A,2);    % out-degree

pcdf = ones(ng,1); % for each gene in the network
for i = 1:numel(pcdf)
    if ins(i) < (ins(i)+outs(i))/2
        % y = binocdf(x,N,p) computes (x,y) that follow binomial dist (N,p)
        pcdf(i)  = binocdf(ins(i), ins(i)+outs(i), 0.5);
    else
        pcdf(i)  = binocdf(outs(i), ins(i)+outs(i), 0.5);
    end
end
fc = log2((outs+1)./(ins+1)); % positive = early; negative = late.
recurrence = sum(~~(MP+MC+MR))'; % frequently appear

fnode = fopen('cyto.TEDG.nodes.txt','w');
fprintf(fnode, 'Gene\tP_CDF\tFC\tOccurrence\n');
for i = 1:ng
    fprintf(fnode, '%s\t%f\t%f\t%f\n', evoGenes{i}, pcdf(i), fc(i), recurrence(i));
end
fclose(fnode);

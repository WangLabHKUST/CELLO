function hsw = plotSwitch(geneTable,savi)
% plot switched genes

if nnz(geneTable.isSwitch) > 10
    disp('The input has more than 10 switched genes. Plotting the first 10 now ...');
end

nfig = 0;
for i = 1:size(geneTable,1)
    if nfig <= 10 && geneTable.isSwitch(i) > 0 % max(saviGene.isSwitch) = 1
        % savi report of the gene in the case.
        for k = 1:geneTable.isSwitch(i)
            gcsavi = savi(strcmp(savi.Gene_Name, geneTable.GeneID{i}) & strcmp(savi.CaseID, geneTable.SwitchCases{i,k}), :);
            nvarpos = size(gcsavi,1);
            hsw = figure;
            hold on
            for v = 1:nvarpos
                if gcsavi.Primary_freq(v) > gcsavi.Recurrent_freq(v)
                    plot([10 90], [gcsavi.Primary_freq(v), gcsavi.Recurrent_freq(v)],'r-o','LineWidth',2);
                else
                    plot([10 90], [gcsavi.Primary_freq(v), gcsavi.Recurrent_freq(v)],'k-o','LineWidth',2)
                end
            end
            xlim([0 100])
            xticks([10 90])
            xticklabels({'Primary', 'Relapsed'})
            ylim([0 100])
            ylabel('Mutation Frequency in Tumor')
            legend(gcsavi.Amino_Acid_Change, 'Location','north')
            title([geneTable.GeneID{i},' in ',geneTable.SwitchCases{i,k}])
            box on
            hold off
            nfig = nfig + 1;
        end
    end
end
function hsw = mutSwitch(savi,gene)
% Plot clonal swithcing event

allgenelist = unique(savi.Gene_Name);

if ~strcmp(gene,allgenelist)
    error('The input gene is not found in the mutant gene list!')
else
    savi = mutStats({gene},savi);
end

genelist = gene;
ng = numel(genelist);

esavi = savi(savi.isKnownDriver > 0,:);
ne = size(esavi,1);
vafcut = 5;
esavi.primary = false(ne,1);
esavi.primary(esavi.Primary_freq >= vafcut) = true;
esavi.recurrence = false(ne,1);
esavi.recurrence(esavi.Recurrent_freq >= vafcut) = true;

[unicase,~,cix] = unique(esavi.CaseID);
ns = numel(unicase);

W = cell(100,1);
nw = 0;
for i = 1:ns
    for k = 1:ng
        gesavi = esavi(cix == i & esavi.isKnownDriver == k,:);
        n = size(gesavi,1);
        if n > 1
            nw = nw + 1;
            W{nw} = gesavi;
        end
    end
end
W = W(1:nw);

f = 0;
for i = 1:nw
    if all(contains(W{i}.Gene_Name,gene)) ...
            && any(W{i}.primary & ~W{i}.recurrence) ...
            && any(~W{i}.primary & W{i}.recurrence)
        f = f + 1;
        hsw = figure('position',[2000, 600, 700, 600]);
        hold on
        for k = 1:size(W{i},1)
            y1 = W{i}.Primary_freq(k);
            y2 = W{i}.Recurrent_freq(k);
            if y1 >= y2
                plot([1, 2], [y1, y2],'ro-','linewidth',1.5,'markersize',10,'markerfacecolor','r')
                text(2.03,y2,W{i}.Amino_Acid_Change{k},'HorizontalAlignment','left', 'color','k','fontsize', 15)
            else
                plot([1, 2], [y1, y2],'ks-','linewidth',1.5,'markersize',10,'markerfacecolor','k')
                text(2.03,y2,W{i}.Amino_Acid_Change{k},'HorizontalAlignment','left', 'color','k','fontsize', 15)
            end
        end
        
        legend(W{i}.Amino_Acid_Change, 'Location', 'north')
        
        xlim([0.8 2.2])
        xticks([1 2])
        xticklabels({'Primary', 'Recurrence'})
        
        ylim([-5 105])
        ylabel('VAF')
        
        title([W{i}.CaseID{1},':',gene])
        
        set(gca,'tickdir','out','TickLength',[0.0075 0.0075],'fontsize',16,'box','off','linewidth',1.5)
        hold off
        
    end
end

if f == 0
    disp('The input gene has no clonal switching event!')
    hsw = 0;
end




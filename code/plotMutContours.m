function nmutPCR = plotMutContours(caseid,pfreq,rfreq,freq_cutoff,doPlot)
% plot contours for all patients, and put the figures into figure_dir

if doPlot == 'plot'
    figdir = 'outMutContours';
    mkdir(figdir)
end

[uniCase,~,ixcase] = unique(caseid);
nCase = max(ixcase);
if numel(uniCase) ~= nCase
    error('unique case ids are wrong!');
end

% no. of muts in Primary, Common, and Recurrent.
PCR = zeros(numel(uniCase),3); 

for i = 1:nCase
    disp(['Plotting Mutation Contours for Patient ',uniCase{i}])
    
    PCR(i,1) = nnz(ixcase == i & pfreq >= freq_cutoff & rfreq < freq_cutoff);   % primary
    PCR(i,2) = nnz(ixcase == i & pfreq >= freq_cutoff & rfreq >= freq_cutoff);  % common
    PCR(i,3) = nnz(ixcase == i & pfreq < freq_cutoff & rfreq >= freq_cutoff);	% recurrent
    
    if strcmp(doPlot,'plot')
        % Start plotting ...
        h = figure('visible','off');
        hold on
        line([0 100],[0 100],'linewidth',1,'color',[.5 .5 .5])
        title(['Mutation Contour of Patient ',uniCase{i}],'FontSize',20,'fontweight','bold')
        line([0 0],[5 100],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')
        line([5 5],[5 100],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')
        line([0 5],[5 5],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')
        line([0 5],[100 100],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')

        line([5 100],[0 0],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')
        line([5 100],[5 5],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')
        line([5 5],[0 5],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')
        line([100 100],[0 5],'linewidth',1,'color',[.5 .5 .5],'linestyle','--')

        prFreq = [pfreq(ixcase == i & ~isnan(pfreq + rfreq)), rfreq(ixcase == i & ~isnan(pfreq + rfreq))];
        plot(prFreq(:,1), prFreq(:,2), '.', 'color', 'k');
        dscatter(prFreq(:,1), prFreq(:,2), 'PLOTTYPE','contour');
        dscatter(prFreq(:,1), prFreq(:,2));

        text(95,2,num2str(PCR(i,1)),'FontSize',16,'fontweight','bold')
        text(95,95,num2str(PCR(i,2)),'FontSize',16,'fontweight','bold')
        text(0,95,num2str(PCR(i,3)),'FontSize',16,'fontweight','bold')

        set(gca, 'xlim',[-5 105], 'ylim',[-5 105], 'linewidth',1, ...
            'tickdir','out', 'FontSize',16, 'FontWeight','bold', ...
            'xtick',[0,15,50,100], 'ytick',[0,15,50,100])

        axis square;
        box on;

        xlabel('VAF of Untreated Tumor', 'FontSize',20, 'FontWeight','bold');
        ylabel('VAF of Recurrent Tumor', 'FontSize',20, 'FontWeight','bold');

        save2pdf([figdir,'/',uniCase{i},'.pdf'])
        close(h);
    end
end

nmutPCR = table(uniCase,PCR);

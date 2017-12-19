function [auc,fpr,tpr] = calcauroc(label,prediction,doPlot,modelname)

[tpr, fpr, c] = newroc(label',prediction' + rand(numel(prediction),1)'*eps);

if strcmp(doPlot,'plot1')
    if strcmp(modelname, 'Random Forest')
        plot(fpr, tpr, 'b-', 'LineWidth', 2)
    elseif strfind(modelname, 'SVM')
        plot(fpr, tpr, 'r-', 'LineWidth', 2)
    else
        error('unknown model and color')
    end
elseif strcmp(doPlot,'plotall')
	figure;
    hold on
    plot(fpr, tpr, 'b-', 'LineWidth', 2)
    plot([0 1], [0 1], 'k--', 'LineWidth', 2)
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    legend({modelname,'Random Control'}, 'Location','southeast')
    hold off
end

auc = 0;
for i = 1:size(c,2)-1
    auc = auc + (fpr(i+1) - fpr(i))*(tpr(i) + tpr(i+1))/2;
end
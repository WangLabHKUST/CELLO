function annotationAxes = CreateAnnotaionAxes(mainAxes,mainParent)

    %Create Annotation Axis, Remove graphics objects, YAxis annotations
    %(except YLabel) and make background transparent
    annotationAxes = copyobj(mainAxes,mainParent);
    
    set(annotationAxes,'XLimMode','Manual');
    
    children = get(annotationAxes,'Children');
    for i = 1:numel(children)
       delete(children(i)); 
    end

    %Save the yLabelpostion because it will move when we delete yAxis
    %ticks
    yLabel = get(annotationAxes,'YLabel');
    yLabelPosition = get(yLabel,'Position');
    
    set(annotationAxes,'YGrid' ,'off', ...
        'YMinorGrid', 'off', ...
        'YMinorTick','off', ...
        'YTick', [], ...
        'YTickLabel', []);
    
    %Restore the pevious label postition
    set(yLabel,'Position',yLabelPosition);
end
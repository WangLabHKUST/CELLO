function list = SetVisibleOff(handle, list)
    if (strcmp(get(handle,'Visible'),'on'))
        set(handle,'Visible','off');
        list = [list handle];
    end
    
    children = get(handle,'Children');
    for i = 1:numel(children)
        list = SetVisibleOff(children(i),list);
    end
end
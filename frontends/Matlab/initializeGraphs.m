function initializeGraphs(handles)
    global baseDir;
    axes(handles.axes1);
    plot(zeros(1));
    title('Rotation (x:blue y:green z:red)');
    xlabel('volume');
    ylabel('radian');
    
    axes(handles.axes4);
    plot(zeros(1));
    title('Translation (x:blue y:green z:red)');
    xlabel('volume');
    ylabel('mm');

    axes(handles.axes5);
    plot(zeros(1));
    title('Root mean square error');
    xlabel('volume');
    ylabel('a.u.');

    axes(handles.axes6);
    imshow(sprintf('%sFix.JPG', baseDir));
end

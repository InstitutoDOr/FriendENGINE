function updateGraphs(handles, friendObj)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global baseDir;
global names;

xmax=max(100, size(friendObj.rotationx, 2));
xmin=max(1, xmax-100);

axes(handles.axes1);
cla;
hold on, plot( friendObj.rotationx, 'b'), plot( friendObj.rotationy, 'g'), plot( friendObj.rotationz, 'r'), hold off
xlim([xmin xmax]);

axes(handles.axes4);
cla;
hold on, plot( friendObj.translationx, 'b'), plot( friendObj.translationy, 'g'), plot( friendObj.translationz, 'r'), hold off
xlim([xmin xmax]);

axes(handles.axes5);
cla;
hold on, plot( friendObj.rms, 'b'), hold off
xlim([xmin xmax]);

if (friendObj.pipelineType == 2)
   axes(handles.axes6);
   if (friendObj.percentage > 1) 
       friendObj.percentage = 1;
   end;
   if (friendObj.percentage < 0) 
       friendObj.percentage = 0;
   end;
   imageNumber=ceil(size(names, 2)*friendObj.percentage);
   if (imageNumber ==0) 
       imageNumber=1;
   end;
   imshow(sprintf('%sfigures/%s', baseDir, names{imageNumber}));
end;

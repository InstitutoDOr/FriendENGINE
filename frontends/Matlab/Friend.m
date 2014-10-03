function varargout = Friend(varargin)
% FRIEND MATLAB code for Friend.fig
%      FRIEND, by itself, creates a new FRIEND or raises the existing
%      singleton*.
%
%      H = FRIEND returns the handle to a new FRIEND or the handle to
%      the existing singleton*.
%
%      FRIEND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRIEND.M with the given input arguments.
%
%      FRIEND('Property','Value',...) creates a new FRIEND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Friend_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Friend_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Friend

% Last Modified by GUIDE v2.5 17-Sep-2014 09:08:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Friend_OpeningFcn, ...
                   'gui_OutputFcn',  @Friend_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function timerFunction(obj, event)
global guiObject;
hObject = guiObject;
handles = guidata(hObject);
processPhase(handles);

% --- Executes just before Friend is made visible.
function Friend_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Friend (see VARARGIN)

% Choose default command line output for Friend
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% graphs variables
global rotationx;
global rotationy;
global rotationz;

global translationx;
global translationy;
global translationz;

global rms;


% internal variable
global baseDir;
global timerObj;
global guiObject;
global names;
global percentage;
global actualVolume;
global phase;

% comunication protocol variables
global mainThread;
global responseThread;
global sessionID;


% initializing variables
phase = 0;
percentage=0;

host = '127.0.0.1';
port = 5678;
timeOut = 10;
mainThread = tcpip(host, port);
set(mainThread, 'TimeOut', timeOut);

responseThread = tcpip(host, port);
set(responseThread, 'TimeOut', timeOut);

sessionID = '';
actualVolume=0;

guiObject = hObject;

[pathstr,name,ext] = fileparts(get(hObject,'FileName'));
baseDir=sprintf('%s/', pathstr);
files = dir(sprintf('%s/figures/*.jpg', baseDir));
names = {files(:).name};

rotationx=[];
rotationy=[];
rotationz=[];

translationx=[];
translationy=[];
translationz=[];

rms=[];
if strcmp(get(hObject,'Visible'),'off')
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
    
    timerObj = timer;
    timerObj.ExecutionMode = 'fixedDelay';
    timerObj.TimerFcn = @timerFunction;
    timerObj.Period = 2/10; % one tenth of the TR
end;


% --- Outputs from this function are returned to the command line.
function varargout = Friend_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create session
function createSession()
global mainThread;
global sessionID;
fprintf(mainThread, 'NEWSESSION');
% reading session id
sessionID = fgetl(mainThread);
% reading ok
response=fgetl(mainThread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sending plugin information
function setPlugInInformation()
global mainThread;

% sending the PLUGIN command and parameters
fprintf(mainThread, 'PLUGIN');
fprintf(mainThread, 'libROI');
fprintf(mainThread, 'no');
fprintf(mainThread, 'processROI');
fprintf(mainThread, 'initializeROIProcessing');
fprintf(mainThread, 'finalizeProcessing');
fprintf(mainThread, 'no');
fprintf(mainThread, 'no');
% getting the acknowledge
response=fgetl(mainThread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get a process terminated ? response
function [response] = getResponse(varargin)
global responseThread;
global sessionID;

fopen(responseThread);
fprintf(responseThread, 'SESSION');
fprintf(responseThread, '%s', sessionID);
response=fgetl(responseThread);
if (strcmp(response, 'OK') == 1)
   fprintf(responseThread, '%s', varargin{1});
   if (nargin > 1)
      fprintf(responseThread, '%d', varargin{2}); 
   end;
   response=fgetl(responseThread);
   if (nargin > 1)
      fgetl(responseThread);
   end;
end;
fclose(responseThread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the feedback value. The feedback value was previously calculated and
% stored in a session variable in Friend Engine. This just opens another
% thread, referencing the session to get these values.
function [class, percentage] = getFeedbackValue()
global responseThread;
global sessionID;
global actualVolume;

    fopen(responseThread);
    fprintf(responseThread, 'SESSION');
    fprintf(responseThread, '%s', sessionID);
    response=fgetl(responseThread);
    fprintf(responseThread, 'TEST'); 
    fprintf(responseThread, '%d', actualVolume);
    class=str2double(fgetl(responseThread));
    percentage=str2double(fgetl(responseThread));
    response=fgetl(responseThread);
    fclose(responseThread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes a phase in the Friend Engine Pipeline : (PREPROC, FEEDBACK,
% query for motion correction parameters and results
function processPhase(handles)
global phase;
global actualVolume;

global rotationx;
global rotationy;
global rotationz;

global translationx;
global translationy;
global translationz;

global rms;
global percentage;
global sessionID;

global mainThread;
global responseThread;
global pipelineType;
global timerObj;

if (phase == 1)
    % sending PREPROC command
    fprintf(mainThread, 'NBPREPROC');
    response=fgetl(mainThread);
    phase = 15;
end

if (phase == 15)
    % sees if PREPROC terminates
    response=getResponse('PREPROC');
    if (strcmp(response, 'OK') == 1)
       phase = 2;
    end;
end

if (phase == 2)
    % sending FEEDBACK command
    fprintf(mainThread, 'NBFEEDBACK');
    response=fgetl(mainThread);
    actualVolume=1;
    phase = 25;
end

if (phase == 25)
    % retrieving Graph Parameters
    response=getResponse('GRAPHPARS', actualVolume);
    if (strcmp(response, 'END') == 1)
       stop(timerObj);
       
       if (pipelineType == 1)
          fprintf(mainThread, 'GLM');
          response=fgetl(mainThread);
           
          fprintf(mainThread, 'FEATURESELECTION');
          response=fgetl(mainThread);
       end;
       fprintf(mainThread, 'ENDSESSION');
       fprintf(mainThread, '%s', sessionID);
       response=fgetl(mainThread);
       phase = 100;
       return;
    end;
    
    if (strcmp(response, 'NOK.') == 0)
       tokens=regexp(response, ';', 'split');
       if (size(tokens, 2) == 9)

            % getting the feedback
            [class, percentage] = getFeedbackValue();

            % incrementing the to be processed volume
            actualVolume = actualVolume + 1;

            % updating graph vars
            rotationx = [rotationx str2double(tokens{3})];
            rotationy = [rotationy str2double(tokens{4})];
            rotationz = [rotationz str2double(tokens{5})];

            translationx = [translationx str2double(tokens{6})];
            translationy = [translationy str2double(tokens{7})];
            translationz = [translationz str2double(tokens{8})];

            rms = [rms str2double(tokens{9})];

            % showing the graph
            updateGraphs(handles);
        end;
    end;
    fclose(responseThread);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Initiates the Feedback
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mainThread;global pipelineType;

global phase;


fopen(mainThread);

% create session
createSession();

% sending plugin information
setPlugInInformation();

phase = 1;

global pipelineType;
pipelineType = 1;

global timerObj;
start(timerObj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Initiates the Feedback
function pushbutton6_Callback(hObject, eventdata, handles)
global mainThread;
global phase;


fopen(mainThread);

% create session
createSession();

% sending plugin information
setPlugInInformation();

phase = 1;

global pipelineType;
pipelineType = 2;

% changing the mask type to native space
fprintf(mainThread, 'SET');
fprintf(mainThread, 'ActivationLevelMaskType');
fprintf(mainThread, '1');
response=fgetl(mainThread);

% changing the mask to the generated by the funcional localizer run
fprintf(mainThread, 'SET');
fprintf(mainThread, 'ActivationLevelMask');
fprintf(mainThread, 'glmdirtstats_features_RUN01_bin');
response=fgetl(mainThread);

% changing the mask to the generated by the funcional localizer run
fprintf(mainThread, 'SET');
fprintf(mainThread, 'Prefix');
fprintf(mainThread, 'outputdirRUN02\DRIN-');
response=fgetl(mainThread);

global timerObj;
start(timerObj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateGraphs(handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rotationx;
global rotationy;
global rotationz;

global translationx;
global translationy;
global translationz;

global rms;

global baseDir;
global percentage;

xmax=max(100, size(rotationx, 2));
xmin=max(1, xmax-100);

axes(handles.axes1);
cla;
hold on, plot( rotationx, 'b'), plot( rotationy, 'g'), plot( rotationz, 'r'), hold off
xlim([xmin xmax]);

axes(handles.axes4);
cla;
hold on, plot( translationx, 'b'), plot( translationy, 'g'), plot( translationz, 'r'), hold off
xlim([xmin xmax]);

axes(handles.axes5);
cla;
hold on, plot( rms, 'b'), hold off
xlim([xmin xmax]);

global pipelineType;
if (pipelineType == 2)
   global names;
   axes(handles.axes6);
   if (percentage > 1) 
       percentage = 1;
   end;
   if (percentage < 0) 
       percentage = 0;
   end;
   imageNumber=ceil(size(names, 2)*percentage);
   if (imageNumber ==0) 
       imageNumber=1;
   end;
   imshow(sprintf('%sfigures/%s', baseDir, names{imageNumber}));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU IMPLEMENTEI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global timerObj;
stop(timerObj);
delete(timerObj);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function varargout = FriendMotor(varargin)
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
global timerObj;
global friendObj;

hObject = guiObject;
handles = guidata(hObject);
friendObj = processPhase(friendObj);
% showing the graph
updateGraphs(handles, friendObj);
if (friendObj.phase == 100)
    stop(timerObj);
end;

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


% internal variable
global timerObj;
global guiObject;

global baseDir;
global names;
global friendObj;


guiObject = hObject;


if strcmp(get(hObject,'Visible'),'off')
    % initialize Friend Object Structure
    host = '127.0.0.1';
    friendObj = initializeFriendObject(host, 5678);

    % initialize path variables
    [pathstr,name,ext] = fileparts(get(hObject,'FileName'));
    baseDir=sprintf('%s/', pathstr);
    files = dir(sprintf('%s/figures/*.jpg', baseDir));
    names = {files(:).name};
    
    initializeGraphs(handles);
    timerObj = timer;
    timerObj.ExecutionMode = 'fixedDelay';
    timerObj.TimerFcn = @timerFunction;
    timerObj.Period = 1; 
end;


% --- Outputs from this function are returned to the command line.
function varargout = Friend_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Initiates the Feedback
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global friendObj;


fopen(friendObj.mainThread);

% create session
friendObj.sessionID = createSession(friendObj);

% changing the MNI mask
configureVariable(friendObj, 'MNIMask', 'studydirhmat_spm_final.nii');   

% changing the MNI template
configureVariable(friendObj, 'MNITemplate', 'studydirMNI152_T1_1mm_brain.nii.gz');

% sending plugin information
friendObj = setPluginInformation(friendObj, 3);

friendObj.phase = 1;
friendObj.pipelineType = 1;

global timerObj;
start(timerObj);

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


% --- Initiates the Feedback
function pushbutton6_Callback(hObject, eventdata, handles)
    
global friendObj;

fopen(friendObj.mainThread);

% create session
friendObj.sessionID = createSession(friendObj);

friendObj.phase = 1;

friendObj.pipelineType = 2;
friendObj.rotationx = [];
friendObj.rotationy = [];
friendObj.rotationz = [];
    
friendObj.translationx = [];
friendObj.translationy = [];
friendObj.translationz = [];
friendObj.rms = [];

% changing the MNI mask
configureVariable(friendObj, 'MNIMask', 'studydirhmat_spm_final.nii');   

% changing the MNI template
configureVariable(friendObj, 'MNITemplate', 'studydirMNI152_T1_1mm_brain.nii.gz');

% changing the mask to the generated by the funcional localizer run
configureVariable(friendObj, 'Prefix', fullfile('outputdirRUN02', 'DRIN-'));

% sending plugin information
friendObj = setPluginInformation(friendObj, 3);

global timerObj;
start(timerObj);


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

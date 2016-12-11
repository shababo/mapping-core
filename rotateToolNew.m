function varargout = rotateToolNew(varargin)
% ROTATETOOLNEW MATLAB code for rotateToolNew.fig
%      ROTATETOOLNEW, by itself, creates a new ROTATETOOLNEW or raises the existing
%      singleton*.
%
%      H = ROTATETOOLNEW returns the handle to a new ROTATETOOLNEW or the handle to
%      the existing singleton*.
%
%      ROTATETOOLNEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROTATETOOLNEW.M with the given input arguments.
%
%      ROTATETOOLNEW('Property','Value',...) creates a new ROTATETOOLNEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rotateToolNew_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rotateToolNew_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rotateToolNew

% Last Modified by GUIDE v2.5 06-Dec-2016 12:47:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rotateToolNew_OpeningFcn, ...
                   'gui_OutputFcn',  @rotateToolNew_OutputFcn, ...
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


% --- Executes just before rotateToolNew is made visible.
function rotateToolNew_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rotateToolNew (see VARARGIN)

% Choose default command line output for rotateToolNew
handles.output = hObject;
handles.merge_image = varargin{1};
showIM(hObject, eventdata, handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rotateToolNew wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rotateToolNew_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in minus1.
function minus1_Callback(hObject, eventdata, handles)
% hObject    handle to minus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta=str2double(get(handles.edit1,'String'));
theta=theta-1;
set(handles.edit1,'String',num2str(theta));
showIM(hObject, eventdata, handles)

% --- Executes on button press in minuspointone.
function minuspointone_Callback(hObject, eventdata, handles)
% hObject    handle to minuspointone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta=str2double(get(handles.edit1,'String'));
theta=theta-0.1;
set(handles.edit1,'String',num2str(theta));
showIM(hObject, eventdata, handles)

% --- Executes on button press in pluspointone.
function pluspointone_Callback(hObject, eventdata, handles)
% hObject    handle to pluspointone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta=str2double(get(handles.edit1,'String'));
theta=theta+0.1;
set(handles.edit1,'String',num2str(theta));
showIM(hObject, eventdata, handles)


% --- Executes on button press in plus1.
function plus1_Callback(hObject, eventdata, handles)
% hObject    handle to plus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta=str2double(get(handles.edit1,'String'));
theta=theta+1;
set(handles.edit1,'String',num2str(theta));
showIM(hObject, eventdata, handles)


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
showIM(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function showIM(hObject, eventdata, handles)

theta=str2double(get(handles.edit1,'String'));
assignin('base','theta',theta)
axes(handles.axes1)
imagesc(imrotate(handles.merge_image,theta,'bilinear','loose'))
axis image


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','moveon',1)
% Hint: delete(hObject) closes the figure
delete(hObject);

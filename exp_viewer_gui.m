function varargout = exp_viewer_gui(varargin)
% EXP_VIEWER_GUI MATLAB code for exp_viewer_gui.fig
%      EXP_VIEWER_GUI, by itself, creates a new EXP_VIEWER_GUI or raises the existing
%      singleton*.
%
%      H = EXP_VIEWER_GUI returns the handle to a new EXP_VIEWER_GUI or the handle to
%      the existing singleton*.
%
%      EXP_VIEWER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXP_VIEWER_GUI.M with the given input arguments.
%
%      EXP_VIEWER_GUI('Property','Value',...) creates a new EXP_VIEWER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exp_viewer_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exp_viewer_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exp_viewer_gui

% Last Modified by GUIDE v2.5 16-Apr-2015 23:04:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exp_viewer_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @exp_viewer_gui_OutputFcn, ...
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


% --- Executes just before exp_viewer_gui is made visible.
function exp_viewer_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exp_viewer_gui (see VARARGIN)

% Choose default command line output for exp_viewer_gui
handles.output = hObject;

handles.data = varargin{1};

draw_plot(handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes exp_viewer_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = exp_viewer_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)plot(data.time,sweeps{trace_ind}(:,1),'k');
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_trial = str2num(get(handles.trial_number,'String'))-1;
if new_trial < 1
    new_trial = handles.data.n_trials;
end

set(handles.trial_number,'String',num2str(new_trial));
guidata(hObject,handles)
draw_plot(handles)


% --- Executes on button press in forward.
function forward_Callback(hObject, eventdata, handles)
% hObject    handle to forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_trial = str2num(get(handles.trial_number,'String'))+1;
if new_trial > handles.data.n_trials
    new_trial = 1;
end

set(handles.trial_number,'String',num2str(new_trial));
guidata(hObject,handles)
draw_plot(handles)


function draw_plot(handles)

trace_ind = str2num(get(handles.trial_number,'String'));

axes(handles.data_axes)
plot(handles.data.time,handles.data.traces(trace_ind,:),'k');
hold on;
scatter(handles.data.time(find(handles.data.spikes(trace_ind,:))),100*ones(1,sum(handles.data.spikes(trace_ind,:))),20*ones(1,sum(handles.data.spikes(trace_ind,:))),'filled');
hold on;
scatter(handles.data.time(find(handles.data.stims(trace_ind,:))),120*ones(1,sum(handles.data.stims(trace_ind,:))),20*ones(1,sum(handles.data.stims(trace_ind,:))),'filled');
hold off;

ylim([-100 150])


function trial_number_Callback(hObject, eventdata, handles)
% hObject    handle to trial_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trial_number as text
%        str2double(get(hObject,'String')) returns contents of trial_number as a double


% --- Executes during object creation, after setting all properties.
function trial_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







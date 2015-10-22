function varargout = trial_by_trial_detection_gui(varargin)
% TRIAL_BY_TRIAL_DETECTION_GUI MATLAB code for trial_by_trial_detection_gui.fig
%      TRIAL_BY_TRIAL_DETECTION_GUI, by itself, creates a new TRIAL_BY_TRIAL_DETECTION_GUI or raises the existing
%      singleton*.
%
%      H = TRIAL_BY_TRIAL_DETECTION_GUI returns the handle to a new TRIAL_BY_TRIAL_DETECTION_GUI or the handle to
%      the existing singleton*.
%
%      TRIAL_BY_TRIAL_DETECTION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIAL_BY_TRIAL_DETECTION_GUI.M with the given input arguments.
%
%      TRIAL_BY_TRIAL_DETECTION_GUI('Property','Value',...) creates a new TRIAL_BY_TRIAL_DETECTION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trial_by_trial_detection_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trial_by_trial_detection_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trial_by_trial_detection_gui

% Last Modified by GUIDE v2.5 14-Jul-2015 12:46:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trial_by_trial_detection_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @trial_by_trial_detection_gui_OutputFcn, ...
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


% --- Executes just before trial_by_trial_detection_gui is made visible.
function trial_by_trial_detection_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trial_by_trial_detection_gui (see VARARGIN)

% Choose default command line output for trial_by_trial_detection_gui
handles.output = hObject;

handles.data.traces = varargin{1};
handles.data.n_trials = size(varargin{1},1);
handles.data.denoised_traces = cell(1,length(varargin{2}));
for i = 1:length(varargin{2})
    handles.data.denoised_traces{i} = varargin{2}(i).trials.curves{varargin{2}(i).min_err_ind};
end

draw_plot(handles)
set(handles.num_traces,'String',['Num Traces: ' num2str(handles.data.n_trials)])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trial_by_trial_detection_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trial_by_trial_detection_gui_OutputFcn(hObject, eventdata, handles) 
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

start_ind = 1; %2450
trace_len = length(handles.data.traces(trace_ind,:));

axes(handles.data_axes)
plot((0:trace_len-1)/20000,handles.data.traces(trace_ind,:),'c');
hold on;
plot((0:trace_len-1)/20000,-(handles.data.denoised_traces{trace_ind}-handles.data.denoised_traces{trace_ind}(1)) + ...
    handles.data.traces(trace_ind,1) + 25,'b','LineWidth',2)
% hold on;
% plot((0:trace_len-1)/20000,handles.data.sweeps{trace_ind}(:,2)/10,'r');
% hold on;
% scatter(handles.data.time(find(handles.data.spikes(trace_ind,:))),100*ones(1,sum(handles.data.spikes(trace_ind,:))),20*ones(1,sum(handles.data.spikes(trace_ind,:))),'filled');
% hold on;
% scatter(handles.data.time(find(handles.data.stims(trace_ind,:))),120*ones(1,sum(handles.data.stims(trace_ind,:))),20*ones(1,sum(handles.data.stims(trace_ind,:))),'filled');
hold off;
axis tight
ylim([-100 40])
ylabel('Current (pA)'); xlabel('Time (sec)')


function trial_number_Callback(hObject, eventdata, handles)
% hObject    handle to trial_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trial_number as text
%        str2double(get(hObject,'String')) returns contents of trial_number as a double

new_trial = str2double(get(hObject,'String'));
if isempty(new_trial) || new_trial < 1 || new_trial > handles.data.n_trials
    set(hObject,'String',1)
end
new_trial = str2double(get(hObject,'String'));

guidata(hObject,handles)
draw_plot(handles)
    


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

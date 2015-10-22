function varargout = target_by_target_gui(varargin)
% TARGET_BY_TARGET_GUI MATLAB code for target_by_target_gui.fig
%      TARGET_BY_TARGET_GUI, by itself, creates a new TARGET_BY_TARGET_GUI or raises the existing
%      singleton*.
%
%      H = TARGET_BY_TARGET_GUI returns the handle to a new TARGET_BY_TARGET_GUI or the handle to
%      the existing singleton*.
%
%      TARGET_BY_TARGET_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TARGET_BY_TARGET_GUI.M with the given input arguments.
%
%      TARGET_BY_TARGET_GUI('Property','Value',...) creates a new TARGET_BY_TARGET_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before target_by_target_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to target_by_target_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help target_by_target_gui

% Last Modified by GUIDE v2.5 24-Jun-2015 08:47:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @target_by_target_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @target_by_target_gui_OutputFcn, ...
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


% --- Executes just before target_by_target_gui is made visible.
function target_by_target_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to target_by_target_gui (see VARARGIN)

% Choose default command line output for target_by_target_gui
handles.output = hObject;

handles.data = varargin{1};

draw_plot(handles)
set(handles.num_traces,'String',['Num Traces: ' num2str(handles.data.n_targets)])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes target_by_target_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = target_by_target_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)plot(data.time,sweeps{target_ind}(:,1),'k');
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_trial = str2num(get(handles.trial_number,'String'))-1;
if new_trial < 1
    new_trial = handles.data.n_targets;
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
if new_trial > handles.data.n_targets
    new_trial = 1;
end

set(handles.trial_number,'String',num2str(new_trial));
guidata(hObject,handles)
draw_plot(handles)


function draw_plot(handles)

target_ind = str2num(get(handles.trial_number,'String'));
these_traces = zeros([size(handles.data.traces{target_ind}) 2]);
these_traces(:,:,1) = handles.data.traces{target_ind};

results = handles.data.results{target_ind};
for i=1:size(these_traces,1)

    these_traces(i,:,2) = -results(i).trials.curves{results(i).min_err_ind} ...
        + results(i).trials.curves{results(i).min_err_ind}(1) + median(these_traces(i,1,1));
end


features = zeros(0,4);
for i=1:length(handles.data.feature_mats{target_ind})
    if ~isempty(handles.data.feature_mats{target_ind}{i})
        features = [features; handles.data.feature_mats{target_ind}{i}];
    end
end

axes(handles.traces_axes)
% PUT MIN ERR CURVES OVERLAID... WITH COLORED EVENTS??
colorsettings = zeros(size(these_traces,1),3,2);
colorsettings(:,3,2) = 1.0;
plot_trace_stack_multi(these_traces,zeros(size(these_traces)),[],colorsettings,[])
axis tight



feature_labels = {'amplitude','tau1','tau2','delay'};

axes(handles.scatter_axes)
[h,ax] = plotmatrix(features,'o');                        % create a 4 x 4 matrix of plots
for i = 1:4                                       % label the plots
  xlabel(ax(4,i), feature_labels{i})
  ylabel(ax(i,1), feature_labels{i})
end

  set(handles.figure1,'NextPlot','add')



% ylim([-150 0])


function trial_number_Callback(hObject, eventdata, handles)
% hObject    handle to trial_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trial_number as text
%        str2double(get(hObject,'String')) returns contents of trial_number as a double

new_trial = str2double(get(hObject,'String'));
if isempty(new_trial) || new_trial < 1 || new_trial > handles.data.n_targets
    set(hObject,'String',1)
end
new_trial = str2double(get(hObject,'String'));

guidata(hObject,handles)
draw_plot(handles)
    


% --- Executes during object creation, after setting all properties.
function trial_number_CreateFcn(hObjoect, eventdata, handles)
% hObject    handle to trial_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

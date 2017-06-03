function varargout = trial_by_trial_gui(varargin)
% TRIAL_BY_TRIAL_GUI MATLAB code for trial_by_trial_gui.fig
%      TRIAL_BY_TRIAL_GUI, by itself, creates a new TRIAL_BY_TRIAL_GUI or raises the existing
%      singleton*.
%
%      H = TRIAL_BY_TRIAL_GUI returns the handle to a new TRIAL_BY_TRIAL_GUI or the handle to
%      the existing singleton*.
%
%      TRIAL_BY_TRIAL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIAL_BY_TRIAL_GUI.M with the given input arguments.
%
%      TRIAL_BY_TRIAL_GUI('Property','Value',...) creates a new TRIAL_BY_TRIAL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trial_by_trial_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trial_by_trial_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trial_by_trial_gui

% Last Modified by GUIDE v2.5 19-Apr-2017 11:37:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trial_by_trial_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @trial_by_trial_gui_OutputFcn, ...
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


% --- Executes just before trial_by_trial_gui is made visible.
function trial_by_trial_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trial_by_trial_gui (see VARARGIN)

% Choose default command line output for trial_by_trial_gui
handles.output = hObject;

handles.data = varargin{1};

draw_plot(handles)
set(handles.num_traces,'String',['Num Traces: ' num2str(handles.data.n_trials)])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trial_by_trial_gui wait for user response (see UIRESUME)
% uiwait(handles.trial_by_trial_gui);


% --- Outputs from this function are returned to the command line.
function varargout = trial_by_trial_gui_OutputFcn(hObject, eventdata, handles) 
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
% hObject    handle to forward (see GCBO)Select Address
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

current_xlim = get(handles.data_axes,'xlim');
current_ylim = get(handles.data_axes,'ylim');

trace_ind = str2num(get(handles.trial_number,'String'));

start_ind = 1; %2450
trace_len = length(handles.data.sweeps{trace_ind}(:,1));

gain_mult = 1;

% axes(handles.data_axes)
% plot((0:trace_len-1)/20000,handles.data.sweeps{trace_ind}(:,1)*gain_mult,'b');
% hold on;
% plot((0:trace_len-1)/20000,handles.data.sweeps{trace_ind}(:,2)*gain_mult,'r');


timebase = (0:trace_len-1)/20000;
if isfield(handles.data,'trial_metadata')
    
    if strcmp(handles.data.trial_metadata(trace_ind).stim_type,'LED')
        stim_sweep = handles.data.sweeps{trace_ind}(:,2);
    elseif strcmp(handles.data.trial_metadata(trace_ind).stim_type,'2P')
        
        stim_sweep = handles.data.sweeps{trace_ind}(:,2);
        if isfield(handles.data.trial_metadata,'lut_used') && handles.data.trial_metadata(trace_ind).lut_used
            stim_sweep = stim_sweep/max(stim_sweep)*handles.data.trial_metadata(trace_ind).pulseamp;
        end
    end
else
    
    stim_sweep = handles.data.sweeps{trace_ind}(:,2)*gain_mult;
    
end
% handles.data_axes = plotyy(handles.data_axes(1),timebase,stim_sweep,timebase,handles.data.sweeps{trace_ind}(:,1)*gain_mult);
axes(handles.data_axes);
if strcmp('current-clamp',handles.data.trial_metadata(trace_ind).cell1_clamp_type)
        gain_mult = 1/20;
end

switch handles.data.trial_metadata(trace_ind).cell1_clamp_type
    case 'current-clamp'
        this_trace = handles.data.sweeps{trace_ind}(:,1)*1;
    case 'voltage-clamp'
        this_trace = handles.data.sweeps{trace_ind}(:,1);
    case 'cell-attached'
        this_trace = -(handles.data.sweeps{trace_ind}(:,1) - median(handles.data.sweeps{trace_ind}(1:10,1)));
end
if get(handles.hpf,'Value')
    this_trace = highpass_filter(this_trace,20000);
end

plot(timebase,this_trace); hold on; plot(timebase,handles.data.sweeps{trace_ind}(:,3)/max(handles.data.sweeps{trace_ind}(:,3))*10 + 10)
hold on;
if get(handles.draw_thresh,'Value')
    thresh = str2double(get(handles.thresh,'String'));
    plot(timebase,thresh*ones(size(timebase)))
    hold on;
    crossings = detect_peaks(this_trace',thresh,30,0,length(this_trace),-Inf,0,0,0);
    crossings = crossings{1};
    scatter(timebase(crossings),thresh*ones(size(crossings)))
    hold on
end

% hold on;
% scatter(handles.data.time(find(handles.data.spikes(trace_ind,:))),100*ones(1,sum(handles.data.spikes(trace_ind,:))),20*ones(1,sum(handles.data.spikes(trace_ind,:))),'filled');
% hold on;
% scatter(handles.data.time(find(handles.data.stims(trace_ind,:))),120*ones(1,sum(handles.data.stims(trace_ind,:))),20*ones(1,sum(handles.data.stims(trace_ind,:))),'filled');
hold off;
if ~get(handles.hold_axes,'Value')
    axis tight
else
    set(handles.data_axes,'xlim',current_xlim)
    set(handles.data_axes,'ylim',current_ylim)
end
% if trace_ind ~= 1
% if isfield(handles.data,'trialtime')
%     title(['Experiment Time: ' num2str(handles.data.trialtime(trace_ind)) ' sec, ' num2str(diff(handles.data.trialtime(max(trace_ind-1,1):trace_ind))) 'sec since prev trial'])
title_string = '';
if isfield(handles.data,'trial_metadata') && isfield(handles.data.trial_metadata(trace_ind),'run_count')
    title_string = num2str(handles.data.trial_metadata(trace_ind).run_count);
end
if isfield(handles.data,'trial_metadata') && isfield(handles.data.trial_metadata(trace_ind),'relative_position')
    title_string = [title_string ': ' mat2str(handles.data.trial_metadata(trace_ind).relative_position)];
end
if isfield(handles.data,'trialtime')
    title_string = [title_string ': Experiment Time: ' num2str(handles.data.trialtime(trace_ind)) ' sec, ' num2str(diff(handles.data.trialtime(max(trace_ind-1,1):trace_ind))) 'sec since prev trial'];
end
title(title_string)


    % else

%     title(['Experiment Time: ' num2str(handles.data.trialtime(trace_ind)) ' sec'])
% end
% ylim([-700 200])


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


% --- Executes on button press in hold_axes.
function hold_axes_Callback(hObject, eventdata, handles)
% hObject    handle to hold_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold_axes

if ~get(hObject,'Value')
    axes(handles.data_axes)
    axis tight
end


% --- Executes on button press in hpf.
function hpf_Callback(hObject, eventdata, handles)
% hObject    handle to hpf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hpf

draw_plot(handles)



function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double

draw_plot(handles)


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in draw_thresh.
function draw_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to draw_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of draw_thresh

draw_plot(handles)

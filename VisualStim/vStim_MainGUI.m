function varargout = vStim_MainGUI(varargin)
% VSTIM_MAINGUI MATLAB code for vStim_MainGUI.fig
%      VSTIM_MAINGUI, by itself, creates a new VSTIM_MAINGUI or raises the existing
%      singleton*.
%
%      H = VSTIM_MAINGUI returns the handle to a new VSTIM_MAINGUI or the handle to
%      the existing singleton*.
%
%      VSTIM_MAINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VSTIM_MAINGUI.M with the given input arguments.
%
%      VSTIM_MAINGUI('Property','Value',...) creates a new VSTIM_MAINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vStim_MainGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes propertyOPpor application
%      stop.  All inputs are passed to vStim_MainGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vStim_MainGUI

% Last Modified by GUIDE v2.5 06-Feb-2017 16:19:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @vStim_MainGUI_OpeningFcn, ...
    'gui_OutputFcn',  @vStim_MainGUI_OutputFcn, ...
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


% --- Executes just before vStim_MainGUI is made visible.
function vStim_MainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vStim_MainGUI (see VARARGIN)

% Choose default command line output for vStim_MainGUI
handles.output = hObject;

% make sure toggle buttons are in the correct state
handles.Status.String = 'No paradigm loaded';
handles.LoadParadigm.String = 'Load selected paradigm';handles.LoadParadigm.Value = false;
handles.PreviewParadigm.Enable = 'off';handles.PreviewParadigm.Value = false;
handles.StartParadigm.Enable = 'off';handles.StartParadigm.Value = false;

%% check for existing paradigms
cPath = fileparts(which(handles.figure1.Name)); %path to imager - should contain paradigm files
Paradigms = dir([cPath filesep 'Settings' filesep 'Default' filesep '*.txt']);
addpath([cPath filesep 'Paradigms']); %make sure paradigms are in the search path
Cnt = 1;

for iParams = 1:length(Paradigms)
    [~,temp{iParams}] = fileparts(Paradigms(iParams).name);
    if ~isempty(dir([cPath filesep 'Paradigms' filesep temp{iParams} '.m'])) %check for .m file
        pNames{Cnt} = temp{iParams};
        Cnt = Cnt +1;
    end
end

handles.SelectParadigm.String = pNames;
handles.SelectParadigm.Value = 1;
clear pNames

SelectParadigm_Callback(handles.SelectParadigm,[],handles) % load variables for current paradigm
StaticVariableNames_Callback(handles.StaticVariableNames, eventdata, handles) % load current static variable
FlexibleVariableNames_Callback(handles.FlexibleVariableNames, eventdata, handles) % load current flexible variable
cFlexibleVariable_Callback(handles.cFlexibleVariable, [], handles) %call cFlexible callback to update the sum of possible cases

%% initialize serial port
IOPort('CloseAll');
handles.SerialPort = [];
handles.Arduino = [];
handles.WavePlayer = [];
info.AvailableSerialPorts = FindSerialPort();
handles.SerialDevices.Value = 1;
handles.trialByte = 102;
handles.stimByte = 101;

if isempty(info.AvailableSerialPorts)
    disp('No serial port found')
    handles.SerialDevices.String = 'none'; %serial device selector
    handles.DAQdeviceName.String = 'No port detected'; %show device name
else
    Ports = info.AvailableSerialPorts; %find all serial ports
    
    for x = 1 : length(Ports)
        try
            % search analog output module
            handles.WavePlayer = BpodWavePlayer(Ports{x});
            fprintf('Found analog output module on port: %s\n', Ports{x})
            Ports(x) = []; %remove port from the list if this worked
        end
    end
    if isempty(handles.WavePlayer)
        disp('Found no analog output module');
    end

    [handles,Found] = CheckArduinoPort(Ports,handles); %check if correct serial device has the right handshake to identify as stim arduino
    if Found ~= 0
        disp(['Found responsive arduino on port ' Ports{Found}]);
        flushinput(handles.Arduino);
        info.AvailableSerialPorts(ismember(info.AvailableSerialPorts,handles.Arduino.Port)) = [];
    else
        disp('No responsive arduino found.');
        handles.Arduino = [];
    end
    
    handles.SerialDevices.String = Ports; %show available ports
    if isempty(handles.Arduino) && isempty(handles.WavePlayer)
        handles.DAQdeviceName.String = 'Arduino: none - WavePlayer: none'; %show port name
    elseif ~isempty(handles.Arduino) && isempty(handles.WavePlayer)
        handles.DAQdeviceName.String = ['Arduino: ' handles.Arduino.Port ' - WavePlayer: none']; %show port name
    elseif isempty(handles.Arduino) && ~isempty(handles.WavePlayer)
        handles.DAQdeviceName.String = ['Arduino: none - WavePlayer: ' handles.WavePlayer.Port.PortName]; %show port name
    elseif ~isempty(handles.Arduino) && ~isempty(handles.WavePlayer)
        handles.DAQdeviceName.String = ['Arduino: ' handles.Arduino.Port ' - WavePlayer: ' handles.WavePlayer.Port.PortName]; %show port name
    end
end

%% Check data path
if ~isdir(handles.ShowPath.String)
    mkdir(handles.ShowPath.String);
end
handles.DataPath = handles.ShowPath.String;

%% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vStim_MainGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vStim_MainGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in SelectParadigm.
function SelectParadigm_Callback(hObject, eventdata, handles)
% hObject    handle to SelectParadigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectParadigm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectParadigm

cPath = fileparts(which(handles.figure1.Name)); %path to imager - should contain paradigm files
cParam = handles.SelectParadigm.String{handles.SelectParadigm.Value}; %current selection
fID = fopen([cPath filesep 'Settings' filesep 'Default' filesep cParam '.txt']); %open param text file to load containing variables
Cnt = 0;
pVals =blanks(100); %allocate variable array

% get contents from text file into pre-filled char array
while ~isempty(Cnt)
    pOut = fgets(fID);
    if pOut ~= -1
        Cnt = Cnt+1;
        pOut = strtrim(pOut);
        pNames(Cnt,:) = blanks(100); %allocate string array
        pVals(Cnt,:) = blanks(100); %allocate string array
        pNames(Cnt,1:strfind(pOut,' ')-1) = pOut(1:strfind(pOut,' ')-1);
        pVals(Cnt,1:length(strtrim(pOut(strfind(pOut,' ')+3:end)))) = strtrim(pOut(strfind(pOut,' ')+3:end));
    else
        Cnt = [];
    end
end
fclose(fID);

%trim to longest variable name/value
pNames(:,106-min(sum(ismember(pNames,' '),2)):end) = [];
pVals = strtrim(pVals);

%find flexible vars and show in according listboxes
ind = ismember(pNames,'f');
handles.StaticVariableNames.String = cellstr([pNames(~ind(:,1),2:size(pNames,2)) pVals(~ind(:,1),:)]);
if isempty(handles.StaticVariableNames.String{1})
    handles.StaticVariableNames.Value = 0;
    handles.StaticVariableNames.String = [];
else
    handles.StaticVariableNames.Value = 1;
    StaticVariableNames_Callback(handles.StaticVariableNames, eventdata, handles)
    cStaticVariable_Callback(handles.cStaticVariable, [], handles)
end

handles.FlexibleVariableNames.String = cellstr([pNames(ind(:,1),2:size(pNames,2)) pVals(ind(:,1),:)]);
if isempty(handles.FlexibleVariableNames.String{1})
    handles.FlexibleVariableNames.Value = 0;
    handles.FlexibleVariableNames.String = [];
else
    handles.FlexibleVariableNames.Value = 1;
    FlexibleVariableNames_Callback(handles.FlexibleVariableNames, eventdata, handles)
end
cFlexibleVariable_Callback(handles.cFlexibleVariable, [], handles)


% --- Executes during object creation, after setting all properties.
function SelectParadigm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectParadigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StaticVariableNames.
function StaticVariableNames_Callback(hObject, eventdata, handles)
% hObject    handle to StaticVariableNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StaticVariableNames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StaticVariableNames

if hObject.Value > 0
    cVal = textscan(hObject.String{hObject.Value},'%s');
    cVal{2} = str2num(hObject.String{hObject.Value}(length(cVal{1}{1})+1:end)); %values for variable
    
    if length(cVal{2}) > 1
        handles.cStaticVariable.String = [cVal{1}{1} ' = [' num2str(cVal{2}) ']'];
    else
        handles.cStaticVariable.String = [cVal{1}{1} ' = ' num2str(cVal{2}) ''];
    end
end

% --- Executes during object creation, after setting all properties.
function StaticVariableNames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StaticVariableNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FlexibleVariableNames.
function FlexibleVariableNames_Callback(hObject, eventdata, handles)
% hObject    handle to FlexibleVariableNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FlexibleVariableNames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FlexibleVariableNames

if hObject.Value > 0
    cVal = textscan(hObject.String{hObject.Value},'%s');
    cVal{2} = str2num(hObject.String{hObject.Value}(length(cVal{1}{1})+1:end)); %values for variable
    
    if length(cVal{2}) > 1
        handles.cFlexibleVariable.String = [cVal{1}{1} ' = [' num2str(cVal{2}) ']'];
    else
        handles.cFlexibleVariable.String = [cVal{1}{1} ' = ' num2str(cVal{2}) ''];
    end
end

% --- Executes during object creation, after setting all properties.
function FlexibleVariableNames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FlexibleVariableNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MoveToFlex.
function MoveToFlex_Callback(hObject, eventdata, handles)
% hObject    handle to MoveToFlex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MoveToFlex

if ~(size(handles.StaticVariableNames.String,1) == 0)
    cValue = handles.StaticVariableNames.String{handles.StaticVariableNames.Value};
    handles.StaticVariableNames.String(handles.StaticVariableNames.Value) = [];
    handles.FlexibleVariableNames.String{end+1} = cValue;
    
    if handles.StaticVariableNames.Value > size(handles.StaticVariableNames.String,1) %make sure value of listbox is in cell range
        handles.StaticVariableNames.Value = size(handles.StaticVariableNames.String,1);
    end
    if handles.FlexibleVariableNames.Value == 0 %make sure value of listbox is in cell range
        handles.FlexibleVariableNames.Value = 1;
    end
end
%call cFlexible callback to update the sum of possible cases
FlexibleVariableNames_Callback(handles.FlexibleVariableNames, [], handles)
cFlexibleVariable_Callback(handles.cFlexibleVariable, [], handles)

StaticVariableNames_Callback(handles.StaticVariableNames, [], handles)
cStaticVariable_Callback(handles.cStaticVariable, [], handles)



% --- Executes on button press in MoveToStatic.
function MoveToStatic_Callback(hObject, eventdata, handles)
% hObject    handle to MoveToStatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~(size(handles.FlexibleVariableNames.String,1) == 0)
    cValue = handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value};
    handles.FlexibleVariableNames.String(handles.FlexibleVariableNames.Value) = [];
    cVal = textscan(cValue,'%s');
    cVal{2} = str2num(cValue(length(cVal{1}{1})+1:end)); %values for variable
    if length(cVal{2})>1 %if flexible variable has multiple values, only use the first one for static
        cValue = strrep(cValue,['[' num2str(cVal{2}) ']'],num2str(cVal{2}(1)));
    end
    handles.StaticVariableNames.String{end+1} = cValue;
    
    if handles.FlexibleVariableNames.Value > size(handles.FlexibleVariableNames.String,1)
        handles.FlexibleVariableNames.Value = size(handles.FlexibleVariableNames.String,1);
    end
    if handles.StaticVariableNames.Value == 0 %make sure value of listbox is in cell range
        handles.StaticVariableNames.Value = 1;
    end
end
%call cFlexible callback to update the sum of possible cases
FlexibleVariableNames_Callback(handles.FlexibleVariableNames, [], handles)
cFlexibleVariable_Callback(handles.cFlexibleVariable, [], handles)

StaticVariableNames_Callback(handles.StaticVariableNames, [], handles)
cStaticVariable_Callback(handles.cStaticVariable, [], handles)

function cStaticVariable_Callback(hObject, eventdata, handles)
% hObject    handle to Tag_CurrentStaticVariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tag_CurrentStaticVariable as text
%        str2double(get(hObject,'String')) returns contents of Tag_CurrentStaticVariable as a double

if ~(size(handles.StaticVariableNames.String,1) == 0)
    cVal = textscan(handles.StaticVariableNames.String{handles.StaticVariableNames.Value},'%s'); %string from listbox
    cVal{2} = str2num(handles.StaticVariableNames.String{handles.StaticVariableNames.Value}(length(cVal{1}{1})+1:end)); %values from listbox
    cSelect = str2num(hObject.String(length([cVal{1}{1} ' = '])+1:end)); %current value
    
    if ~isempty(cSelect)
        cSelect = cSelect(1); %make sure variable has only a single value
        if length(cVal{2}) > 1
            handles.StaticVariableNames.String{handles.StaticVariableNames.Value} = strrep(handles.StaticVariableNames.String{handles.StaticVariableNames.Value},['[' num2str(cVal{2}) ']'],num2str(cSelect));
        else
            handles.StaticVariableNames.String{handles.StaticVariableNames.Value} = strrep(handles.StaticVariableNames.String{handles.StaticVariableNames.Value},num2str(cVal{2}),num2str(cSelect));
        end
        hObject.String = [cVal{1}{1} ' = ' num2str(cSelect)];
    else
        if length(cVal{2}) > 1
            hObject.String = [cVal{1}{1} ' = [' num2str(cVal{2}) ']'];
        else
            hObject.String = [cVal{1}{1} ' = ' num2str(cVal{2})];
        end
    end
else
    hObject.String = 'No variable loaded';
end

function cFlexibleVariable_Callback(hObject, eventdata, handles)
% hObject    handle to Tag_CurrentFlexibleVariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tag_CurrentFlexibleVariable as text
%        str2double(get(hObject,'String')) returns contents of Tag_CurrentFlexibleVariable as a double

if ~(size(handles.FlexibleVariableNames.String,1) == 0)
    cVal = textscan(handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value},'%s'); %string from listbox
    cVal{2} = str2num(handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value}(length(cVal{1}{1})+1:end)); %values from listbox
    cSelect = str2num(hObject.String(length([cVal{1}{1} ' = '])+1:end)); %current value
    
    if ~isempty(cSelect)
        if length(cSelect) > 1
            hObject.String = [cVal{1}{1} ' = [' num2str(cSelect) ']'];
            if length(cVal{2}) > 1
                handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value} = strrep(handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value},num2str(cVal{2}),num2str(cSelect));
            else
                handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value} = strrep(handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value},num2str(cVal{2}),['[' num2str(cSelect) ']']);
            end
        else
            if length(cVal{2}) > 1
                handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value} = strrep(handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value},['[' num2str(cVal{2}) ']'],num2str(cSelect));
            else
                handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value} = strrep(handles.FlexibleVariableNames.String{handles.FlexibleVariableNames.Value},num2str(cVal{2}),num2str(cSelect));
            end
            hObject.String = [cVal{1}{1} ' = ' num2str(cSelect)];
        end
    else
        if length(cVal{2}) > 1
            hObject.String = [cVal{1}{1} ' = [' num2str(cVal{2}) ']'];
        else
            hObject.String = [cVal{1}{1} ' = ' num2str(cVal{2})];
        end
    end
    % Compute number of possible cases, based on combinations of flexible variables
    for iEntries = 1:size(handles.FlexibleVariableNames.String,1)
        cVal = textscan(handles.FlexibleVariableNames.String{iEntries},'%s'); %string from listbox
        FlexVals{iEntries} = str2num(handles.FlexibleVariableNames.String{iEntries}(length(cVal{1}{1})+1:end)); %value from listbox
    end
    handles.CaseSum.String = num2str(size(CombVec(FlexVals{:}),2)); %get number of possible combinations
else
    hObject.String = 'No variable loaded';
    handles.CaseSum.String = '1';
end

% --- Executes during object creation, after setting all properties.
function Tag_CurrentFlexibleVariable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tag_CurrentFlexibleVariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScreenSize_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ScreenSize = textscan(hObject.String,'%f%s%f'); %ScreenSize in cm
ScreenDist = str2double(handles.EyeDistance.String); %ScreenDistance in cm
ScreenSize = 2*atan([ScreenSize{1} ScreenSize{3}]/(2*ScreenDist))./pi*180; %Screen size in visual degrees
ScreenSize = round(ScreenSize,2);
handles.ScreenSizeAngle.String = [num2str(ScreenSize(1)) ' x ' num2str(ScreenSize(2))];


% --- Executes during object creation, after setting all properties.
function ScreenSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScreenSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScreenRes_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScreenRes as text
%        str2double(get(hObject,'String')) returns contents of ScreenRes as a double


% --- Executes during object creation, after setting all properties.
function ScreenRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScreenRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EyeDistance_Callback(hObject, eventdata, handles)
% hObject    handle to EyeDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EyeDistance as text
%        str2double(get(hObject,'String')) returns contents of EyeDistance as a double

ScreenSize = textscan(handles.ScreenSize.String,'%f%s%f'); %ScreenSize in cm
ScreenDist = str2double(handles.EyeDistance.String); %ScreenDistance in cm
ScreenSize = 2*atan([ScreenSize{1} ScreenSize{3}]/(2*ScreenDist))./pi*180; %Screen size in visual degrees
ScreenSize = round(ScreenSize,2);
handles.ScreenSizeAngle.String = [num2str(ScreenSize(1)) ' x ' num2str(ScreenSize(2))];
ScreenAngle_Callback(handles.ScreenAngle, [], handles);

% --- Executes during object creation, after setting all properties.
function EyeDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EyeDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SerialDevices.
function SerialDevices_Callback(hObject, eventdata, handles)
% hObject    handle to SerialDevices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.Arduino)
    close(handles.Arduino);
end

Ports = FindSerialPort;
Ports(ismember(Ports,handles.WavePlayer.Port.PortName)) = []; %ignore waveplayer port

[handles,Found] = CheckArduinoPort(Ports,handles); %check if correct serial device has the right handshake to identify as stim arduino
if Found ~= 0
    disp(['Found responsive arduino on port ' Ports{Found}]);
    flushinput(handles.Arduino);
else
    disp('No responsive arduino found.');
    handles.Arduino = [];
end

if isempty(handles.Arduino) && isempty(handles.WavePlayer)
    handles.DAQdeviceName.String = 'Arduino: none - WavePlayer: none'; %show port name
elseif ~isempty(handles.Arduino) && isempty(handles.WavePlayer)
    handles.DAQdeviceName.String = ['Arduino: ' handles.Arduino.Port ' - WavePlayer: none']; %show port name
elseif isempty(handles.Arduino) && ~isempty(handles.WavePlayer)
    handles.DAQdeviceName.String = ['Arduino: none - WavePlayer: ' handles.WavePlayer.Port.PortName]; %show port name
elseif ~isempty(handles.Arduino) && ~isempty(handles.WavePlayer)
    handles.DAQdeviceName.String = ['Arduino: ' handles.Arduino.Port ' - WavePlayer: ' handles.WavePlayer.Port.PortName]; %show port name
end

% Hints: contents = cellstr(get(hObject,'String')) returns SerialDevices contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SerialDevices


% --- Executes during object creation, after setting all properties.
function SerialDevices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SerialDevices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DAQdeviceName_Callback(hObject, eventdata, handles)
% hObject    handle to DAQdeviceName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DAQdeviceName as text
%        str2double(get(hObject,'String')) returns contents of DAQdeviceName as a double


% --- Executes during object creation, after setting all properties.
function DAQdeviceName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DAQdeviceName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NrTrials_Callback(hObject, eventdata, handles)
% hObject    handle to NrTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NrTrials as text
%        str2double(get(hObject,'String')) returns contents of NrTrials as a double


% --- Executes during object creation, after setting all properties.
function NrTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NrTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RandTrials.
function RandTrials_Callback(hObject, eventdata, handles)
% hObject    handle to RandTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to NrTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NrTrials as text
%        str2double(get(hObject,'String')) returns contents of NrTrials as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NrTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RandomTrials.
function RandomTrials_Callback(hObject, eventdata, handles)
% hObject    handle to RandomTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadParadigm.
function LoadParadigm_Callback(hObject, eventdata, handles)
% hObject    handle to LoadParadigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject.Value
    % Wait for user confirmation and run paradigm
    uiwait(msgbox('Press enter or space to continue','Wait','modal'));
    pause(1);
    hObject.String = 'Release current paradigm';
    eval([handles.SelectParadigm.String{handles.SelectParadigm.Value} '(handles);'],'sca; disp(''Stimulation error - aborted !'')')
else
    hObject.String = 'Load selected paradigm';
    handles.Status.String = 'No paradigm loaded';
    handles.StartParadigm.Enable = 'off';
    handles.StartParadigm.Value = false;
    handles.PreviewParadigm.Enable = 'off';
    handles.PreviewParadigm.Value = false;
    handles.StartParadigm.String = 'Start paradigm';
    handles.PreviewParadigm.String = 'Preview paradigm';
end


% --- Executes on button press in PreviewParadigm.
function PreviewParadigm_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewParadigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PreviewParadigm


% --- Executes on button press in StartParadigm.
function StartParadigm_Callback(hObject, eventdata, handles)
% hObject    handle to StartParadigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StartParadigm



function ITI_Callback(hObject, eventdata, handles)
% hObject    handle to ITI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ITI as text
%        str2double(get(hObject,'String')) returns contents of ITI as a double


% --- Executes during object creation, after setting all properties.
function ITI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ITI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CaseSum_Callback(hObject, eventdata, handles)
% hObject    handle to CaseSum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CaseSum as text
%        str2double(get(hObject,'String')) returns contents of CaseSum as a double


% --- Executes during object creation, after setting all properties.
function CaseSum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CaseSum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cFlexibleVariable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cFlexibleVariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadSettings.
function LoadSettings_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Custom' filesep handles.SelectParadigm.String{handles.SelectParadigm.Value}]; %path to custom settings of a given paradigm
if ~isdir(cPath)
    mkdir(cPath)
end

settingsFiles = dir(cPath); %find setting files
if size(ls(cPath)',1) > 2
    ind = listdlg('PromptString','Select a settings file',... %list dialog to select files
        'SelectionMode','single', 'ListString',{settingsFiles(3:end).name});
    if isempty(ind)
        cParam = 0; %no settings file selected
    else
        cParam = settingsFiles(ind+2).name;
    end
else
    cParam = 0; %no file found
end

if cParam == 0
    warning('No settings file to load. Save settings for current paradigm first.')
else
    fID = fopen([cPath filesep cParam]); %open param text file to load containing variables
    Cnt = 0;
    
    % get contents from text file into pre-filled char array
    while ~isempty(Cnt)
        pOut = fgets(fID);
        if pOut ~= -1
            Cnt = Cnt+1;
            pOut = strtrim(pOut);
            pNames(Cnt,:) = blanks(100); %allocate string array
            pVals(Cnt,:) = blanks(100); %allocate string array
            pNames(Cnt,1:strfind(pOut,' ')-1) = pOut(1:strfind(pOut,' ')-1);
            pVals(Cnt,1:length(strtrim(pOut(strfind(pOut,' ')+3:end)))) = strtrim(pOut(strfind(pOut,' ')+3:end));
        else
            Cnt = [];
        end
    end
    fclose(fID);
    
    %trim to longest variable name/value
    pNames(:,106-min(sum(ismember(pNames,' '),2)):end) = [];
    pVals = strtrim(pVals);
    
    %find flexible vars and show in according listboxes
    ind = ismember(pNames,'f');
    handles.StaticVariableNames.String = cellstr([pNames(~ind(:,1),2:size(pNames,2)) pVals(~ind(:,1),:)]);
    handles.FlexibleVariableNames.String = cellstr([pNames(ind(:,1),2:size(pNames,2)) pVals(ind(:,1),:)]);
    if isempty(handles.FlexibleVariableNames.String{1})
        handles.FlexibleVariableNames.Value = 0;
        handles.FlexibleVariableNames.String = [];
    else
        handles.FlexibleVariableNames.Value = 1;
        FlexibleVariableNames_Callback(handles.FlexibleVariableNames, eventdata, handles)
    end
end

% --- Executes on button press in SaveToFile.
function SaveToFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveToFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sVars = handles.StaticVariableNames.String;
for iVars = 1:size(sVars,1)
    sVars{iVars} = ['s' sVars{iVars}]; %add 's' to identify as static var in the file
end

fVars = handles.FlexibleVariableNames.String;
for iVars = 1:size(fVars,1)
    fVars{iVars} = ['f' fVars{iVars}]; %add 'f' to identify as flexible var in the file
end

%place settings file in custom folder for current paradigm
cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Custom' filesep handles.SelectParadigm.String{handles.SelectParadigm.Value} filesep]; %path to custom settings of a given paradigm
if ~isdir(cPath)
    mkdir(cPath)
end

% get name for settings file and check to avoid overwrite
dPrompt = {'Enter the name of the settings file'};
pName = 'Save settings file';
cParam = inputdlg(dPrompt,pName,1,{'CustomSettings'});
if isempty(cParam)
    return;
end
cName = avoidOverwrite([cParam{1} '.txt'],cPath);

%save file
fID = fopen([cPath cName],'wt');
Vars = [sVars;fVars];
for iRows = 1:size(Vars,1)
    fprintf(fID,[Vars{iRows} char(10)],'%s');
end
fclose(fID);


% --- Executes on button press in SaveAsDefault.
function SaveAsDefault_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sVars = handles.StaticVariableNames.String;
for iVars = 1:size(sVars,1)
    sVars{iVars} = ['s' sVars{iVars}]; %add 's' to identify as static var in the file
end

fVars = handles.FlexibleVariableNames.String;
for iVars = 1:size(fVars,1)
    fVars{iVars} = ['f' fVars{iVars}]; %add 'f' to identify as flexible var in the file
end

cParam = handles.SelectParadigm.String{handles.SelectParadigm.Value}; %current selection
cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Default']; %path to default settings of a given paradigm

%open param text and replace settings with current selection
fID = fopen([cPath filesep cParam '.txt'],'wt');
Vars = [sVars;fVars];
for iRows = 1:size(Vars,1)
    fprintf(fID,[Vars{iRows} char(10)],'%s');
end
fclose(fID);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sca; %clean up psychtoolbox windows
IOPort('CloseAll');
if ~isempty(handles.Arduino) && isobject(handles.Arduino)
    delete(handles.Arduino);
end
% savefig(handles.figure1,handles.figure1.FileName); %save changes in the figure back to .fig file

% Hint: delete(hObject) closes the figure
delete(hObject);


function SubjectName_Callback(hObject, eventdata, handles)
% hObject    handle to SubjectName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SubjectName as text
%        str2double(get(hObject,'String')) returns contents of SubjectName as a double


% --- Executes during object creation, after setting all properties.
function SubjectName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubjectName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ExperimentNr_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExperimentNr as text
%        str2double(get(hObject,'String')) returns contents of ExperimentNr as a double


% --- Executes during object creation, after setting all properties.
function ExperimentNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChangePath.
function ChangePath_Callback(hObject, eventdata, handles)
% hObject    handle to ChangePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.DataPath = uigetdir;
handles.ShowPath.String = handles.DataPath;
guidata(hObject, handles);

% --- Executes on button press in ShowPath.
function ShowPath_Callback(hObject, eventdata, handles)
% hObject    handle to ShowPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in TriggerMode.
function TriggerMode_Callback(hObject, eventdata, handles)
% hObject    handle to TriggerMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TriggerMode

if hObject.Value
    hObject.String = 'Triggermode: Slave';
else
    hObject.String = 'Triggermode: Master';
end

function ScreenSizeAngle_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenSizeAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScreenSizeAngle as text
%        str2double(get(hObject,'String')) returns contents of ScreenSizeAngle as a double


% --- Executes during object creation, after setting all properties.
function ScreenSizeAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScreenSizeAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunBatchSeq.
function RunBatchSeq_Callback(hObject, eventdata, handles)
% hObject    handle to RunBatchSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RunBatchSeq

if hObject.Value
    cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Batch' filesep ]; %path to batch files
    if ~isdir(cPath)
        mkdir(cPath)
    end
    
    settingsFiles = dir(cPath); %find setting files
    if size(ls(cPath)',1) > 2
        ind = listdlg('PromptString','Select a file',... %list dialog to select files
            'SelectionMode','single', 'ListString',{settingsFiles(3:end).name});
        if isempty(ind)
            cFile = 0; %no batch file selected
        else
            cFile = settingsFiles(ind+2).name;
        end
    else
        cFile = 0; %no batch file found
    end
    
    if cFile == 0
        warning('No batch files found. Create a batch sequence first.')
        hObject.Value = false;
        return;
    else
        %% Wait for user confirmation before moving on
        uiwait(msgbox('Press enter or space to continue','Wait','modal'));
        pause(1);
        
        % load batch file and run through sequence
        load([cPath cFile])
        try
            for iRuns = 1:length(cParadigm)
                
                handles.NrTrials.String = num2str(cParam(iRuns,1)); %update trialcount for current paradigm
                handles.ITI.String = num2str(cParam(iRuns,2)); %update ITI for current paradigm
                
                %% get new settings if a custom settings file should be loaded
                if ~isempty(cSettings{iRuns})
                    cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Custom' filesep cParadigm{iRuns} filesep]; %path to custom settings of a given paradigm
                    fID = fopen([cPath cSettings{iRuns}]); %open param text file to load containing variables
                    Cnt = 0;
                    % get contents from text file into pre-filled char array
                    while ~isempty(Cnt)
                        pOut = fgets(fID);
                        if pOut ~= -1
                            Cnt = Cnt+1;
                            pOut = strtrim(pOut);
                            pNames(Cnt,:) = blanks(100); %allocate string array
                            pVals(Cnt,:) = blanks(100); %allocate string array
                            pNames(Cnt,1:strfind(pOut,' ')-1) = pOut(1:strfind(pOut,' ')-1);
                            pVals(Cnt,1:length(strtrim(pOut(strfind(pOut,' ')+3:end)))) = strtrim(pOut(strfind(pOut,' ')+3:end));
                        else
                            Cnt = [];
                        end
                    end
                    fclose(fID);
                    
                    %trim to longest variable name/value
                    pNames(:,106-min(sum(ismember(pNames,' '),2)):end) = [];
                    pVals = strtrim(pVals);
                    
                    %find flexible vars and show in according listboxes
                    ind = ismember(pNames,'f');
                    handles.StaticVariableNames.String = cellstr([pNames(~ind(:,1),2:size(pNames,2)) pVals(~ind(:,1),:)]);
                    handles.FlexibleVariableNames.String = cellstr([pNames(ind(:,1),2:size(pNames,2)) pVals(ind(:,1),:)]);
                    if isempty(handles.FlexibleVariableNames.String{1})
                        handles.FlexibleVariableNames.Value = 0;
                        handles.FlexibleVariableNames.String = [];
                    else
                        handles.FlexibleVariableNames.Value = 1;
                        FlexibleVariableNames_Callback(handles.FlexibleVariableNames, eventdata, handles)
                    end
                    clear pNames pVals
                end
                
                %% launch current paradigm
                hObject.String = 'Release current sequence';
                disp(['Running paradigm ' int2str(iRuns) '/' int2str(length(cParadigm))])
                eval([cParadigm{iRuns} '(handles);']) %run paradigm
                pause(cParam(iRuns,3)); %take a pause after paradigm ran through
                
            end
        catch
            disp('Stimulation error - aborted !')
        end
    end
    hObject.Value = false;
    hObject.String = 'Load batch sequence';
    
else
    hObject.String = 'Load batch sequence';
    handles.Status.String = 'No paradigm loaded';
    handles.StartParadigm.Enable = 'off';
    handles.StartParadigm.Value = false;
    handles.PreviewParadigm.Enable = 'off';
    handles.PreviewParadigm.Value = false;
    handles.StartParadigm.String = 'Start paradigm';
    handles.PreviewParadigm.String = 'Preview paradigm';
end


% --- Executes on button press in MakeBatchSeq.
function MakeBatchSeq_Callback(hObject, eventdata, handles)
% hObject    handle to MakeBatchSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

checker = true;
Cnt = 1; %counter for paradigms to run
while checker
    
    ind = listdlg('PromptString','Select a paradigm to add',... %list dialog to select files
        'SelectionMode','single', 'ListString',handles.SelectParadigm.String);
    cParadigm{Cnt} = handles.SelectParadigm.String{ind};
    
    cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Custom' filesep cParadigm{Cnt}]; %path to custom settings of a given paradigm
    if ~isdir(cPath)
        mkdir(cPath)
    end
    
    settingsFiles = dir(cPath); %find setting files
    if size(ls(cPath)',1) > 2
        ind = listdlg('PromptString','Select a settings file',... %list dialog to select files
            'SelectionMode','single', 'ListString',{settingsFiles(3:end).name});
        cSettings{Cnt} = settingsFiles(ind+2).name;
    else
        cSettings{Cnt} = []; %no file found
        disp('Selected paradigm has no custom setting file, using default instead.')
    end
    
    dPrompt = {'Enter the number of trials'};
    pName = 'Trialcount';
    temp = inputdlg(dPrompt,pName,1,{'60'});
    cParam(Cnt,1) = str2num(temp{1});
    
    dPrompt = {'Enter the ITI in seconds'};
    pName = 'ITI';
    temp = inputdlg(dPrompt,pName,1,{'5'});
    cParam(Cnt,2) = str2num(temp{1});
    
    dPrompt = {'Enter pause after paradigm in seconds'};
    pName = 'Pause after paradigm';
    temp = inputdlg(dPrompt,pName,1,{'5'});
    cParam(Cnt,3) = str2num(temp{1});
    
    % check if another experiment should be added
    cResponse = questdlg('Add another experiment to batch?','Proceed?','Yes','No','No');
    if strcmp(cResponse,'No')
        checker = false;
    else
        Cnt = Cnt + 1;
    end
end

cPath = [fileparts(which(handles.figure1.Name)) filesep 'Settings' filesep 'Batch' filesep]; %path to batch file
dPrompt = {'Enter name of batch file'};
pName = 'Batch name';
cName = inputdlg(dPrompt,pName,1,{'DefaultBatch'});
cName = avoidOverwrite([cName{1} '.mat'],cPath); %check if file exists already
save([cPath cName],'cParadigm','cSettings','cParam')


% --- Executes on button press in ArduinoControl.
function ArduinoControl_Callback(hObject, eventdata, handles)
% hObject    handle to ArduinoControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.Arduino,handles.ArduinoSettings] = vStim_ArduinoControl(handles.Arduino); %get arduino handle and settings
guidata(hObject, handles); % update handles structure

%% additional code
function CandidatePorts = FindSerialPort
%code to find serial ports on either windows or linux systems
if ispc
    [~,RawString]=system('mode');
    CandidatePorts=regexp(RawString,'COM\d+','match');
elseif isunix
    %     RawSerialPortList = dir('/dev/ttyUSB*');
    %     RawSerialPortList = [ dir('/dev/ttyACM*'); RawSerialPortList];
    %     RawSerialPortList = [ dir('/dev/ttyS*'); RawSerialPortList];
    [~,RawSerialPortList] = system('ls /dev/ttyUSB*');
    [~,PortList] = system('ls /dev/ttyACM*');
    RawSerialPortList = [RawSerialPortList PortList];
    string = strtrim(RawSerialPortList);
    PortStringPositions = strfind(string, '/dev/tty');
    nPorts = length(PortStringPositions);
    CandidatePorts = cell(1,nPorts);
    for x = 1:nPorts-1
        CandidatePort = strtrim(string(PortStringPositions(x):PortStringPositions(x+1)-1));
        CandidatePorts{x} = CandidatePort;
    end
    CandidatePorts{end} = strtrim(string(PortStringPositions(end):end));
end

function [handles,Found] = CheckArduinoPort(Ports,handles)
% short code to find the correct port that is connected to the stimulator arduino.
Found = 0; x = 0;
while (Found == 0) && (x < length(Ports)) && ~isempty(Ports{1})
    x = x + 1;
    try
        disp(['Searching Arduino on port ' Ports{x}]);
        handles.Arduino = serial(Ports{x}, 'BaudRate', 115200, 'DataBits', 8, 'StopBits', 1, 'Timeout', 1, 'DataTerminalReady', 'off'); %create object for serial communication
        set(handles.Arduino, 'OutputBufferSize', 8000); %adjust buffer limits for serial communication
        set(handles.Arduino, 'InputBufferSize', 50000); %adjust buffer limits for serial communication
        if ~strcmp(handles.Arduino.Status,'open')
            fopen(handles.Arduino); % open serial port
            pause(1);
        end
        flushinput(handles.Arduino);
        fwrite(handles.Arduino,50); %check arduino
        
        tic
        while handles.Arduino.BytesAvailable == 0  % wait for handshake for a max of 1s
            if toc > 1
                break
            end
        end
        Byte = fread(handles.Arduino,1); %correct magic number is 10
        if Byte == 10
            Found = x; %correct arduino
        else
            fclose(handles.Arduino);
            delete(handles.Arduino);
            handles.Arduino = [];
        end
    catch
        delete(handles.Arduino);
        handles.Arduino = [];
    end
end

function outFile = avoidOverwrite(inFile,inPath,numDigits,startCount)
% Function to check if a given file already exists and create a new
% filename that can be used to avoid overwriting. New file names are
% created by appending a sequence like '_001' to the file name. The number
% is hereby increased until an unused file name is found.
%
% Usage: outFile = avoidOverwrite(inFile,inPath,numDigits)
%
%        inFile: The name of the file that should be checked. Make sure
%        to add filetype if required.
%        inPath: The path of the file. Can be omitted to only check in
%        the Matlab path.
%        numDigits: The number of digits used for enumration. Default is 2
%        digits, so filenames are set up as e.g. 'inFile_01.mat' and so on.
%        startCount: Defines the starting number for enumeration. Standard
%        is 0 so the first file is created as 'inFile_00.mat'.
%        outFile: The filename that can be used to save a file without
%        overwriting existing data.

%% check input
if ~exist('inPath','var')
    inPath = [];     %number of digits when adding numbers to a filename
end

if ~exist('numDigits','var')
    numDigits = 2;     %number of digits when adding numbers to a filename
end
numDigits = num2str(numDigits);

if ~exist('startCount','var')
    startCount = 0;     %first value for counter. This determines the first filename during enumeration.
end

if ~strcmp(inPath(end),filesep) %make sure path has a seperator at the end
    inPath = [inPath filesep];
end

%% check if file exists already and enumerate
if exist([inPath inFile],'file') == 2 %file exists already, check for alternative
    [~,shortFile,fileType] = fileparts(inFile); %exclude the file type
    checker = true; %check for alternate file names
    Cnt = startCount; %counter for file name
    
    while checker
        testPath = [inPath shortFile '_' num2str(Cnt, ['%0' numDigits 'i']) fileType];
        
        if exist(testPath,'file') == 2
            Cnt = Cnt + 1; %increase counter until a non-existing file name is found
        else
            checker = false;
        end
        
        if Cnt == 10^numDigits-1 && checker
            numDigits = numDigits+1;
            warning(['No unused file found at given number of digits. Number of digits increased to ' num2str(numDigits) '.']);
        end
    end
    outFile = [shortFile '_' num2str(Cnt, ['%0' numDigits 'i']) fileType];
    
else
    outFile = inFile;
end



function EyeLineDistBottom_Callback(hObject, eventdata, handles)
% hObject    handle to EyeLineDistBottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EyeLineDistBottom as text
%        str2double(get(hObject,'String')) returns contents of EyeLineDistBottom as a double


% --- Executes during object creation, after setting all properties.
function EyeLineDistBottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EyeLineDistBottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EyeLineDistTop_Callback(hObject, eventdata, handles)
% hObject    handle to EyeLineDistTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EyeLineDistTop as text
%        str2double(get(hObject,'String')) returns contents of EyeLineDistTop as a double


% --- Executes during object creation, after setting all properties.
function EyeLineDistTop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EyeLineDistTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EyeDistTop_Callback(hObject, eventdata, handles)
% hObject    handle to EyeDistTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EyeDistTop as text
%        str2double(get(hObject,'String')) returns contents of EyeDistTop as a double


% --- Executes during object creation, after setting all properties.
function EyeDistTop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EyeDistTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EyeDistLeft_Callback(hObject, eventdata, handles)
% hObject    handle to EyeDistLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EyeDistLeft as text
%        str2double(get(hObject,'String')) returns contents of EyeDistLeft as a double


% --- Executes during object creation, after setting all properties.
function EyeDistLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EyeDistLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScreenAngle_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

screenAngle = str2double(hObject.String);

if screenAngle > -90 || screenAngle < 90
    screenDistance = str2double(handles.EyeDistance.String);
    screenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
    screenHeight = screenSize{3};  % height of screen, in cm
    distanceTop = str2double(handles.EyeDistTop.String); %distance from screen top to animal fixation spot
    
    handles.EyeLineDistTop.String = num2str(screenDistance - (cos((90-screenAngle)/360*pi*2) * distanceTop));
    handles.EyeLineDistBottom.String = num2str(screenDistance + (cos((90-screenAngle)/360*pi*2) * (screenHeight-distanceTop)));
else
    hObject.string = '0';
    disp('Screen angle can only between -90 and 90 degrees');
end

% Hints: get(hObject,'String') returns contents of ScreenAngle as text
%        str2double(get(hObject,'String')) returns contents of ScreenAngle as a double


% --- Executes during object creation, after setting all properties.
function ScreenAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScreenAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveScreenSettings.
function SaveScreenSettings_Callback(hObject, eventdata, handles)
% hObject    handle to SaveScreenSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userResp = questdlg('This will save the current user settings into the figure ! Continue ?', ...
    'Save figure handles?', ...
    'Continue','Abort', 'Abort');

if strcmp(userResp,'Continue')
    userResp =  questdlg('Seriously - This can get the GUI in trouble if you dont know what you are doing. Really continue ?', ...
        'Save figure handles?', ...
        'Yes - I know what Im doing','Abort', 'Abort');
    
    if strcmp(userResp,'Yes - I know what Im doing')
        copyfile(handles.figure1.FileName,[handles.figure1.FileName(1:end-4) '_Backup.fig']); %make a backup of the original .fig file - just to be save
        savefig(handles.figure1,handles.figure1.FileName); %save changes in the figure back to .fig file
    end
end

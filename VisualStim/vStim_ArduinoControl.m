function varargout = vStim_ArduinoControl(varargin)
% VSTIM_ARDUINOCONTROL MATLAB code for vStim_ArduinoControl.fig
%      VSTIM_ARDUINOCONTROL, by itself, creates a new VSTIM_ARDUINOCONTROL or raises the existing
%      singleton*.
%
%      H = VSTIM_ARDUINOCONTROL returns the handle to a new VSTIM_ARDUINOCONTROL or the handle to
%      the existing singleton*.
%
%      VSTIM_ARDUINOCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VSTIM_ARDUINOCONTROL.M with the given input arguments.
%
%      VSTIM_ARDUINOCONTROL('Property','Value',...) creates a new VSTIM_ARDUINOCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vStim_ArduinoControl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vStim_ArduinoControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vStim_ArduinoControl

% Last Modified by GUIDE v2.5 04-Nov-2016 14:47:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vStim_ArduinoControl_OpeningFcn, ...
                   'gui_OutputFcn',  @vStim_ArduinoControl_OutputFcn, ...
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


% --- Executes just before vStim_ArduinoControl is made visible.
function vStim_ArduinoControl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vStim_ArduinoControl (see VARARGIN)

% Choose default command line output for vStim_ArduinoControl
handles.output = hObject;

% set up data table if currently empty
if isempty(handles.controlTable.Data{1})
    data = [15;1;1;20;60];
    if size(handles.controlTable.RowName,1) > length(data) %if more rows as pre-given values are present
        data(size(handles.controlTable.RowName,1)) = 0; %fill up with zeros
    end
    set(handles.controlTable, 'Data', num2cell(data));
    
    temp = handles.controlTable.Position(4);
    handles.controlTable.Position(3:4) = [handles.controlTable.Extent(3) handles.controlTable.Extent(4)];
    handles.controlTable.Position(2) = handles.controlTable.Position(2) + temp - (handles.controlTable.Extent(4));
end

%% check if serial object was provided
check = false;
if length(varargin)>=1 && ~isempty(varargin{1})
    if strcmp(varargin{1}.Type,'serial') && varargin{1}.isvalid
        handles.Arduino = varargin{1};
        flushinput(handles.Arduino);
        
        sCenter(2) = handles.controlTable.Data{1}; %get value for center of the stimulator
        sCenter(1) = floor(sCenter(2));sCenter(2) = (sCenter(2)-sCenter(1))*100; %modify for arduino transfer
        NegAmp(2) = handles.controlTable.Data{1} - handles.controlTable.Data{2}; %get negative value for amplitude of the stimulator
        NegAmp(1) = floor(NegAmp(2));NegAmp(2) = (NegAmp(2)-NegAmp(1))*100; %modify for arduino transfer
        PosAmp(2) = handles.controlTable.Data{1} + handles.controlTable.Data{2}; %get positive value for amplitude of the stimulator
        PosAmp(1) = floor(PosAmp(2));PosAmp(2) = (PosAmp(2)-PosAmp(1))*100; %modify for arduino transfer
        
        fwrite(handles.Arduino,[50 sCenter NegAmp PosAmp 1 num2str(1000/handles.controlTable.Data{5}) 'a' ...
            num2str(1000/handles.controlTable.Data{4}) 'a' num2str(handles.controlTable.Data{3}*1000)]); %send stim data to arduino
           
        tic
        while handles.Arduino.BytesAvailable == 0  % wait for handshake for a max of 1s
            if toc > 1
                break
            end
        end
        Byte = fread(handles.Arduino,1); %correct magic number is 10
        if Byte == 10
            check = true;
            handles.ArduinoDeviceName.String = ['Arduino port:  ' handles.Arduino.Port];
        else %wrong response. delete object and seach all ports
            fclose(handles.Arduino);
            delete(handles.Arduino);
            handles.Arduino = [];
        end
    end
end

if ~check %no valid serial object found, search for arduino yourself
    info.AvailableSerialPorts = FindSerialPort();
    Ports = info.AvailableSerialPorts; %find all serial ports
    [handles,Found] = CheckArduinoPort(Ports,handles); %check if correct serial device has the right handshake to identify as stim arduino
    if Found ~= 0
        disp(['Found responsive arduino on port ' Ports{Found}]);
        flushinput(handles.Arduino);
        info.AvailableSerialPorts(ismember(info.AvailableSerialPorts,handles.Arduino.Port)) = [];
        handles.ArduinoDeviceName.String = ['Arduino port:  ' Ports{Found}];
    else
        disp('No responsive arduino found.');
        handles.ArduinoDeviceName.String = 'No response Arduino found';
        handles.Arduino = [];
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vStim_ArduinoControl wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vStim_ArduinoControl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Arduino; %output arduino handle
varargout{2} = cell2mat(handles.controlTable.Data); %output arduino settings


function StimCenter_Callback(hObject, eventdata, handles)
% hObject    handle to StimCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimCenter as text
%        str2double(get(hObject,'String')) returns contents of StimCenter as a double
a

% --- Executes during object creation, after setting all properties.
function StimCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StimAmp_Callback(hObject, eventdata, handles)
% hObject    handle to StimAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimAmp as text
%        str2double(get(hObject,'String')) returns contents of StimAmp as a double


% --- Executes during object creation, after setting all properties.
function StimAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StimDuration_Callback(hObject, eventdata, handles)
% hObject    handle to StimDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimDuration as text
%        str2double(get(hObject,'String')) returns contents of StimDuration as a double


% --- Executes during object creation, after setting all properties.
function StimDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StimTest.
function StimTest_Callback(hObject, eventdata, handles)
% hObject    handle to StimTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.Arduino)
    flushinput(handles.Arduino);
    sCenter(2) = handles.controlTable.Data{1}; %get value for center of the stimulator
    sCenter(1) = floor(sCenter(2));sCenter(2) = (sCenter(2)-sCenter(1))*100; %modify for arduino transfer
    NegAmp(2) = handles.controlTable.Data{1} - handles.controlTable.Data{2}; %get negative value for amplitude of the stimulator
    NegAmp(1) = floor(NegAmp(2));NegAmp(2) = (NegAmp(2)-NegAmp(1))*100; %modify for arduino transfer
    PosAmp(2) = handles.controlTable.Data{1} + handles.controlTable.Data{2}; %get positive value for amplitude of the stimulator
    PosAmp(1) = floor(PosAmp(2));PosAmp(2) = (PosAmp(2)-PosAmp(1))*100; %modify for arduino transfer
    
    fwrite(handles.Arduino,[50 sCenter NegAmp PosAmp 1 num2str(handles.controlTable.Data{5}) 'a' ...
        num2str(1000/handles.controlTable.Data{4}) 'a' num2str(handles.controlTable.Data{3}*1000)]); %send stim data to arduino
                  
    tic
    while handles.Arduino.BytesAvailable == 0  % wait for handshake for a max of 5s
        if toc > 5
            break
        end
    end
    if handles.Arduino.BytesAvailable == 0
        disp('Arduino did not respond')
    else
        Byte = fread(handles.Arduino,1); %correct number is 10
        if Byte ~= 10
            disp(['Arduino sent wrong response: ' num2str(Byte) ' instead of 10'])
        end
    end
    
    Stim = get(handles.StimType,'value'); %get right stimulus to be tested.
    
    fwrite(handles.Arduino,Stim); %start stimulus sequence
    pause(0.2);
    Byte = fread(handles.Arduino,1); %correct number is 11
    if Byte(end) ~= (Stim + 10)
        disp(['Arduino sent wrong response: ' num2str(Byte) ' instead of ' num2str(Stim + 10)])
    else
        disp('Stimulation successful')
    end
end

% --- Executes on selection change in remSerialDevices.
function remSerialDevices_Callback(hObject, eventdata, handles)
% hObject    handle to remSerialDevices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns remSerialDevices contents as cell array
%        contents{get(hObject,'Value')} returns selected item from remSerialDevices


% --- Executes during object creation, after setting all properties.
function remSerialDevices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to remSerialDevices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ArduinoDeviceName_Callback(hObject, eventdata, handles)
% hObject    handle to ArduinoDeviceName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ArduinoDeviceName as text
%        str2double(get(hObject,'String')) returns contents of ArduinoDeviceName as a double


% --- Executes during object creation, after setting all properties.
function ArduinoDeviceName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ArduinoDeviceName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StimType.
function StimType_Callback(hObject, eventdata, handles)
% hObject    handle to StimType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StimType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StimType


% --- Executes during object creation, after setting all properties.
function StimType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in WhiskerSettings.
function WhiskerSettings_Callback(hObject, eventdata, handles)
% hObject    handle to WhiskerSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = [15;1;1;20;60]; % standard settings for whisker stimulation
if size(handles.controlTable.RowName,1) > length(data) %if more rows as pre-given values are present
    data(size(handles.controlTable.RowName,1)) = 0; %fill up with zeros
end
set(handles.controlTable, 'Data', num2cell(data));

% --- Executes on button press in TrunkSettings.
function TrunkSettings_Callback(hObject, eventdata, handles)
% hObject    handle to TrunkSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = [15;5;0.1;20;48]; % standard settings for whisker stimulation
if size(handles.controlTable.RowName,1) > length(data) %if more rows as pre-given values are present
    data(size(handles.controlTable.RowName,1)) = 0; %fill up with zeros
end
set(handles.controlTable, 'Data', num2cell(data));

% --- Executes on button press in WheelSettings.
function WheelSettings_Callback(hObject, eventdata, handles)
% hObject    handle to WheelSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.controlTable, 'Data');
data{3} = 1; data{4} = 1; data{5} = 60;
set(handles.controlTable, 'Data', data);



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
    [~,RawSerialPortList] = system('ls /dev/ttyACM*');
    [~,PortList] = system('ls /dev/ttyUSB*');
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
        handles.Arduino = serial(Ports{x}, 'BaudRate', 115200, 'DataBits', 8, 'StopBits', 1, 'Timeout', 1, 'DataTerminalReady', 'off'); %create object for serial communication
        set(handles.Arduino, 'OutputBufferSize', 8000); %adjust buffer limits for serial communication
        set(handles.Arduino, 'InputBufferSize', 50000); %adjust buffer limits for serial communication
        if ~strcmp(handles.Arduino.Status,'open')
            fopen(handles.Arduino); % open serial port
            pause(1);
        end
        flushinput(handles.Arduino);
        sCenter(2) = handles.controlTable.Data{1}; %get value for center of the stimulator
        sCenter(1) = floor(sCenter(2));sCenter(2) = (sCenter(2)-sCenter(1))*100; %modify for arduino transfer
        NegAmp(2) = handles.controlTable.Data{1} - handles.controlTable.Data{2}; %get negative value for amplitude of the stimulator
        NegAmp(1) = floor(NegAmp(2));NegAmp(2) = (NegAmp(2)-NegAmp(1))*100; %modify for arduino transfer
        PosAmp(2) = handles.controlTable.Data{1} + handles.controlTable.Data{2}; %get positive value for amplitude of the stimulator
        PosAmp(1) = floor(PosAmp(2));PosAmp(2) = (PosAmp(2)-PosAmp(1))*100; %modify for arduino transfer
        
        fwrite(handles.Arduino,[50 sCenter NegAmp PosAmp 1 num2str(1000/handles.controlTable.Data{5}) 'a' ...
            num2str(1000/handles.controlTable.Data{4}) 'a' num2str(handles.controlTable.Data{3}*1000)]); %send stim data to arduino
                          
        tic
        while handles.Arduino.BytesAvailable == 0  % wait for handshake for a max of 1s
            if toc > 2
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

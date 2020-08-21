function vStim_ArduinoMap(handles)
% Stimulation code to produce regular stimulus patterns with an arduino or 
% through Psychtoolbox directly. This is meant to be analyzed with an FFT
% approach later and some settings are comparable to vStim_DriftingNoise.
%
% StimTypes are the same as set in arduino control:
% 1 = Whisker stim
% 2 = Trunk stim
% 3 = Wheel stim
% 4 = Wheel jitter
% 5 = Audio stim
% 6 = Visual stim
% 7 = Audio + Visual stim
%
% Ask Simon if you have specific questions on how to use this code

%% Isolate all relevant variables from MainGUI
StatNames = {}; FlexNames = {};
StatVals = []; FlexVals = {};

for iVars = 1:size(handles.StaticVariableNames.String,1)
    temp = textscan(handles.StaticVariableNames.String{iVars},'%s%f');
    StatNames{iVars} = temp{1}{1};
    StatVals(iVars) = temp{2};
    clear temp
end
for iVars = 1:size(handles.FlexibleVariableNames.String,1)
    temp = textscan(handles.FlexibleVariableNames.String{iVars},'%s%f');
    FlexNames{iVars} = temp{1}{1};
    FlexVals{iVars} = str2num(handles.FlexibleVariableNames.String{iVars}(length(temp{1}{1})+1:end));
    clear temp
end

BasicVarNames = [StatNames FlexNames];  %all names basic variables that are required for correct function
FlexCases = combvec(FlexVals{:}); %get combinations for flexible variables
if isempty(FlexCases)
    FlexCases = 1; %at least one case, even if there are no flexible variables
end
BasicVarVals = zeros(length(BasicVarNames),size(FlexCases,2)); %values of each basic variable that is required. these can change based on the amount of cases that are produced.
    
for x = 1:length(FlexVals)
    BasicVarVals(ismember(BasicVarNames,FlexNames{x}),:) = FlexCases(x,:); % fill flexibe variable values from flexcases
end
for x = 1:length(StatNames)
    BasicVarVals(ismember(BasicVarNames,StatNames{x}),:) = repmat(StatVals(x),1,size(FlexCases,2)); % fill static variable values from StatVals
end

BasicVarVals = repmat(BasicVarVals,1,ceil(str2double(handles.NrTrials.String)/size(BasicVarVals,2))); %produce enough cases to cover all trials
if handles.RandomTrials.Value %randomize order of trials in each block of cases
    ind = [];
    for x = 1:ceil(str2double(handles.NrTrials.String)/size(FlexCases,2))
        ind = [ind randperm(size(FlexCases,2))+size(FlexCases,2)*(x-1)];
    end
    BasicVarVals = BasicVarVals(:,ind);
end

%% initialize Psychtoolbox and open screen
PsychDefaultSetup(1);
screenNumber = max(Screen('Screens')); % Draw to the external screen if avaliable
Background = mean(BasicVarVals(ismember(BasicVarNames,'Background'),:)); %background color. 
if Background < 0 || Background > 1; Background = 0; end % has be a number between 0 and 1 (defaults to black otherwise).
window = Screen('OpenWindow', screenNumber, Background); %open ptb window and save handle in pSetting
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Get the size of the current window
handles.ScreenRes.String = [num2str(screenXpixels) ' x ' num2str(screenYpixels)];
ifi=Screen('GetFlipInterval', window); %refresh rate
StimData.Paradigm = 'vStim_ArduinoMap'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals(:,1:str2double(handles.NrTrials.String)) = BasicVarVals(:,1:str2double(handles.NrTrials.String)); %basic variable values
StimData.ScreenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.handles = GetFigureData(handles); %get information from the figure handles
StimData.sRate = ifi;
HideCursor(window);

%% Produce required textures
handles.Status.String = 'Loading textures ...'; drawnow();
visPulseDur = max(BasicVarVals(ismember(BasicVarNames,'visPulseDur'),:))/1000; %max. visual pulse duration in seconds

% produce required textures when using different spatial/temporal frequencies or Stimlengths
tempFreq = ceil(visPulseDur/(1/max(BasicVarVals(ismember(BasicVarNames,'tempFreq'),:)))); %nr of textures per stimulus
nrCases = max(BasicVarVals(ismember(BasicVarNames,'trialDuration'),:))/(1/max(BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),:))); %nr of cycles per trial
vScreenWidth = (2*atan([StimData.ScreenSize{1} StimData.ScreenSize{3}]/(str2num(handles.EyeDistance.String)*2)))/pi*180; %screen size in vis. angle

if any(BasicVarVals(ismember(BasicVarNames,'useVisualNoise'),:)) && any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),:),[6 7])) % if visual stimulation + noise is requested
    for iCases = 1:(tempFreq*nrCases) %compute enough textures so that every pulse has a different pattern
        cTex = spatialPattern([screenXpixels, screenYpixels],-1,0.05,0.12,vScreenWidth,1)';
        cTex = cTex+abs(min(min(cTex)));cTex = cTex./max(max(cTex));
        noiseTex(iCases) =  Screen('MakeTexture', window, uint8(cTex*255));
        StimData.NoiseTexture{iCases} = cTex;
    end
end

%% Prepare audio data
if any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),:),[5 7])) || any(BasicVarVals(ismember(BasicVarNames,'BackgroundNoise'),:)) % if audio stimulation is required
    if any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),:),[5 7]))
        %make beep sound. One beep is 'audioPulseDur' ms long and repeated to match 'visPulseDur' at 'tempFreq'
        SamplingFreq = 192000;
        PsychToolBoxSound('init', SamplingFreq); %initialze soundcard
        
        if any(BasicVarVals(ismember(BasicVarNames,'useAudioNoise'),:)) %use noise bursts instead of click sounds
            beep = (rand(1,round(SamplingFreq*(min(BasicVarVals(ismember(BasicVarNames,'audioPulseDur'),:))/1000)))*2)-1; %white noise
        else %multifrequency convolved click
            tt = -pi:2*pi*1/(SamplingFreq*min(BasicVarVals(ismember(BasicVarNames,'audioPulseDur'),:))/1000):pi;
            tt=tt(1:end-1);
            beep = (1+cos(tt)).*(sin(2*tt)+sin(4*tt)+sin(6*tt)+sin(8*tt)+sin(16*tt)); clear tt
        end
        beep = beep./max(beep);
        beep = beep.*min(BasicVarVals(ismember(BasicVarNames,'audioAmp'),:)); %set loudness based on minimum audioAmp value
        
        gap = 1/max(BasicVarVals(ismember(BasicVarNames,'tempFreq'),:)); %gap between beeps
        if gap < visPulseDur
            beep(gap*SamplingFreq) = 0;
            beep = repmat(beep,1,floor(visPulseDur/gap));
        end
    else %background white noise
        SamplingFreq = 44100;
        PsychToolBoxSound('init', SamplingFreq); %initialze soundcard
        
        beep = (rand(1,round(SamplingFreq*(max(BasicVarVals(ismember(BasicVarNames,'trialDuration'),:)+1))))*2)-1; %background white noise. Duration of a trial + 1s. Noise is started at trial onset.
        beep = beep./max(beep);
        beep = beep.*min(BasicVarVals(ismember(BasicVarNames,'audioAmp'),:)); %set loudness based on minimum audioAmp value
    end
    
    PsychToolBoxSound('Load', [], beep); % upload sound to soundcard buffer
    
end


%% check arduino settings
[handles.Arduino,handles.ArduinoSettings] = vStim_ArduinoControl(handles.Arduino); %get arduino handle and settings
if ~isempty(handles.Arduino) && ~isempty(handles.ArduinoSettings)
    flushinput(handles.Arduino);
    sCenter(2) = handles.ArduinoSettings(1); %get center of the stimulator
    sCenter(1) = floor(sCenter(2));sCenter(2) = (sCenter(2)-sCenter(1))*100; %modify for arduino transfer
    NegAmp(2) = handles.ArduinoSettings(1) - min(BasicVarVals(ismember(BasicVarNames,'arduinoStimAmp'),:)); %get negative amplitude
    NegAmp(1) = floor(NegAmp(2));NegAmp(2) = (NegAmp(2)-NegAmp(1))*100; %modify for arduino transfer
    PosAmp(2) = handles.ArduinoSettings(1) + min(BasicVarVals(ismember(BasicVarNames,'arduinoStimAmp'),:)); %get positive amplitude
    PosAmp(1) = floor(PosAmp(2));PosAmp(2) = (PosAmp(2)-PosAmp(1))*100; %modify for arduino transfer
    
    fwrite(handles.Arduino,[50 sCenter NegAmp PosAmp 1 num2str(min(BasicVarVals(ismember(BasicVarNames,'arduinoPulseDur'),:))) 'a' ... %set up arduino single pulse amp and duration
        num2str(1000/min(BasicVarVals(ismember(BasicVarNames,'arduinoStimFreq'),:))) 'a' ... %set up arduino stim freq
        num2str(min(BasicVarVals(ismember(BasicVarNames,'arduinoStimDur'),:))*1000)]); %set up arduino stim sequence duration
    
    tic
    while handles.Arduino.BytesAvailable == 0  % wait for handshake for a max of 5s
        if toc > 5
            break
        end
    end
    if handles.Arduino.BytesAvailable == 0
        disp('Arduino did not respond')
        disp('Aborting paradigm')
        ResetHandles;
        return;
    else
        Byte = fread(handles.Arduino,1); %correct number is 10
        if Byte ~= 10
            disp(['Arduino sent wrong response on initialize: ' num2str(Byte) ' instead of 10'])
        end
    end    
else
    disp('Arduino not found or settings missing.')
    disp('Aborting paradigm')
    ResetHandles;
    return;
end
StimData.handles.ArduinoSettings = handles.ArduinoSettings; %add arduino settings to StimData

%% At this point the code should be ready to start stimulation.
%make sure trigger lines are set to false if serial port is present
if ~isempty(handles.SerialPort)
    IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is on during stimulation, second line goes on and off on each frame switch
end

%% Run paradigm
if isfield(handles,'DataPath')
    dPath = [handles.DataPath filesep]; %set data path
else
    dPath = [fileparts(which(handles.figure1.Name)) filesep]; %path to imager
end
dPath = [dPath handles.SubjectName.String filesep];
if ~isdir(dPath)
    mkdir(dPath);
end

% start audio background noise before starting trials for trialDur/2 or up to 10s, to give the mouse some time to get used to it.
if any(BasicVarVals(ismember(BasicVarNames,'BackgroundNoise'),:)) && ~any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),1),[5 7]))
    PsychToolBoxSound('Play'); %play background noise
    pause(min([10 max(BasicVarVals(ismember(BasicVarNames,'trialDuration'),:))/2]));
end

for iTrials = 1:str2double(handles.NrTrials.String)
    handles.Status.String = ['Running - Trial ' num2str(iTrials) ' ...'];drawnow();
    
    % start audio background noise before and after trial (to gap iti or time between two paradigms
    if BasicVarVals(ismember(BasicVarNames,'BackgroundNoise'),iTrials) && ~any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),iTrials),[5 7]))
        PsychToolBoxSound('Play'); %play background noise
    end
    
    % wait for baseline if value is assigned
    if BasicVarVals(ismember(BasicVarNames,'BaselineDur'),iTrials) > 0;
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=1,RTS=0'); %switch on stimulus line before actual sequence to allow to record some baseline data
        end
        tic; while BasicVarVals(ismember(BasicVarNames,'BaselineDur'),iTrials) > toc; end %wait for baseline to pass
    end
    
    StimData.TimeStamps{iTrials} = RunTrial(iTrials); %run current trial cycle
    
    if iTrials ~= str2double(handles.NrTrials.String)
        save([dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'],'StimData');
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...'];drawnow();
        try %blank screen after stimulus presentation
            Screen('FillRect', window,Background);
            if any(BasicVarVals(ismember(BasicVarNames,'BackgroundNoise'),:)) && ~any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),1),[5 7]))
                PsychToolBoxSound('Play'); %play background noise
            end
        catch
            ResetHandles
            return;
        end
        pause(str2double(handles.ITI.String)); %inter-trial interval
    else
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
        end
        if any(BasicVarVals(ismember(BasicVarNames,'BackgroundNoise'),:)) && ~any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),1),[5 7]))
            PsychToolBoxSound('Play'); %play background noise
        end
        save([dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'],'StimData');
    end
end

handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
handles.Status.String = 'vStim_ArduinoMap completed';
ResetHandles

%% nested function
function timeStamps = RunTrial(cTrial)
    % Animate moving bar
    if IsLinux
        Priority(9); % Set to realtime priority to ensure high drawing speed
    else
        Priority(MaxPriority(window));
    end
    
    trialDuration = BasicVarVals(ismember(BasicVarNames,'trialDuration'),cTrial);
    TrigSize = BasicVarVals(ismember(BasicVarNames,'visTriggerSize'),cTrial);
        
    Cnt = 1; %counter for timeStamps (increases with each refresh)
    fCnt = 1; %counter for frames (increases with time difference between frames divided by ifi)
    sCnt = 1; %counter for arduino stimulation (increases every cycle a serial command is sent)
    aCnt = 1; %counter for audio stimulation (increases every cycle a audio command is sent)
    tCnt = 1; %counter for textures (only increases with new texture)
    tCntLimit = round(1/BasicVarVals(ismember(BasicVarNames,'tempFreq'),cTrial)/ifi); %limit for when to use the next texture based on tempFreq
    visPulseDur = ceil(BasicVarVals(ismember(BasicVarNames,'visPulseDur'),cTrial)/1000/ifi); %duration of a visual pulse in #frames.
    cycleDuration = round((1/BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),cTrial))*round(1/ifi)); %duration of a cycle in #frames.
    timeStamps = zeros(1,trialDuration*round(1/ifi));
    FrameCounts = zeros(1,trialDuration*round(1/ifi));
    
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + trialDuration; %absolute time when trial should be ended
    
    while(cTime < absStimTime)
        
        stimOn = '0'; %flag to switch ttl line of there is ongoing stimulation        
         
        % select textures if visual stimulation is required
        cTex = []; %background screen if empty
        if rem(fCnt,cycleDuration) <= (cycleDuration/2+visPulseDur) && rem(fCnt,cycleDuration) > cycleDuration/2 %show visual stimulus at the middle of each cycle
            if any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),cTrial),[6 7])) % if visual stimulation is requested
                if BasicVarVals(ismember(BasicVarNames,'useVisualNoise'),cTrial) == 1
                    cTex = noiseTex(tCnt); % noise texture.
                    if rem(rem(fCnt,cycleDuration),tCntLimit) == 0
                        tCnt = tCnt +1; %switch noise texture if more than 'tCntLimit' of the same texture were presented already.
                    end
                else
                    cTex = NaN; %white screen if NaN
                end
            end
        else
            tCnt = 1; %reset texture counter
        end
        
        %determine screen content for current frame
        if isempty(cTex)
            Screen('FillRect', window, Background);
        elseif isnan(cTex)
            Screen('FillRect', window, 255);
            stimOn = '1'; %switch TTL if visual stimulus is being shown
        else
            Screen('DrawTexture', window, cTex);
            stimOn = '1'; %switch TTL if visual stimulus is being shown
        end
        
        %add visual indicator if ShowVisTrigger is true
        if BasicVarVals(ismember(BasicVarNames,'ShowVisTrigger'),cTrial) 
            if rem(fCnt,2) == 1 %show white square on even frame count
%                 Screen('FillRect', window,255,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
                Screen('FillRect', window,255,[0 0 TrigSize TrigSize])
            elseif rem(fCnt,2) == 0 %show black square on uneven frame count
%                 Screen('FillRect', window,0,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
                Screen('FillRect', window,0,[0 0 TrigSize TrigSize])
            end
        end
        
        %show frame
        cTime = Screen('Flip', window, cTime + 0.5 * ifi);
        timeStamps(Cnt) = cTime; 
        FrameCounts(Cnt) = fCnt;
        
        if any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),cTrial),[1:4 8:14])) % if arduino-based stimulation is requested
            if round(fCnt/cycleDuration) >= sCnt %serial command not sent in current cycle. Should trigger stimulation in the middle of each cycle.
                fwrite(handles.Arduino,BasicVarVals(ismember(BasicVarNames,'stimType'),cTrial)); %get arduino to start stimulus
                sCnt = sCnt+1; %increase counter for arduino stimulation
                stimOn = '1'; %switch TTL to indicate arduino stim being sent
            end
        end
        
        if any(ismember(BasicVarVals(ismember(BasicVarNames,'stimType'),cTrial),[5 7])) % if audio stimulation is requested
            if round(fCnt/cycleDuration) >= aCnt %audio not presented in current cycle. Should trigger stimulation in the middle of each cycle.
                PsychToolBoxSound('Play'); %play sound
                aCnt = aCnt+1; %increase counter for arduino stimulation
                stimOn = '1'; %switch TTL to indicate arduino stim being sent
            end
        end
        
        IOPort('ConfigureSerialPort',handles.SerialPort,['DTR=1,RTS=' stimOn]); %first line is trial running, second line is for ongoing stimulation

        %increase counter for timestamps and textures
        if Cnt == 1
            fCnt = fCnt + round(diff([sTime cTime])/ifi); %difference between first frame and start time makes counter value
        else
            fCnt = fCnt + round(diff([timeStamps(Cnt-1) timeStamps(Cnt)])/ifi); %difference between last and current frame makes counter value
        end
        Cnt = Cnt +1;
            
        drawnow(); %get current gui settings

        [keyIsDown, ~, keyCode, ~] = KbCheck; %check keyboard for button press
        if keyIsDown && strcmpi(KbName(find(keyCode)),'ESCAPE') %abort presentation
            handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
            ResetHandles
            return;
        end
    end
    
    %reset screen to background
    Screen('FillRect', window,Background);
    Screen('Flip', window);
    
    %produce triggers if serial port is present
    if ~isempty(handles.SerialPort)
        IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
    end

    Priority(0); % Set to normal priority
end

function ResetHandles
    handles.LoadParadigm.String = 'Load selected paradigm';
    handles.LoadParadigm.Value = false;
    handles.Status.String = 'No paradigm loaded';
    handles.StartParadigm.Enable = 'off';
    handles.PreviewParadigm.Enable = 'off';
    handles.StartParadigm.Value = false;
    handles.PreviewParadigm.Value = false;
    handles.StartParadigm.String = 'Start paradigm';
    handles.PreviewParadigm.String = 'Preview paradigm';
    Priority(0);
    drawnow();
    sca
    try
        PsychToolBoxSound('close'); %initialze soundcard
    catch
        disp('Error when closing sound server')
    end
end

function dataOut = GetFigureData(dataIn)
    %% Set variables in main GUI
    hFields = fieldnames(dataIn); %fields in main GUI
    for iFields = 1:length(hFields) %go through control panel elements
        if isprop(dataIn.(hFields{iFields}),'Type')
            if strcmp(get(dataIn.(hFields{iFields}),'Type'),'uicontrol') %if a field is a user control
                if strcmp(get(dataIn.(hFields{iFields}),'Style'),'edit') %if control is a textbox
                    dataOut.(hFields{iFields}) = dataIn.(hFields{iFields}).String;
                elseif strcmp(get(dataIn.(hFields{iFields}),'Style'),'togglebutton') %if control is a binary input
                    dataOut.(hFields{iFields}) = dataIn.(hFields{iFields}).Value;
                elseif strcmp(get(dataIn.(hFields{iFields}),'Style'),'popupmenu') %if control is a popup menu
                    for iCells = 1:length(dataIn.(hFields{iFields}).String)
                        dataOut.(hFields{iFields}){iCells} = dataIn.(hFields{iFields}).String{iCells};
                    end
                end
            end
        end
    end
end

function x = spatialPattern(DIM,BETA,fc,cutoff,visDegs,binary)
% function x = spatialPattern(DIM,BETA,fc,cutoff,visDegs)
%
% This function generates 1/f spatial noise, with a normal error 
% distribution (the grid must be at least 10x10 for the errors to be normal). 
% 1/f noise is scale invariant, there is no spatial scale for which the 
% variance plateaus out, so the process is non-stationary.

%% generate initial values and power spectrum
inc = visDegs ./ DIM;
for x = 1:2
    temp = (ceil(cutoff*visDegs(x))*cutoff^-1)-visDegs(x); %size in degrees to be added to approximate the cutoff correctly
    pixSize(x) = DIM(x) + round(temp / inc(x)); %pixel size for this dimension to give good result
    degSize(x) = visDegs(x) + round(temp / inc(x)) * inc(x); %size in degrees for extended image
end

%only use larger dimension to compute spatial noise
[pixSize,ind] = max(pixSize); 
degSize = degSize(ind);
cInc = inc(ind);

% produce grid of frequencies
u = [(0:floor(pixSize/2)) -(ceil(pixSize/2)-1:-1:1)]'/pixSize; %first dimension (collumns)
u = repmat(u,1,pixSize); % Reproduce these frequencies along every row
v = u'; %second dimension (rows)

% Generate the power spectrum
S_f = (sqrt(u.^2 + v.^2)+round(fc/cInc)).^(BETA);

%% remove frequencies above cutoff value
cutFreq = round(cutoff * degSize); %index for cutoff in pixels
if cutFreq > pixSize/2 %cutoff cant be more than half the image size
    cutFreq = pixSize/2;
end
S_f(cutFreq+1:end-cutFreq,:) = 0; %remove in between frequencies
S_f(:,cutFreq+1:end-cutFreq) = 0; %remove in between frequencies
S_f(S_f==inf) = 0; %make sure there are no inf values (usually at 0)

%create circle and remove frequencies along its radius
[x,y] = meshgrid(0:2*(cutFreq-1));
C = sqrt((x-ceil(cutFreq-1)).^2+(y-ceil(cutFreq-1)).^2);
C = C>C(cutFreq,end);

temp = S_f(1:cutFreq,1:cutFreq);
temp(C(cutFreq:end,cutFreq:end)) = 0;
S_f(1:cutFreq,1:cutFreq) = temp;

temp = S_f(end-cutFreq+1:end,1:cutFreq);
temp(C(1:cutFreq,cutFreq:end)) = 0;
S_f(end-cutFreq+1:end,1:cutFreq) = temp;

temp = S_f(1:cutFreq,end-cutFreq+1:end);
temp(C(cutFreq:end,1:cutFreq)) = 0;
S_f(1:cutFreq,end-cutFreq+1:end) = temp;

temp = S_f(end-cutFreq+1:end,end-cutFreq+1:end);
temp(C(1:cutFreq,1:cutFreq)) = 0;
S_f(end-cutFreq+1:end,end-cutFreq+1:end) = temp;
        
    
%% produce random data            
% Generate a grid of random phase shifts
phi = rand(pixSize);

% Inverse Fourier transform to obtain the the spatial pattern
x = ifft2(S_f .* (cos(2*pi*phi)+1i*sin(2*pi*phi)));

% Pick just the real component
x = real(x);
x = x(1:DIM(1),1:DIM(2));

%make binary mask if reqested
if binary
    x = x > median(x(:));
end
end

function PsychToolBoxSound(Function, SF, data)
% Modified version of the PsychToolboxSoundServer function that is part of
% BPod. This one is a bit simplified but has the same purpose. This code is
% made for use with the Xonar Soundcard and needs to be adjusted when
% using other hardware under Linux. 
% Only uses 2 instead of 8 channels and a single slave.

Function = lower(Function);
switch Function
    case 'init'
        PsychPortAudio('Verbosity', 0);
        InitializePsychSound(1);
        PsychPortAudio('Close');
        AudioDevices = PsychPortAudio('GetDevices');
        nDevices = length(AudioDevices);
        CandidateDevices = []; nCandidates = 0;
        if ispc
            for x = 1:nDevices
                if strcmp(AudioDevices(x).DeviceName, 'Primary Sound Driver')
                    nCandidates = nCandidates + 1;
                    CandidateDevices(nCandidates) = AudioDevices(x).DeviceIndex;
                    CandidateChannels(nCandidates) = AudioDevices(x).NrOutputChannels;
                end
            end
        else
            for x = 1:nDevices
                DeviceName = AudioDevices(x).DeviceName;
                if ~isempty(strfind(DeviceName, 'Xonar DX: Multichannel')) || ~isempty(strfind(DeviceName, 'Xonar U7: USB Audio')) % Assumes ASUS Xonar DX or U7 Soundcard
                    nCandidates = nCandidates + 1;
                    CandidateDevices(nCandidates) = AudioDevices(x).DeviceIndex;
                end
            end
        end
        if nCandidates > 0
            for x = 1:nCandidates
                disp(['Candidate device found! Trying candidate ' num2str(x) ' of ' num2str(nCandidates)])
                try
                    CandidateDevice = PsychPortAudio('Open', CandidateDevices(x), 9, 4, SF, 2);
                    PsychPortAudio('Close', CandidateDevice);
                    handles.SoundServer.SoundDeviceID = CandidateDevices(x);
                    disp('Success! A compatible sound card was detected.')
                catch
                    disp('Error: Detected sound card failed to open.')
                end
            end
        else
            disp('Error: no compatible sound subsystem detected. On Windows, ensure ASIO drivers are installed.')
        end
        
        handles.SoundServer.MasterOutput = PsychPortAudio('Open', handles.SoundServer.SoundDeviceID, 9, 4, SF, 2);
        PsychPortAudio('Start', handles.SoundServer.MasterOutput, 0, 0, 1);
        handles.SoundServer.SlaveOutput = PsychPortAudio('OpenSlave', handles.SoundServer.MasterOutput);
        PsychPortAudio('FillBuffer',handles.SoundServer.SlaveOutput(1), zeros(2,round(SF/1000)));
        PsychPortAudio('Start',handles.SoundServer.SlaveOutput(1));
        disp('PsychToolbox sound server successfully initialized.')
        
    case 'close'
        PsychPortAudio('Close');
        disp('PsychToolbox sound server successfully closed.')
        
    case 'load'
        Siz = size(data);
        if Siz(1) > 2
            error('Sound data must be a row vector and have not more than 2 channels');
        end
        
        %if audio data is not stereo, duplicate for both channels
        if Siz(1) == 1
            data(2,:) = data(1,:);
        end
        PsychPortAudio('FillBuffer', handles.SoundServer.SlaveOutput, data);
        
    case 'play'
        PsychPortAudio('Start', handles.SoundServer.SlaveOutput);
        
    case 'stop'
        PsychPortAudio('Stop', handles.SoundServer.SlaveOutput);
end
end
end
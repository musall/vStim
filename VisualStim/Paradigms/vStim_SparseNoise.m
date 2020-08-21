function vStim_SparseNoise(handles)
%code to produce sparse noise stimuli. Shows a set of black/white squares
%at multiple locations on screen to identify receptive fields for visually 
%responsive neurons. 
% Relevant inputs : 
%  - squareSize:        Size of each square in visual angles
%  - squareDistance:    Distance between squares in visual angles.
%  - stimDuration:      Duration of each individual pattern in miliseconds.
%  - trialDuration:     Duration of a single trial in seconds.
%
% Talk to Simon for more details.

%% Isolate all relevant variables from MainGUI
[BasicVarNames,BasicVarVals] = getGuiVars(handles);
BasicVarVals = max(BasicVarVals,[],2); %only use single set of values
nrTrials = str2double(handles.NrTrials.String); %total nr of trials
nrTex = 1000/BasicVarVals(strcmpi(BasicVarNames,'stimDuration'))*BasicVarVals(strcmpi(BasicVarNames,'trialDuration'))*nrTrials; %total nr of required textures
texPerTrial = nrTex/nrTrials; % nr of textures per trial

%% initialize Psychtoolbox and open screen
PsychDefaultSetup(1);
screenNumber = max(Screen('Screens')); % Draw to the external screen if avaliable
Background = round(max(BasicVarVals(strcmpi(BasicVarNames,'Background')))*255); % background color (0 = black, .5 = gray, 1 = white) 
screenInfo = Screen('Resolution',screenNumber); %get some information on the screen

if max([screenInfo.width, screenInfo.height]) > 1024
    StimData.resX = screenInfo.height / ceil(max([screenInfo.width, screenInfo.height])/1024) ;  % texture resolution, height
    StimData.resY = screenInfo.width / ceil(max([screenInfo.width, screenInfo.height])/1024) ;  % texture resolution, width
else
    StimData.resX = screenInfo.height;  % texture resolution, height
    StimData.resY = screenInfo.width;  % texture resolution, width
end

StimData.Paradigm = 'vStim_SparseNoise'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals = BasicVarVals; %basic variable values
StimData.handles = GetFigureData(handles); %get information from the figure handles

%stim variables
StimData.squareSize = max(BasicVarVals(strcmpi(BasicVarNames,'squareSize'))); %use largest enter square size for stimuli
StimData.degPerPixel = StimData.squareSize;  % Number of degrees per pixel of the undistorted image (real number)
StimData.minDist = max(BasicVarVals(strcmpi(BasicVarNames,'squareSpace'))); % Distance between squares (center to corner) in degrees
StimData.bright = [255 0 Background]; % Brightness of stimuli [white black background]
screenSize = textscan(handles.ScreenSizeAngle.String,'%f%c%f'); %get screen size in visual angles
StimData.anglesWidth = screenSize{1}; % width of screen, in visual angles
StimData.anglesHeight = screenSize{3}; % height of screen, in visual angles
screenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.physWidth  = screenSize{3};  % width of screen, in cm
StimData.physHeight = screenSize{1};  % height of screen, in cm

%eye position on screen
StimData.zDistBottom = str2double(handles.EyeLineDistBottom.String); % Distance to bottom of screen, along the horizontal eye line, in cm
StimData.zDistTop = str2double(handles.EyeLineDistTop.String); % Distance to top of screen, along the horizontal eye line, in cm
StimData.eyeX =  str2double(handles.EyeDistLeft.String);   % eye X location from left, in cm
StimData.eyeY = str2double(handles.EyeDistTop.String);   % eye Y location from top, in cm

%% compute some variables for stimulus presentation
minZ = min(StimData.zDistBottom, StimData.zDistTop);
wInDeg = atand(StimData.eyeX / minZ) + atand((StimData.physWidth - StimData.eyeX) / minZ);
wInSq = ceil(wInDeg / StimData.squareSize);
hInDeg = atand(StimData.eyeY / minZ) + atand((StimData.physHeight - StimData.eyeY) / minZ);
hInSq = ceil(hInDeg / StimData.squareSize);
screenSize = [hInDeg wInDeg];

% Compute where on the image to call the center
% Center this point in the undistorted image on the mouse's eye (degrees)
% Add +1 to account for padding in gridFromCoordsSingle
StimData.cartCX = wInSq * atand(StimData.eyeX / minZ) / wInDeg + 1;
StimData.cartCY = hInSq * atand(StimData.eyeY / minZ) / hInDeg + 1;

%% compute coordinates for sparse noise pattern and initialize PTB window
h = msgbox(['Generating ' num2str(nrTex) ' coordinates. Screen size is ' num2str(screenSize(2)) ' x ' num2str(screenSize(1)) '.'],'Creating coordinates','help');
delete(h.Children(1)); drawnow; %no 'ok' button
tic
StimData.nCoords = sparseNoisePatterns(screenSize, StimData.squareSize, StimData.minDist, nrTex);
toc
if ishandle(h);close(h);end
uiwait(msgbox('Coordinates generated. Press enter or space to continue.','Wait','modal'));

window = Screen('OpenWindow', screenNumber, Background); %open ptb window and save handle in pSettings
HideCursor(window);
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
[xRes, yRes] = Screen('WindowSize', window); % Get the size of the current window
handles.ScreenRes.String = [num2str(xRes) ' x ' num2str(yRes)];
ifi=Screen('GetFlipInterval', window); %refresh rate
StimData.sRate = ifi;

%make sure trigger lines are set to false if serial port is present
if ~isempty(handles.SerialPort)
    IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
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

for iTrials = 1:str2double(handles.NrTrials.String)
    
    % generate textures for current trial
    firstTex = (iTrials-1)*texPerTrial + 1; %first coordinate in current trial
    rawIms = gridFromCoordsSingle(ceil(screenSize ./ StimData.squareSize), StimData.nCoords(firstTex:firstTex + texPerTrial -1), StimData.bright);
    imStack = sphereDistortImage(rawIms, StimData); % Run spherical distortion
    imStack = rot90(imStack, 1);
%     imStack = permute(imStack,[2 1 3]);

    noiseTex = zeros(1,texPerTrial);
    for iTex = 1:texPerTrial
        noiseTex(iTex) = Screen('MakeTexture', window, imStack(:,:,iTex)); %create texture
    end
    
    %set up baseline if requested
    handles.Status.String = ['Running - Trial ' num2str(iTrials) ' ...'];drawnow();
    if BasicVarVals(ismember(BasicVarNames,'BaselineDur')) > 0;
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=1,RTS=0'); %switch on stimulus line before actual sequence to allow to record some baseline data
        end
        tic; while BasicVarVals(ismember(BasicVarNames,'baselineDur')) > toc; end %wait for baseline to pass
    end
    
    %run stimulation
    [StimData.TimeStamps{iTrials}, StimData.texID{iTrials}] = RunTrial;
    Screen('Close'); %delete textures of current trial from vRam
    
    %save stimdata to local and server path
    lPath = [dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat']; %local file path
    if iTrials ~= str2double(handles.NrTrials.String)
        save(lPath,'StimData');
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...'];drawnow();
        try %black screen after stimulus presentation
            Screen('FillRect', window,Background);
        catch
            ResetHandles
            return;
        end
        pause(str2double(handles.ITI.String)); %inter-trial interval
    else
        if ~isempty(handles.SerialPort) %produce triggers if serial port is present
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
        end
        save(lPath,'StimData');
    end
    
    %check if server is present and copy data to open folder if imager is running
    fPath = ['/mnt/Churchland/data/WidefieldImager/Animals/' handles.SubjectName.String '/PhaseMap/']; %path to server
    if IsLinux
        if isdir(fPath)
            if ~isempty(dir([fPath '*' date '*_open'])) %check for open folder. This should be indiciate for active session of the widefield system.
                a = dir([fPath '*' date '*_open']);
                fPath = [fPath a.name filesep handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'];
                copyfile(lPath,fPath); %copy local file to server
            end
        end
    end
end

if ~isempty(handles.SerialPort)
    IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
end
save([dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'],'StimData'); %save StimData

handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
handles.Status.String = 'vStim_SparseNoise completed';
ResetHandles

%% stim function
function [timeStamps, texID] = RunTrial % Animate sparse noise stimulation

    trigSize = BasicVarVals(strcmpi(BasicVarNames,'visTriggerSize'));
    trialDuration = BasicVarVals(strcmpi(BasicVarNames,'trialDuration'));
    stimDuration = BasicVarVals(strcmpi(BasicVarNames,'stimDuration'))/1000;
    timeStamps = zeros(1,trialDuration*round(1/ifi)); %timestamp for each frame
    texID = zeros(1,trialDuration*round(1/ifi)); %identify current texture for later analysis
    texTime = round(stimDuration/ifi); %duration of each stimulus in nr. of frames
    
    Cnt = 1;    %counter for timeStamps
    texCnt = 1; %counter for textures
    tCnt = 1;   %counter for textureID
    
    if IsWin
        Priority(2); % Set to realtime priority to ensure high drawing speed
    else
        Priority(9); % Set to realtime priority to ensure high drawing speed
    end
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + trialDuration;     

    while(cTime < absStimTime)
        
        % Draw sparse noise texture
        Screen('DrawTexture', window, noiseTex(texCnt), [], [0 0 xRes yRes]);
        
        %draw frame indicator
        if BasicVarVals(strcmpi(BasicVarNames,'showVisTrigger')) %add visual indicator if showVisTrigger is true
            if rem(Cnt,2) == 1 %show white square on even frame count
                Screen('FillRect', window,255 ,[0 0 trigSize trigSize])
%                 Screen('FillRect', window,255,[xRes-trigSize yRes-trigSize xRes yRes])
            else %show black square on uneven frame count
                Screen('FillRect', window,0 ,[0 0 trigSize trigSize])
%                 Screen('FillRect', window,0,[xRes-trigSize yRes-trigSize xRes yRes])
            end
        end
        
        %show frame
        cTime = Screen('Flip', window, cTime + 0.5 * ifi);
        timeStamps(Cnt) = cTime; 
        texID(Cnt) = (iTrials-1)*texPerTrial + texCnt; %identify current texture

        %increase counter for texture ID
        if Cnt == 1
            tCnt = tCnt + round(diff([sTime cTime])/ifi); %difference between first frame and start time makes counter value
        else
            tCnt = tCnt + round(diff([timeStamps(Cnt-1) timeStamps(Cnt)])/ifi); %difference between last and current frame makes counter value
        end
        Cnt = Cnt +1;
        
        %check if texture should be changed
        if texCnt < ceil(tCnt / texTime) && texCnt < length(noiseTex)
            texCnt = texCnt + 1;
        end
        
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,['DTR=1,RTS=' int2str(rem(Cnt,2))]); %first line is one during stimulation, second line goes on and off on each frame switch
        end

        [keyIsDown, ~, keyCode, ~] = KbCheck;
        if keyIsDown && any(strcmpi(KbName(find(keyCode)),'ESCAPE')) %abort presentation
            handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
            ResetHandles
            return;
        end
    end
    
    %produce triggers if serial port is present
    if ~isempty(handles.SerialPort)
        IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
    end
    
     %bank screen after stimulus presentation
    Screen('FillRect', window,Background);
    Screen('Flip', window); %
    Priority(0);
end

%% other small functions
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
end

function dataOut = GetFigureData(dataIn)
    % Set variables in main GUI
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
end
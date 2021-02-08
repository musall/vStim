function vStim_DriftingBar(handles)
% function for whole field mapping using drifting bar stimulus.
% Talk to Simon for details.

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

%% Open PTB screen window
PsychDefaultSetup(1);
screenNumber = max(Screen('Screens')); % Draw to the external screen if avaliable

Background = mean(BasicVarVals(ismember(BasicVarNames,'Background'),:)); %background color. 
if Background < 0 || Background > 1; Background = 0; end % has be a number between 0 and 1 (defaults to black otherwise).
window = Screen('OpenWindow', screenNumber, Background * 255); %open ptb window and save handle in pSettings
HideCursor(window);
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Get the size of the current window
handles.ScreenRes.String = [num2str(screenXpixels) ' x ' num2str(screenYpixels)];

% get some variables that will be saved later
ifi=Screen('GetFlipInterval', window); %refresh rate
StimData.Paradigm = 'vStim_DriftingBar'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals(:,1:str2double(handles.NrTrials.String)) = BasicVarVals(:,1:str2double(handles.NrTrials.String)); %basic variable values
StimData.ScreenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.handles = GetFigureData(handles); %get information from the figure handles
StimData.sRate = ifi;
 
%% Produce required textures
handles.Status.String = 'Loading textures ...'; drawnow();

% Convert square sizes from visual angles to pixel values
SquareSize = BasicVarVals(ismember(BasicVarNames,'SquareSize'),:)./180*pi; %convert visual angle to radians
SquareSize = tan(SquareSize./2).*(str2num(handles.EyeDistance.String)*2); %SquareSize in cm (or units of eye distance)
BasicVarVals(ismember(BasicVarNames,'SquareSize'),:) = ceil(SquareSize./(StimData.ScreenSize{1}./screenXpixels)); %barwidth in pixels
   
% Create square matrix and convert into texture. Second texture is an inverse so background can flicker around.
SquareSize = unique(BasicVarVals(ismember(BasicVarNames,'SquareSize'),:)); %get square sizes
for iSquares = 1:length(SquareSize)
    numSquares = ceil(max([screenXpixels screenYpixels])/SquareSize(iSquares));
    SquareMat = repmat([true(SquareSize(iSquares),SquareSize(iSquares));false(SquareSize(iSquares),SquareSize(iSquares))],ceil(numSquares/2),1); %
    SquareMat = repmat([SquareMat ~SquareMat],1,ceil(numSquares/2));
    SquareMat = SquareMat(1:screenYpixels,1:screenXpixels);
    
    SquaresTex(iSquares,1) = Screen('MakeTexture',window,uint8(SquareMat)*255); %create 8bit textures 
    
    if any(unique(BasicVarVals(ismember(BasicVarNames,'BarStyle'),:)) == 3) %if there is a case where BarStyle = 3, create an inverse texture as well
        SquaresTex(iSquares,2) = Screen('MakeTexture',window,uint8(~SquareMat)*255); %inverse of previous map
    end
end
WhiteTex = Screen('MakeTexture', window, 255); %small texture to save some memory. Can be streched into bar shape later.
 
%produce triggers if serial port is present
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

lPath = [dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat']; %path to local file

for iTrials = 1:str2double(handles.NrTrials.String)
    handles.Status.String = ['Running - Trial ' num2str(iTrials) ' ...'];drawnow();
    
    if BasicVarVals(ismember(BasicVarNames,'BaselineDur'),iTrials) > 0;
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=1,RTS=0'); %switch on stimulus line before actual sequence to allow to record some baseline data
        end
        if ~isempty(handles.Arduino)
            fwrite(handles.Arduino, handles.trialByte)
        end
        tic; while BasicVarVals(ismember(BasicVarNames,'BaselineDur'),iTrials) > toc; end %wait for baseline to pass
    end
    
    [StimData.TimeStamps{iTrials},StimData.barShift{iTrials},StimData.FrameCounts{iTrials}] = RunTrial(iTrials); %run current trial cycle
    
    if iTrials ~= str2double(handles.NrTrials.String)
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...'];drawnow();
        try %black screen after stimulus presentation
            Screen('FillRect', window,Background);
        catch
            ResetHandles
            return;
        end
        save(lPath,'StimData'); %this save wil not be included if the current trial was aborted prematurely
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

handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
handles.Status.String = 'vStim_DriftingBar completed';
ResetHandles

%% nested function
function [timeStamps,barShift,FrameCounts] = RunTrial(cTrial)
    % Animate moving bar
    if IsLinux
        Priority(9); % Set to realtime priority to ensure high drawing speed
    else
        Priority(MaxPriority(window));
    end
    
    BarStyle = BasicVarVals(ismember(BasicVarNames,'BarStyle'),cTrial);
    StimDuration = BasicVarVals(ismember(BasicVarNames,'StimDuration'),cTrial);
    TrigSize = BasicVarVals(ismember(BasicVarNames,'VisTriggerSize'),cTrial);
    
    BarWidth = BasicVarVals(ismember(BasicVarNames,'BarWidth'),cTrial)/180*pi; %convert visual angle to radians
    BarWidth = tan(BarWidth/2)*(str2num(handles.EyeDistance.String)*2); %barwidth in cm (or units of eye distance)
    BarWidth = ceil(BarWidth/(StimData.ScreenSize{1}/screenXpixels)); %barwidth in pixels
        
    destRect = [0 0 screenXpixels screenYpixels]; %barposition on the screen
    if BasicVarVals(ismember(BasicVarNames,'BarOrient'),cTrial) == 0 %check current bar orientation
        barInd = [1 3];
        screenWidth = screenXpixels + BarWidth;
    else
        barInd = [2 4];
        screenWidth = screenYpixels + BarWidth;
    end
    
    screenWidth = ceil(screenWidth / round(1/ifi))* round(1/ifi); %make sure requested screenwidth is a divider of refresh rate for even bar movement
    pixPerFrame = screenWidth/round(1/ifi)* BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),cTrial); %how many pixels per frame should be shifted to get the bar speed
    barShift = (0:(screenWidth/pixPerFrame)-1)*pixPerFrame; %index for where to shift bar in the next frame
    
    tCnt = 1; tCntLimit = round((1/BasicVarVals(ismember(BasicVarNames,'tempFreq'),cTrial))/ifi); %counter for textures
    timeStamps = zeros(1,StimDuration*round(1/ifi));
    FrameCounts = zeros(1,StimDuration*round(1/ifi));
    Cnt = 1; %counter for timeStamps
    
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + StimDuration;
    
    while(cTime < absStimTime)
        if tCnt > length(barShift)
            tCnt = tCnt-length(barShift); %reset texture count if end of screen is reached
        end
        
        if BarStyle == 1 %white bar
            cTex = WhiteTex; %use white bar texture
            srcRect = [0 0 1 1];
        elseif BarStyle == 2 %fixed squares
            cTex = SquaresTex(1); %use checkerboard
            srcRect = [0 0 screenXpixels screenYpixels];
            srcRect(barInd(2)) = BarWidth;
        elseif BarStyle == 3 %flickering squares
            cTex = SquaresTex(rem(ceil(tCnt/tCntLimit),2)+1); %switch texture every 'tCntLimit' cycles
            srcRect = [0 0 screenXpixels screenYpixels];
            srcRect(barInd(2)) = BarWidth;
        else
            error('Unknown bar style')
        end
        
        if BasicVarVals(ismember(BasicVarNames,'BarDirection'),cTrial) == 0
            destRect(barInd) = [barShift(tCnt)-BarWidth barShift(tCnt)]; %change bar destination coordinates
        else
            destRect(barInd) = [barShift(end-tCnt+1)-BarWidth barShift(end-tCnt+1)]; %change bar destination coordinates
        end
        
        if BasicVarVals(ismember(BasicVarNames,'UseAperture'),cTrial) && BarStyle ~= 1; srcRect = destRect; end %if true, bar is like an aperture on the grating instead of a moving bar with a texture on it
        
        %draw texture and flip screen
        Screen('DrawTexture', window, cTex,srcRect,destRect,0);
        if BasicVarVals(ismember(BasicVarNames,'ShowVisTrigger'),cTrial) %add visual indicator if ShowVisTrigger is true
            if rem(tCnt,2) == 1 %show white square on even frame count
%                 Screen('FillRect', window,255,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
                Screen('FillRect', window,255,[0 0 TrigSize TrigSize])
            elseif rem(tCnt,2) == 0 %show black square on uneven frame count
%                 Screen('FillRect', window,0,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
                Screen('FillRect', window,0,[0 0 TrigSize TrigSize])
            end
        end
        
        %show frame
        cTime = Screen('Flip', window, cTime + 0.5 * ifi);
        timeStamps(Cnt) = cTime; 
        FrameCounts(Cnt) = tCnt;
        
        % produce camByte on first stimulus
        if Cnt == 1 && ~isempty(handles.Arduino)
            fwrite(handles.Arduino, handles.camByte)
        end
        
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,['DTR=1,RTS=' int2str(rem(Cnt,2))]); %first line is one during stimulation, second line goes on and off on each frame switch
        end
        if rem(Cnt,2) == 1 && ~isempty(handles.Arduino)
            fwrite(handles.Arduino, handles.stimByte)
        end
        
        %increase counter for timestamps and textures
        if Cnt == 1
            tCnt = tCnt + round(diff([sTime cTime])/ifi); %difference between first frame and start time makes counter value
        else
            tCnt = tCnt + round(diff([timeStamps(Cnt-1) timeStamps(Cnt)])/ifi); %difference between last and current frame makes counter value
        end
        Cnt = Cnt +1;
            
        drawnow(); %get current gui settings

        [keyIsDown, ~, keyCode, ~] = KbCheck;
        if keyIsDown && strcmpi(KbName(find(keyCode)),'ESCAPE') %abort presentation
            handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
            ResetHandles
            return;
        end
    end    
    Screen('FillRect', window,Background);
    cTime = Screen('Flip', window);
    
    % stop camera trigger after sequence
    if Cnt == 1 && ~isempty(handles.Arduino)
        fwrite(handles.Arduino, handles.stopCamByte)
    end
        
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
    drawnow();
    Priority(0); % Set to normal priority
    %produce triggers if serial port is present
    if ~isempty(handles.SerialPort)
        IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
    end
    sca
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
end
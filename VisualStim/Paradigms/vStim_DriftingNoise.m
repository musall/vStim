function vStim_DriftingNoise(handles)
% function for whole field mapping using visual noise instead of checker
% board. Talk to Simon for details.

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

tempFreq = BasicVarVals(ismember(BasicVarNames,'tempFreq'),:) <= 0; %tempFreq <= 0 is not allowed. Show single texture instead.
BasicVarVals(ismember(BasicVarNames,'tempFreq'),tempFreq) = BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),tempFreq);

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
StimData.Paradigm = 'vStim_DriftingNoise'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals(:,1:str2double(handles.NrTrials.String)) = BasicVarVals(:,1:str2double(handles.NrTrials.String)); %basic variable values
StimData.ScreenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.handles = GetFigureData(handles); %get information from the figure handles
StimData.sRate = ifi;
 
%% Compute some variables about stimulus and screen size
vScreenWidth = (2*atan([StimData.ScreenSize{1} StimData.ScreenSize{3}]/(str2num(handles.EyeDistance.String)*2)))/pi*180; %screen size in vis. angle
screenWidth = max([screenXpixels,screenYpixels]); %find maximum screen width

% convert barwidth to visal angles and get values that will require the most textures
vBarWidth = max(BasicVarVals(ismember(BasicVarNames,'BarWidth'),:)); %find maximum barwidth in vis. angles
BarWidth = BasicVarVals(ismember(BasicVarNames,'BarWidth'),:)/180*pi; %convert visual angle to radians
BarWidth = tan(BarWidth/2)*(str2num(handles.EyeDistance.String)*2); %barwidth in cm (or units of eye distance)
BasicVarVals(ismember(BasicVarNames,'BarWidth'),:) = ceil(BarWidth/(StimData.ScreenSize{1}/screenXpixels)); %barwidth in pixels
minSpeed = min(BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),:)); %find minimum bar speed to produce enough textures
tempFreq = max(BasicVarVals(ismember(BasicVarNames,'tempFreq'),:)); %defines how long a given noise pattern is shown (in Hz). If tempFreq = screen refresh rate, a new noise pattern is shown on every frame.

%% Produce required textures
if any(BasicVarVals(ismember(BasicVarNames,'UseAperture'),:) == 0 & BasicVarVals(ismember(BasicVarNames,'WhiteBar'),:) == 0) %use spatial noise. don't need textures otherwise.
    handles.Status.String = 'Loading textures ...'; drawnow();
    barWidth = max(BasicVarVals(ismember(BasicVarNames,'BarWidth'),:)); %find maximum barwidth

    % compute noise textures
    for iTex = 1:tempFreq/minSpeed %get enough textures to change at 'tempFreq' during one cycle (e.g. 5 textures for 5Hz at 1s cycleDuration)
        cTex = spatialPattern([barWidth,screenWidth],-1,0.05,0.12,[vBarWidth,max(vScreenWidth)],1)';
        cTex = cTex+abs(min(min(cTex)));cTex = cTex./max(max(cTex));
        noiseTex(1,iTex) =  Screen('MakeTexture', window, uint8(cTex*255));
        noiseTex(2,iTex) =  Screen('MakeTexture', window, uint8(cTex'*255));
        StimData.NoiseTexture{iTex} = cTex;
    end
end

if any(BasicVarVals(ismember(BasicVarNames,'UseAperture'),:) == 1 & BasicVarVals(ismember(BasicVarNames,'WhiteBar'),:) == 0) %if aperture is requested, produce a full screen noise texture
    for iTex = 1:tempFreq/minSpeed %get enough textures to change at 'tempFreq' during one cycle (e.g. 5 textures for 5Hz at 1s cycleDuration)
        cTex = spatialPattern([screenXpixels,screenYpixels],-1,0.05,0.12,vScreenWidth,1)'; %full screen noise pattern
        cTex = cTex+abs(min(min(cTex)));cTex = cTex./max(max(cTex));
        fullNoiseTex(iTex) =  Screen('MakeTexture', window, uint8(cTex*255));
        StimData.FullNoiseTexture{iTex} = cTex;
    end
end

if any(BasicVarVals(ismember(BasicVarNames,'WhiteBar'),:) == 1)
    WhiteTex = Screen('MakeTexture', window, 255); %small texture to save some memory. Can be streched into bar shape later.
end
 
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

for iTrials = 1:str2double(handles.NrTrials.String)
    handles.Status.String = ['Running - Trial ' num2str(iTrials) ' ...'];drawnow();
    
    if BasicVarVals(ismember(BasicVarNames,'BaselineDur'),iTrials) > 0;
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=1,RTS=0'); %switch on stimulus line before actual sequence to allow to record some baseline data
        end
        tic; while BasicVarVals(ismember(BasicVarNames,'BaselineDur'),iTrials) > toc; end %wait for baseline to pass
    end
    
    RunTrial(iTrials); %run current trial cycle
        
    lPath = [dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat']; %local file path
    if iTrials ~= str2double(handles.NrTrials.String)
        save(lPath,'StimData');
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...'];drawnow();
        try %black screen after stimulus presentation
            Screen('FillRect', window,Background*255);
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

handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
handles.Status.String = 'vStim_DriftingNoise completed';
ResetHandles

%% nested function
function RunTrial(cTrial)
    % Animate moving bar
    if IsLinux
        Priority(9); % Set to realtime priority to ensure high drawing speed
    else
        Priority(MaxPriority(window));
    end
    
    %get some variables for current trial
    StimDuration = BasicVarVals(ismember(BasicVarNames,'StimDuration'),cTrial);
    TrigSize = BasicVarVals(ismember(BasicVarNames,'VisTriggerSize'),cTrial);
    WhiteBar = logical(BasicVarVals(ismember(BasicVarNames,'WhiteBar'),cTrial));
    UseAperture = logical(BasicVarVals(ismember(BasicVarNames,'UseAperture'),cTrial));
    UseFullField = logical(BasicVarVals(ismember(BasicVarNames,'UseFullField'),cTrial));
    BarWidth = BasicVarVals(ismember(BasicVarNames,'BarWidth'),cTrial);
    BarDirection = logical(BasicVarVals(ismember(BasicVarNames,'BarDirection'),cTrial));
    BarOrient = logical(BasicVarVals(ismember(BasicVarNames,'BarOrient'),cTrial));
    tempFreq = BasicVarVals(ismember(BasicVarNames,'tempFreq'),cTrial);
    
    %set up bar  position
    destRect = [0 0 screenXpixels screenYpixels]; %barposition on the screen
    if BarOrient %check current bar orientation
        barInd = [2 4];
        screenWidth = screenYpixels + BarWidth;
    else
        barInd = [1 3];
        screenWidth = screenXpixels + BarWidth;
    end
    
    screenWidth = ceil(screenWidth / round(1/ifi))* round(1/ifi); %make sure requested screenwidth is a divider of refresh rate for even bar movement
    pixPerFrame = screenWidth/round(1/ifi)* BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),cTrial); %how many pixels per frame should be shifted to get the bar speed
    barShift = (0:(screenWidth/pixPerFrame)-1)*pixPerFrame; %index for where to shift bar in the next frame
    
    tCnt = 1; tCntLimit = round((1/tempFreq)/ifi); %counter for textures
    timeStamps = zeros(1,StimDuration*round(1/ifi));
    StimData.FrameCounts{iTrials} = zeros(1,StimDuration*round(1/ifi));
    StimData.srcRect{iTrials} = zeros(StimDuration*round(1/ifi),4);
    StimData.destRect{iTrials} = zeros(StimDuration*round(1/ifi),4);
    Cnt = 1; %counter for timeStamps
    
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + StimDuration;
    
    while(cTime < absStimTime)
        if tCnt > length(barShift); 
            tCnt = tCnt-length(barShift); %reset texture count if end of screen is reached
        end
        
        if BarDirection
            destRect(barInd) = [barShift(tCnt)-BarWidth barShift(tCnt)]; %change bar destination coordinates
        else
            destRect(barInd) = [barShift(end-tCnt+1)-BarWidth barShift(end-tCnt+1)]; %change bar destination coordinates
        end
        
        if UseFullField
            destRect = [0 0 screenXpixels screenYpixels]; %in case of full field, the stimulus should cover the whole screen
            UseAperture = true;
        end
                
        %pick texture and source coordinates
        if WhiteBar %use white bar
            cTex = WhiteTex; %use white bar texture
            srcRect = [0 0 1 1];
        else
            if UseAperture || UseFullField
                cTex = fullNoiseTex(ceil(tCnt/round((1/ifi)/tempFreq))); %full field noise stimulus. pick texture according to current tempFreq.
                srcRect = destRect;
            else
                cTex = noiseTex(BarOrient+1,ceil(tCnt/round((1/ifi)/tempFreq))); %bar noise stimulus. pick texture according to current tempFreq.
                srcRect = [0 0 screenXpixels screenYpixels];
                srcRect(barInd(2)) = BarWidth;
            end
        end
         
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
        
        %save source and target coordinates of bar texture to stimdata
        StimData.srcRect{iTrials}(Cnt,:) = srcRect;
        StimData.destRect{iTrials}(Cnt,:) = destRect;
        
        %show frame
        cTime = Screen('Flip', window, cTime + 0.5 * ifi);
        timeStamps(Cnt) = cTime; 
        StimData.FrameCounts{iTrials}(Cnt) = tCnt;
        
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,['DTR=1,RTS=' int2str(rem(Cnt,2))]); %first line is one during stimulation, second line goes on and off on each frame switch
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
    Screen('FillRect', window,Background*255);
    Screen('Flip', window);
    
    %produce triggers if serial port is present
    if ~isempty(handles.SerialPort)
        IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
    end
    
    % set up trialdata to be saved
    StimData.barShift{iTrials} = barShift;
    StimData.TimeStamps{iTrials} = timeStamps;
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
u = repmat(u,1,pixSize); % Reproduce these frequencies along ever row
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
end
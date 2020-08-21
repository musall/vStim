function vStim_DriftingMorphedBar(handles)
% background, cyclesPerSecond and stimDuration can't be flexible variables

%% Isolate all relevant variables from MainGUI
[BasicVarNames,BasicVarVals] = getGuiVars(handles);

%% initialize Psychtoolbox and open screen
PsychDefaultSetup(1);
screenNumber = max(Screen('Screens')); % Draw to the external screen if avaliable
background = round(BasicVarVals(strcmpi(BasicVarNames,'background'),1)*255); % background color (0 = black, .5 = gray, 1 = white) 
screenInfo = Screen('Resolution',screenNumber); %get some information on the screen
texPerTrial = 1/BasicVarVals(strcmpi(BasicVarNames,'cyclesPerSecond'),1)*screenInfo.hz; %total nr of required textures
StimData.sRate = screenInfo.hz;
handles.ScreenRes.String = [num2str(screenInfo.width) ' x ' num2str(screenInfo.height)];

if max([screenInfo.width, screenInfo.height]) > 1024
    StimData.resY = screenInfo.height / ceil(max([screenInfo.width, screenInfo.height])/1024) ;  % texture resolution, height
    StimData.resX = screenInfo.width / ceil(max([screenInfo.width, screenInfo.height])/1024) ;  % texture resolution, width
else
    StimData.resY = screenInfo.height;  % texture resolution, height
    StimData.resX = screenInfo.width;  % texture resolution, width
end

StimData.Paradigm = 'vStim_DriftingMorphedBar'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals = BasicVarVals; %basic variable values
StimData.handles = GetFigureData(handles); %get information from the figure handles
barOrientations = unique(BasicVarVals(strcmpi(BasicVarNames,'barOrient'),:));

%stim variables
StimData.squareSize = max(BasicVarVals(strcmpi(BasicVarNames,'squareSize'))); %use largest enter square size for stimuli
screenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.physWidth  = screenSize{1};  % width of screen, in cm
StimData.physHeight = screenSize{3};  % height of screen, in cm
barWidth = BasicVarVals(ismember(BasicVarNames,'barWidth'),1);

%eye position on screen
StimData.zDistBottom = str2double(handles.EyeLineDistBottom.String); % Distance to bottom of screen, along the horizontal eye line, in cm
StimData.zDistTop = str2double(handles.EyeLineDistTop.String); % Distance to top of screen, along the horizontal eye line, in cm
StimData.eyeX =  str2double(handles.EyeDistLeft.String);   % eye X location from left, in cm
StimData.eyeY = str2double(handles.EyeDistTop.String);   % eye Y location from top, in cm

%% compute some variables for stimulus presentation
StimData.degPerPixel = 0.1; %resolution of visal angle image
minZ = min(StimData.zDistTop, StimData.zDistBottom);
wInDeg = atand(StimData.eyeX / minZ) + atand((StimData.physWidth - StimData.eyeX) / minZ);
wInSq = ceil(wInDeg / StimData.degPerPixel);
hInDeg = atand(StimData.eyeY / minZ) + atand((StimData.physHeight - StimData.eyeY) / minZ);
hInSq = ceil(hInDeg / StimData.degPerPixel);
screenSize = [hInDeg wInDeg];

% Compute where on the image to call the center
% Center this point in the undistorted image on the mouse's eye (degrees)
% Add +1 to account for padding in gridFromCoordsSingle
StimData.cartCX = wInSq * atand(StimData.eyeX / minZ) / wInDeg + 1;
StimData.cartCY = hInSq * atand(StimData.eyeY / minZ) / hInDeg + 1;

%% initialize PTB window
window = Screen('OpenWindow', screenNumber, background); %open ptb window and save handle in pSettings
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
ifi=Screen('GetFlipInterval', window); %refresh rate
HideCursor(window);
[xRes, yRes] = Screen('WindowSize', window); % Get the size of the current window

%make sure trigger lines are set to false if serial port is present
if ~isempty(handles.SerialPort)
    IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
end

%% compute textures and upload to graphics card
degPerFrame = ceil((max(screenSize+barWidth)/screenInfo.hz*BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),1))*10)/10; %how many degrees per frame should be shifted to get the right bar speed. ceiled after 2 decimals
texPerCycle = screenInfo.hz/BasicVarVals(ismember(BasicVarNames,'cyclesPerSecond'),1); %textures required per cycle
barTex = zeros(length(barOrientations),texPerCycle);

% create bar texture if using white or checkerboard stimulus
if BasicVarVals(ismember(BasicVarNames,'barStyle'),1) <= 0
    cBar = ones(ceil([max(screenSize) barWidth]./StimData.degPerPixel)+2); %white bar
elseif BasicVarVals(ismember(BasicVarNames,'barStyle'),1) >= 1
    squareSize = BasicVarVals(ismember(BasicVarNames,'squareSize'),1)/StimData.degPerPixel; %square size in pixels
    numSquares = ceil(((ceil(max(screenSize)./StimData.degPerPixel) + 2)./2)/squareSize); %maximum number of squares
    cBar = checkerboard(squareSize,numSquares,numSquares) > .5;
end

for iOrients = 1:length(barOrientations)
    for iTex = 1:texPerCycle %number of textures, needed for 1 cyle
        if rem(iTex,round((1/BasicVarVals(ismember(BasicVarNames,'tempFreq'),1))/(1/screenInfo.hz))) == 1 %create new bar pattern
            if BasicVarVals(ismember(BasicVarNames,'useSpatialNoise'),1) == 1
                cBar = spatialPattern(ceil([max(screenSize) barWidth]./StimData.degPerPixel)+2,-1,0.05,0.2,[max(screenSize) barWidth],1); %noise pattern
            else
                cBar = ~cBar; %flip checkerboard/white texture
            end
        end
        
        cTex = ones(ceil(screenSize./StimData.degPerPixel)+2).*BasicVarVals(strcmpi(BasicVarNames,'background'),1); %blank texture
        cShift = ((iTex-1)*degPerFrame)/StimData.degPerPixel + 1; %current offset for bar motion
        cShift = round(cShift - barWidth/StimData.degPerPixel + 1:cShift);
        cShift(cShift < 1) = [];
        
        if iOrients == 1
            cShift(cShift > size(cTex,2)) = [];
            cTex(:,cShift) = cBar(1:size(cTex,1),1:length(cShift));
        else
            cShift(cShift > size(cTex,1)) = [];
            cTex(cShift,:) = cBar(1:size(cTex,2),1:length(cShift))';
        end
        
        cTex = sphereDistort(cTex,StimData); % perform spherical distortion
        barTex(iOrients,iTex) =  Screen('MakeTexture', window, uint8(cTex*255)); %upload texture to VRAM

    end
end
clear cBar cTex

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
    %set up baseline if requested
    handles.Status.String = ['Running - Trial ' num2str(iTrials) ' ...'];drawnow();
    if BasicVarVals(ismember(BasicVarNames,'baselineDur'),1) > 0;
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=1,RTS=0'); %switch on stimulus line before actual sequence to allow to record some baseline data
        end
        if ~isempty(handles.Arduino)
            fwrite(handles.Arduino, handles.trialByte)
        end
        tic; while BasicVarVals(ismember(BasicVarNames,'baselineDur'),1) > toc; end %wait for baseline to pass
    end
    
    %run stimulation
    [StimData.TimeStamps{iTrials}, StimData.texID{iTrials}] = RunTrial;
    
    %save stimdata to local and server path
    lPath = [dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat']; %local file path
    if iTrials ~= str2double(handles.NrTrials.String)
        save(lPath,'StimData');
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...'];drawnow();
        try %blank screen after stimulus presentation
            Screen('FillRect', window,background);
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
handles.Status.String = 'vStim_DriftingMorphedBar completed';
ResetHandles

%% stim function
function [timeStamps, texID] = RunTrial % Animate drifting bar stimulation

    trigSize = BasicVarVals(strcmpi(BasicVarNames,'visTriggerSize'),iTrials);
    trialDuration = BasicVarVals(strcmpi(BasicVarNames,'trialDuration'),iTrials);
    timeStamps = zeros(1,trialDuration*round(1/ifi)); %timestamp for each frame
    texID = zeros(1,trialDuration*round(1/ifi)); %identify current texture for later analysis
    Cnt = 1;    %counter for timeStamps
    tCnt = 1;    %counter for textures

    if BasicVarVals(strcmpi(BasicVarNames,'barDirection'),iTrials) == 0
        barShift = 1:size(barTex,2); %forward motion
    else
        barShift = size(barTex,2):-1:1; %reverse motion
    end
    
    if IsWin
        Priority(2); % Set to realtime priority to ensure high drawing speed
    else
        Priority(9); % Set to realtime priority to ensure high drawing speed
    end
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + trialDuration;     

    while(cTime < absStimTime)
        
        if tCnt > length(barShift) % end of sequence is reached
            tCnt = tCnt - length(barShift); %reset texture counter
        end
        
        % Draw current texture
        if BasicVarVals(strcmpi(BasicVarNames,'barOrient'),iTrials) == 0
            Screen('DrawTexture', window, barTex(1,barShift(tCnt)), [], [0 0 xRes yRes]);
        else
            Screen('DrawTexture', window, barTex(2,barShift(tCnt)), [], [0 0 xRes yRes]);
        end
               
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
        texID(Cnt) = (iTrials-1)*texPerTrial + tCnt; %identify current texture

        %increase counter for texture ID
        if Cnt == 1
            tCnt = tCnt + round(diff([sTime cTime])/ifi); %difference between first frame and start time makes counter value
        else
            tCnt = tCnt + round(diff([timeStamps(Cnt-1) timeStamps(Cnt)])/ifi); %difference between last and current frame makes counter value
        end
        Cnt = Cnt +1;
        
        %produce triggers if serial port is present
        if ~isempty(handles.SerialPort)
            IOPort('ConfigureSerialPort',handles.SerialPort,['DTR=1,RTS=' int2str(rem(Cnt,2))]); %first line is one during stimulation, second line goes on and off on each frame switch
        end
        if rem(Cnt,2) == 1 && ~isempty(handles.Arduino)
            fwrite(handles.Arduino, handles.stimByte)
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
    Screen('FillRect', window,background);
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
    drawnow();
    if ~isempty(handles.SerialPort)
        IOPort('ConfigureSerialPort',handles.SerialPort,'DTR=0,RTS=0'); %first line is one during stimulation, second line goes on and off on each frame switch
    end
    Priority(0);
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

function cTex = sphereDistort(cTex,StimData) % perform spherical distortion
    top = StimData.physHeight - StimData.eyeY;
    left = -StimData.eyeX;

    % Project monitor pixels into spherical coordinate
    [xi, yi] = meshgrid(linspace(0, 1, StimData.resX), linspace(1, 0, StimData.resY));
    cartPointsXFull = left + StimData.physWidth .* xi;
    cartPointsYFull = top - StimData.physHeight .* yi;
    cartPointsZFull = StimData.zDistBottom + (StimData.zDistTop - StimData.zDistBottom) .* yi;
    [sphrPointsTh, sphrPointsPh] = cart2sph(cartPointsZFull, cartPointsXFull, cartPointsYFull);
    
    % Recenter image pixels in spherical coordinates
    [initResX, initResY] = size(cTex);
    
    [xi, yi] = meshgrid(0:initResY-1, 0:initResX-1);
    imX = pi / 180 * StimData.degPerPixel * (xi - StimData.cartCX + 0.5);
    imY = pi / 180 * StimData.degPerPixel * (yi - StimData.cartCY + 0.5);
    
    % Use interp2 to find nearest image pixel for each monitor pixel
    cTex = interp2(imX, imY, cTex, sphrPointsTh, sphrPointsPh, 'nearest', background);
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
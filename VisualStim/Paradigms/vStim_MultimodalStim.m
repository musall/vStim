function vStim_MultimodalStim(handles)
% code to produce sensory stimuli for multisensory stimulation.
% UseAperture creates a circular aperture. Otherwise makes a full field stimulus. 
% This paradigm needs a detected analog outout module. Talk to Simon for details.

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
FlexCases = CombVec(FlexVals{:}); %get combinations for flexible variables
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

Screen('Preference', 'SkipSyncTests', 0);
Background = mean(BasicVarVals(ismember(BasicVarNames,'Background'),:))*255; %background color. 
window = Screen('OpenWindow', screenNumber, Background); %open ptb window and save handle in pSettings
TrigSize = BasicVarVals(ismember(BasicVarNames,'VisTriggerSize'),1);
Screen('FillRect', window, 0, [0 0 TrigSize TrigSize]); %make indicator black
HideCursor(window);
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Get the size of the current window
handles.ScreenRes.String = [num2str(screenXpixels) ' x ' num2str(screenYpixels)];
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % Enable alpha blending
ifi=Screen('GetFlipInterval', window); %refresh rate
StimData.Paradigm = 'vStim?MultimodalStim'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals(:,1:str2double(handles.NrTrials.String)) = BasicVarVals(:,1:str2double(handles.NrTrials.String)); %basic variable values
StimData.ScreenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.handles = GetFigureData(handles); %get information from the figure handles
StimData.sRate = ifi;

% modify aperture size from visual angle to pixel
ScreenSize = textscan(handles.ScreenSize.String,'%f'); %get horizontal screen size in cm
texsize = BasicVarVals(ismember(BasicVarNames,'ApertureSize'),:)./180*pi; %convert visual angle to radians
texsize = tan(texsize./2)*str2num(handles.EyeDistance.String); %texsize in cm
BasicVarVals(ismember(BasicVarNames,'ApertureSize'),:) = round(texsize./(ScreenSize{1}/screenXpixels)); %aperture size in pixels

%% Produce required textures
% produce gratings. This is modified from the 'DriftDemo2' code. Check there for explanations
handles.Status.String = 'Loading textures ...'; drawnow();

%convert spatial frequency input from cpd (cycles per degree) to onscreen pixels per cycle
SpatialFreqs = 1./BasicVarVals(ismember(BasicVarNames,'SpatialFreq'),:); %spatial frequency in visual degrees per cycle
SpatialFreqs = SpatialFreqs./180*pi; %convert visual angle to radians
SpatialFreqs = tan(SpatialFreqs./2)*str2num(handles.EyeDistance.String); %cm per cycle
BasicVarVals(ismember(BasicVarNames,'SpatialFreq'),:) = round(SpatialFreqs./(ScreenSize{1}/screenXpixels)); %pixels per cycle

% produce required textures when using different spatial/temporal frequencies or Stimlengths
SpatialFreqs = unique(BasicVarVals(ismember(BasicVarNames,'SpatialFreq'),:));
TemporalFreqs = unique(BasicVarVals(ismember(BasicVarNames,'TemporalFreq'),:));
UseApertures = unique(BasicVarVals(ismember(BasicVarNames,'UseAperture'),:));
ApertureSizes = unique(BasicVarVals(ismember(BasicVarNames,'ApertureSize'),:));
FlexNames = {'SpatialFreq','TemporalFreq','UseAperture','ApertureSize'};
FlexCases = CombVec(SpatialFreqs,TemporalFreqs,UseApertures,ApertureSizes); %possible combinations from above variables


for iCases = 1:size(FlexCases,2)
    p=ceil(FlexCases(1,iCases));
    fr=1/FlexCases(1,iCases)*2*pi;
    
    if FlexCases(3,iCases) == 1
        texsize = round(FlexCases(4,iCases)/2); %half of texsize in pixels
    else
        texsize = round(sqrt(screenXpixels^2+screenYpixels^2)/2);
    end
    x = meshgrid(-texsize:texsize + p, 1);
    grating= 127 + 127*cos(fr*x);
    gratingtex(iCases) = Screen('MakeTexture', window, grating);
end

if any(UseApertures == 0)
    ApertureSizes = [ApertureSizes sqrt(screenXpixels^2+screenYpixels^2)]; %add mask for full field stimuli if needed
end
for iCases = 1:length(ApertureSizes)
    texsize = round(ApertureSizes(iCases)/2);
    mask=ones(2*texsize+1, 2*texsize+1, 2) * Background;
    [x,y] = meshgrid(0:2*texsize);
    C = sqrt((x-ceil(texsize)).^2+(y-ceil(texsize)).^2);
    mask(:, :, 2) = ~(C<=C(texsize+1,end))*255;
    masktex(iCases)=Screen('MakeTexture', window, mask);
end

%% load stimuli to analog output module
optoDur = unique(BasicVarVals(ismember(BasicVarNames,'optoDur'),:));

handles.WavePlayer.loadWaveform(1,ones(1,5000));
handles.WavePlayer.loadWaveform(2,rand(1,5000));

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
    StimData.TimeStamps{iTrials} = RunTrial(iTrials); %run current trial cycle
    
    if iTrials ~= str2double(handles.NrTrials.String)
        save([dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'],'StimData');
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...']; drawnow();
        try %blank screen after stimulus presentation
            Screen('FillRect', window, Background);
        catch
            ResetHandles
            return;
        end
        pause(str2double(handles.ITI.String)); %inter-trial interval
    else
        save([dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'],'StimData');
    end
end

handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
handles.Status.String = 'vStim_DriftingGradients completed';
ResetHandles

%% nested function
function timeStamps = RunTrial(cTrial) % Animate drifting gradients
    %identify case for grating texture
    temp = zeros(length(FlexNames),size(FlexCases,2));
    for x = 1:length(FlexNames)
        temp(x,:) = FlexCases(x,:) == BasicVarVals(ismember(BasicVarNames,FlexNames{x}),cTrial);
    end
    [~,iCase] = max(sum(temp)); %correct case should be the one where FlexCases and BasicVarVals most agree
    if BasicVarVals(ismember(BasicVarNames,'UseAperture'),cTrial) == 1
        mCase = ApertureSizes == BasicVarVals(ismember(BasicVarNames,'ApertureSize'),cTrial); %mask for current aperture size
    else
        mCase = length(ApertureSizes);
    end
    % translate visual angle position into screen coordinates
    xPosition = BasicVarVals(ismember(BasicVarNames,'xPosition'),cTrial); %xPosition in pixels
    yPosition = BasicVarVals(ismember(BasicVarNames,'yPosition'),cTrial); %yPosition in pixels
    temp = [xPosition yPosition]./180*pi; %convert visual angle to radians
    temp = tan(temp./2)*str2num(handles.EyeDistance.String); %positions in cm
    temp = round(temp./(ScreenSize{1}/screenXpixels)); %positions of stimulus center in pixels
    xPosition = temp(1); yPosition = temp(2);
    
    TrigSize = BasicVarVals(ismember(BasicVarNames,'VisTriggerSize'),cTrial);
    StimDuration = BasicVarVals(ismember(BasicVarNames,'StimDuration'),cTrial);
    TemporalFreq = BasicVarVals(ismember(BasicVarNames,'TemporalFreq'),cTrial);
    SpatialFreq = BasicVarVals(ismember(BasicVarNames,'SpatialFreq'),cTrial);
    pixPerFrame = SpatialFreq / round(1/ifi) * TemporalFreq; %how many pixels per frame should be shifted to get the bar speed

    if BasicVarVals(ismember(BasicVarNames,'UseAperture'),cTrial) == 1
        visibleSize = round(BasicVarVals(ismember(BasicVarNames,'ApertureSize'),cTrial)/2)*2+1;
    else
        visibleSize = sqrt(screenXpixels^2+screenYpixels^2)+1;
    end
    srcRect = [0 0 visibleSize visibleSize]; %rectangle from source texture
    [temp(1),temp(2)] = RectCenter([0 0 screenXpixels screenYpixels]);
    destRect = CenterRectOnPointd(srcRect,xPosition+temp(1),-yPosition+temp(2)); %rectangle on screen
       
    timeStamps = zeros(1,StimDuration*round(1/ifi));
    Cnt = 0; %counter for timeStamps
    if IsWin
        Priority(2); % Set to realtime priority to ensure high drawing speed
    else
        Priority(9); % Set to realtime priority to ensure high drawing speed
    end
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + StimDuration;     
    
    while(cTime < absStimTime)
        
        xoffset = mod(Cnt*pixPerFrame,SpatialFreq);
        srcRect=[xoffset 0 xoffset + visibleSize visibleSize];
        Cnt = Cnt+1;
        
        % Draw grating texture, rotated by "angle":
        Screen('DrawTexture', window, gratingtex(iCase), srcRect, destRect, BasicVarVals(ismember(BasicVarNames,'StimAngle'),cTrial),[],[],[],[],kPsychUseTextureMatrixForRotation);
        %Draw aperture
        Screen('DrawTexture', window, masktex(mCase), [0 0 visibleSize visibleSize], destRect);
        
        %draw frame indicator
        if BasicVarVals(ismember(BasicVarNames,'ShowVisTrigger'),cTrial) %add visual indicator if ShowVisTrigger is true
            if rem(Cnt,2) == 1 %show white square on even frame count
%                 Screen('FillRect', window,255,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
                Screen('FillRect', window,255,[0 0 TrigSize TrigSize])
            else %show black square on uneven frame count
%                 Screen('FillRect', window,0,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
                Screen('FillRect', window,0,[0 0 TrigSize TrigSize])
            end
        end
        
        %show frame
        cTime = Screen('Flip', window, cTime + 0.5 * ifi);
        timeStamps(Cnt) = cTime;
        
        if Cnt == 1
            handles.WavePlayer.play(1,1);
            handles.WavePlayer.play(2,2);
        end

        [keyIsDown, ~, keyCode, ~] = KbCheck;
        if keyIsDown && any(strcmpi(KbName(find(keyCode)),'ESCAPE')) %abort presentation
            handles.ExperimentNr.String = num2str(str2double(handles.ExperimentNr.String)+1); %update experiment counter
            ResetHandles
            return;
        end
    end
    
     %blank screen after stimulus presentation
    Screen('FillRect', window,Background);
    Screen('FillRect', window, 0, [0 0 TrigSize TrigSize]); %make indicator black
    Screen('Flip', window); %
    Priority(0);
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
                    if iscell(dataIn.(hFields{iFields}))
                        for iCells = 1:length(dataIn.(hFields{iFields}).String)
                            dataOut.(hFields{iFields}){iCells} = dataIn.(hFields{iFields}).String{iCells};
                        end
                    else
                        dataOut.(hFields{iFields}){1} = dataIn.(hFields{iFields}).String;
                    end
                end
            end
        end
    end
end
end
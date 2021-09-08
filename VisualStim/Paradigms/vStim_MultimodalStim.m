function vStim_MultimodalStim(handles)
% code to produce sensory stimuli for multisensory stimulation.
% UseAperture creates a circular aperture. Otherwise makes a full field stimulus. 
% This paradigm needs a detected analog outout module. 
% Stimtype defines the type of sensory stimulation.
%
% Stimtype 1 only vision, no analog output needed
% Stimtype 2; only audio, channel 1
% Stimtype 3; only tactile, channel 2
% Stimtype 4; vision-audio, channel 1
% Stimtype 5; vision-tactile, channel 2
% Stimtype 6; audio-tactile, channel 1+2
% Stimtype 7; vision-audio-tactile, channel 1+2
% Stimtype 8; empty trial. writes one '0' to channel 0.
% Stimtype 9; control tactile stimulus. creates a tactile stimulus on a second valve that is not directed at the whisker pad.

% Variables that begin with 'optoCases' identify cases where optogenetic
% stimulation should be combiend with sensory stimulation. The values of
% these variables identify the Stimtype for which an optogenetic
% stimulation case should be added.
%
% Talk to Simon for more details.

ignoreEmptyTrials = true; %flag to not produce trials where nothing is shown (Stimtype 8 without laser)

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

%% isolate optogenetic cases and put in a separate category
optoID = {'OptoCases' 'RedPower' 'BluePower'}; %identifier
idx1 = contains(StatNames, optoID); %find static optocases
idx2 = contains(FlexNames, optoID); %find flexible optocases

optoNames = [StatNames(idx1), FlexNames(idx2)];
optoNames = cellfun(@(x) strrep(x,'OptoCases',''), optoNames, 'UniformOutput',false); %remove identifer from variable names
optoVars = [num2cell(StatVals(idx1)), FlexVals(idx2)];

% remove flagged cases from old arrays
StatNames(idx1) = [];
StatVals(idx1) = [];
FlexNames(idx2) = [];
FlexVals(idx2) = [];

%% generate cases based on other variables
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

% add opto variables as zeros
BasicVarNames = [BasicVarNames optoNames]; 
BasicVarVals = [BasicVarVals; zeros(length(optoNames), size(BasicVarVals,2))];

%% add specific optogenetic cases
stimIdx = ismember(BasicVarNames, 'StimType'); %index for stimtype
nonOptoCases = 1:size(BasicVarVals,2);
for x = 1 : length(optoNames)
    if ~contains(optoNames{x}, {'RedPower' 'BluePower'})
        % find cases and check if corresponding stimtype exists
        if any(ismember(unique(optoVars{x}), unique(BasicVarVals(stimIdx,nonOptoCases))))
            
            %find all cases with matching stimtypes and duplicate
            cIdx = ismember(BasicVarVals(stimIdx,nonOptoCases), unique(optoVars{x}));
            newCases = BasicVarVals(:, nonOptoCases(cIdx));
            newCases(size(newCases,1) - length(optoNames) + x, :) = 1;
            
            % check for optoPower
            switch optoNames{x}
                case 'BlueOne'; powerIdx = find(contains(optoNames, 'BluePowerOne'));
                case 'RedOne';  powerIdx = find(contains(optoNames, 'RedPowerOne'));
                case 'BlueTwo'; powerIdx = find(contains(optoNames, 'BluePowerTwo'));
                case 'RedTwo';  powerIdx = find(contains(optoNames, 'RedPowerTwo'));
            end
            
            powerVals = optoVars{powerIdx}; %amplitudes for current laser
            newCases = repmat(newCases, 1, length(powerVals)); %make more new cases if there is more than one amplitude
            powerVals = repmat(powerVals, size(newCases,2)/length(powerVals),1);
            newCases(size(newCases,1) - length(optoNames) + powerIdx, :) = powerVals(:); %asign amplitude to new cases
            
            BasicVarVals = [BasicVarVals, newCases]; %add these to possible cases
        end
    end
end

%% remove empty case without laser
if ignoreEmptyTrials
    useIdx = ~(BasicVarVals(stimIdx,:) == 8 & sum(BasicVarVals(end-length(optoNames)+1:end, :)) == 0);
    BasicVarVals = BasicVarVals(:, useIdx);
end
nrCases = size(BasicVarVals,2); %nr of possible cases

%% make enough cases to match requested trialcount
BasicVarVals = repmat(BasicVarVals,1,ceil(str2double(handles.NrTrials.String)/size(BasicVarVals,2))); %produce enough cases to cover all trials
if handles.RandomTrials.Value %randomize order of trials in each block of cases
    ind = [];
    for x = 1:ceil(str2double(handles.NrTrials.String)/nrCases)
        ind = [ind randperm(nrCases)+nrCases*(x-1)];
    end
    BasicVarVals = BasicVarVals(:,ind);
end

%% load stimuli to analog output module
W = BpodWavePlayer(handles.WavePlayerPort);
W.OutputRange = '-5V:5V'; % make sure output range is correct
W.TriggerMode = 'Master';
W.TriggerProfileEnable = 'On'; % use trigger profiles to produce different waveforms across channels
W.TriggerMode = 'Master'; %output can be interrupted by new stimulus triggers
W.LoopDuration(1:8) = 0; %keep on for a up to 10 minutes

analogRate = unique(BasicVarVals(ismember(BasicVarNames,'AnalogRate'),:)); analogRate = analogRate(1); %sampling rate of analog signals
pulseCount = unique(BasicVarVals(ismember(BasicVarNames,'PulseCount'),:)); pulseCount = pulseCount(1); %number of sensory event
pulseDur = unique(BasicVarVals(ismember(BasicVarNames,'PulseDur'),:)); pulseDur = pulseDur(1); %duration of single sensory event
pulseGap = unique(BasicVarVals(ismember(BasicVarNames,'PulseGap'),:)); pulseGap = pulseGap(1); %gap between sensory events

% auditory pulses
% noise = rand(1,pulseDur*analogRate)-0.5; %single burst
clickLength = 3; %click is usually 3ms long
tt = -pi:2*pi*1000/(analogRate*clickLength):pi; tt(end) = [];
click = (1+cos(tt)).*(sin(2*tt)+sin(4*tt)+sin(6*tt)+sin(8*tt)+sin(16*tt));
noise = zeros(1,pulseDur*analogRate);
noise(1:length(click)) = click / max(click) / 2;
noise([1 end]) = 0;

audioAmp = unique(BasicVarVals(ismember(BasicVarNames,'AudioAmp'),:)); audioAmp = audioAmp(1); %amplitude if auditory noise
puffDur = unique(BasicVarVals(ismember(BasicVarNames,'PuffDur'),:)); puffDur = puffDur(1); %duration of air puffs (only used if shorter as pulsedur)
makePuff = unique(BasicVarVals(ismember(BasicVarNames,'MakePuff'),:)); makePuff = makePuff(1) == 1; %flag to use air puffs or not

% tactile
if makePuff % tactile air puffs
    puffDur = min([puffDur, pulseDur]);
    puff = ones(1,puffDur*analogRate); puff(end) = 0;
else
    % tactile buzzes from actuator
    resFreq = 205; %resonant frequency of tactile actuator
    temp = -pi:2*pi/(analogRate*pulseDur):pi; temp(end) = [];
    puff = (cos(temp)+1)/2;
    puff = (puff(1,:) - min(puff)) / max((puff(1,:) - min(puff)));
    cycles = ceil(2*pi*resFreq*pulseDur);
    tt = 0:cycles/round(analogRate*pulseDur):cycles; tt(end-1) = [];
    puff(length(tt)) = 0; %make sure puff is long enough
    puff = puff.*sin(tt);
end

% load stimuli to output module
W.loadWaveform(1,noise * 5 * audioAmp); %signal 1 is auditory
W.loadWaveform(2,puff * 5); %signal 2 is tactile

% optogenetic stimulus
optoPyramid = unique(BasicVarVals(ismember(BasicVarNames,'OptoPyramid'),:)); optoPyramid = optoPyramid(1); %flag to make excitatory pyramid instead of square wave stimuli
optoDur = unique(BasicVarVals(ismember(BasicVarNames,'OptoDur'),:)); optoDur = optoDur(1); %duration of optogenetic stimulus
optoRamp = unique(BasicVarVals(ismember(BasicVarNames,'OptoRamp'),:)); optoRamp = optoRamp(1); %duration of ramp after square wave
optoFreq = unique(BasicVarVals(ismember(BasicVarNames,'OptoFreq'),:)); optoFreq = optoFreq(1); %frequency of square wave stimulus
redPower1 = unique(BasicVarVals(ismember(BasicVarNames,'RedPowerOne'),:)); redPower1 = redPower1(1); %power of red laser 1
redPower2 = unique(BasicVarVals(ismember(BasicVarNames,'RedPowerTwo'),:)); redPower2 = redPower2(1); %power of red laser 2
bluePower1 = unique(BasicVarVals(ismember(BasicVarNames,'BluePowerOne'),:)); bluePower1 = bluePower1(1); %power of blue laser 1
bluePower2 = unique(BasicVarVals(ismember(BasicVarNames,'BluePowerTwo'),:)); bluePower2 = bluePower2(1); %power of blue laser 2
redPower1 = max([min([redPower1(1),1]), 0]).* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
redPower2 = max([min([redPower2(1),1]), 0]).* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
bluePower1 = max([min([bluePower1(1),1]), 0]).* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
bluePower2 = max([min([bluePower2(1),1]), 0]).* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V

optoStim = vStim_getOptoStim(analogRate, optoDur, optoRamp, optoFreq, optoRamp); %get waveforms for optogenetics
optoStim = optoStim .* 3.3; %scale to 3.3V output
% make sure pluse duration is a divider of optoFreq when using pulsed stimuli
% optoPulseDur = ceil(pulseDur / (1/optoFreq)) * (1/optoFreq);
optoPulseDur = pulseDur;

% load stimuli to output module
redSeqSignals = [3 4]; %keep which signals are used for red sequences
W.loadWaveform(3,optoStim(1,:) * redPower1); %signal 3 is red 1 
W.loadWaveform(4,optoStim(1,:) * redPower2); %signal 4 is red 2

blueSeqSignals = [5 6]; %keep which signals are used for blue sequences
W.loadWaveform(5,optoStim(2,:) * bluePower1); %signal 5 is blue 1
W.loadWaveform(6,optoStim(2,:) * bluePower2); %signal 6 is blue 2

if optoPyramid == 1
    step = 1/round(analogRate*optoPulseDur/2);
    pyraPulse = [step : step : 1 1 - step: -step : 0]; %pyramid
else
    pyraPulse = ones(1,round(analogRate*optoPulseDur)); %square wave
end

redPulseSignals = [7 8]; %keep which signals are used for red pulses
W.loadWaveform(redPulseSignals(1), pyraPulse * redPower1); %signal 7 is red pulse or pyramid 1
W.loadWaveform(redPulseSignals(2), pyraPulse * redPower2); %signal 8 is red pulse or pyramid 2

seqEnable = 9; %keep which signal is used to enable lasers sequence
W.loadWaveform(seqEnable,ones(1,size(optoStim,2)) .* 3.3); %signal 9 is enable trigger

pulseEnable = 10; %keep which signal is used to enable laser pulse
W.loadWaveform(pulseEnable,ones(1,optoPulseDur * analogRate) .* 3.3); %signal 10 is pulse enable trigger

W.loadWaveform(11,0); %signal 11 is empty

% create trigger profiles that match different stimype
redCases = {[3 7] [5 8] [3 5 7 8]};
for stimTypes = 1 : 9
    switch stimTypes
        case 1 %only vision, no analog output needed
        case 2; W.TriggerProfiles(stimTypes, 1) = 1; %only audio, channel 1
        case 3; W.TriggerProfiles(stimTypes, 2) = 2; %only tactile, channel 2
        case 4; W.TriggerProfiles(stimTypes, 1) = 1; %vision-audio, channel 1
        case 5; W.TriggerProfiles(stimTypes, 2) = 2; %vision-tactile, channel 2
        case 6; W.TriggerProfiles(stimTypes, 1:2) = 1:2; %audio-tactile, channel 1+2
        case 7; W.TriggerProfiles(stimTypes, 1:2) = 1:2; %vision-audio-tactile, channel 1+2
        case 8; W.TriggerProfiles(stimTypes, 1) = 11; %empty trial
        case 9; W.TriggerProfiles(stimTypes, 2) = 2; %tactile control trial, channel 2
    end
    
    for redLEDs = 1:3
        % make extra cases where red light is triggered with sensory
        % stimuli. make the same triggerprofile but add red pulses. This
        % will presumably range between profiles 11 and 37.
        if redLEDs == 3
            redSignal = [redPulseSignals, repmat(pulseEnable, 1, length(redCases{redLEDs})/2)];
        else
            redSignal = [redPulseSignals(redLEDs), repmat(pulseEnable, 1, length(redCases{redLEDs})/2)];
        end
        switch stimTypes
            case 1; W.TriggerProfiles(redLEDs*10 + stimTypes, redCases{redLEDs}) = redSignal; %only vision
            case 2; W.TriggerProfiles(redLEDs*10 + stimTypes, [1 redCases{redLEDs}]) = [1 redSignal]; %only audio, channel 1
            case 3; W.TriggerProfiles(redLEDs*10 + stimTypes, [2 redCases{redLEDs}]) = [2 redSignal]; %only tactile, channel 2
            case 4; W.TriggerProfiles(redLEDs*10 + stimTypes, [1 redCases{redLEDs}]) = [1 redSignal]; %vision-audio, channel 1
            case 5; W.TriggerProfiles(redLEDs*10 + stimTypes, [2 redCases{redLEDs}]) = [2 redSignal]; %vision-tactile, channel 2
            case 6; W.TriggerProfiles(redLEDs*10 + stimTypes, [1:2 redCases{redLEDs}]) = [1:2 redSignal]; %audio-tactile, channel 1+2
            case 7; W.TriggerProfiles(redLEDs*10 + stimTypes, [1:2 redCases{redLEDs}]) = [1:2 redSignal]; %vision-audio-tactile, channel 1+2
            case 8; W.TriggerProfiles(redLEDs*10 + stimTypes, redCases{redLEDs}) = redSignal; %red pulses only
            case 9; W.TriggerProfiles(redLEDs*10 + stimTypes, [2 redCases{redLEDs}]) = [2 redSignal]; %tactile control trial, channel 2
        end
    end
end

% define trigger profiles for optogenetic cases
W.TriggerProfiles(51, [3 7]) = [redSeqSignals(1), seqEnable]; %red laser on location 1 (RedOne)
W.TriggerProfiles(52, [4 7]) = [blueSeqSignals(1), seqEnable]; %blue laser on location 1 (BlueOne)
W.TriggerProfiles(53, [5 8]) = [redSeqSignals(2), seqEnable]; %red laser on location 2 (RedTwo)
W.TriggerProfiles(54, [6 8]) = [blueSeqSignals(2), seqEnable]; %blue laser on location 2 (BlueTwo)
W.TriggerProfiles(55, [3 5 7 8]) = [redSeqSignals, seqEnable, seqEnable]; %red laser on location 1+2 (RedBoth)
W.TriggerProfiles(56, [4 6 7 8]) = [blueSeqSignals, seqEnable, seqEnable]; %blue laser on location 1+2 (BlueBoth)

handles.WavePlayer = W; %make sure this is the same

%% initialize Psychtoolbox and open screen
PsychDefaultSetup(1);
screenNumber = max(Screen('Screens')); % Draw to the external screen if avaliable

TrigSize = BasicVarVals(ismember(BasicVarNames,'VisTriggerSize'),1);
Screen('Preference', 'SkipSyncTests', 1);
Background = mean(BasicVarVals(ismember(BasicVarNames,'Background'),:))*255; %background color. 
window = Screen('OpenWindow', screenNumber, Background); %open ptb window and save handle in pSettings
% window = Screen('OpenWindow', screenNumber, Background, [0 0 TrigSize TrigSize]*2); %open ptb window and save handle in pSettings
Screen('FillRect', window, 0, [0 0 TrigSize TrigSize]); %make indicator black
HideCursor(window);
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Get the size of the current window
handles.ScreenRes.String = [num2str(screenXpixels) ' x ' num2str(screenYpixels)];
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % Enable alpha blending
ifi=Screen('GetFlipInterval', window); %refresh rate

% adjust optoshift to be a divider of the video framerate
BasicVarVals(ismember(BasicVarNames,'OptoShift'),:) = round(BasicVarVals(ismember(BasicVarNames,'OptoShift'),:) ./ ifi) * ifi; %

StimData.Paradigm = 'vStim_MultimodalStim'; %name of the paradigm
StimData.VarNames = BasicVarNames; %basic variable names
StimData.VarVals(:,1:str2double(handles.NrTrials.String)) = BasicVarVals(:,1:str2double(handles.NrTrials.String)); %basic variable values
StimData.ScreenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
StimData.handles = GetFigureData(handles); %get information from the figure handles
ifi = 1 / round(1/ifi); %get more precise number on inter-frame interval
StimData.sRate = ifi;

% modify aperture size from visual angle to pixel
ScreenSize = textscan(handles.ScreenSize.String,'%f'); %get horizontal screen size in cm
texsize = BasicVarVals(ismember(BasicVarNames,'ApertureSize'),:)./180*pi; %convert visual angle to radians
texsize = tan(texsize./2)*str2num(handles.EyeDistance.String); %texsize in cm
BasicVarVals(ismember(BasicVarNames,'ApertureSize'),:) = round(texsize./(ScreenSize{1}/screenXpixels)); %aperture size in pixels

% load noise texture if needed
useVisNoise = unique(BasicVarVals(ismember(BasicVarNames,'UseVisNoise'),:)); useVisNoise = useVisNoise(1) == 1; %flag to use pre-computed noise stimulus
if useVisNoise
    noiseFile = [fileparts(which(handles.figure1.Name)) filesep 'Misc' filesep 'gaussian_noise.mat';];
    load(noiseFile, 'gaussian_noise');
    StimData.noiseFile = noiseFile; %save path to stimulus file for later reference
    srcNoise = [0 0 size(gaussian_noise,2) size(gaussian_noise,1)];
    for x = 1 : size(gaussian_noise,3)
        noiseTex(1,x) =  Screen('MakeTexture', window, uint8(gaussian_noise(:,:,x)));
    end
end

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

Screen('FillRect', window, 0, [0 0 TrigSize TrigSize]); %make indicator black
Screen('Flip', window);
java.lang.Thread.sleep(5000); %wait for 5 seconds before moving forward
for iTrials = 1:str2double(handles.NrTrials.String)
    handles.Status.String = ['Running - Trial ' num2str(iTrials) ' ...'];drawnow();
    StimData.TimeStamps{iTrials} = RunTrial(iTrials); %run current trial cycle
    
    if iTrials ~= str2double(handles.NrTrials.String)
        save([dPath handles.SubjectName.String '_' date '_' handles.ExperimentNr.String '_settings.mat'],'StimData');
        handles.Status.String = ['Trial ' num2str(iTrials) ' done - Waiting...']; drawnow();
        try %blank screen after stimulus presentation
            Screen('FillRect', window, Background);
            Screen('FillRect', window, 0, [0 0 TrigSize TrigSize]); %make indicator black
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
handles.Status.String = 'vStim_MultimodalStim completed';
ResetHandles

%% nested function
function timeStamps = RunTrial(cTrial) % Animate drifting gradients
    
    % get current stimtype and optocase and determine analog outputs
    cStim = BasicVarVals(ismember(BasicVarNames,'StimType'),cTrial); %get stimtype
    visualOn = ismember(cStim, [1 4 5 7]); %flag to present visual stimulus
    
    % get laser power
    redPower1 = unique(BasicVarVals(ismember(BasicVarNames,'RedPowerOne'),cTrial)); %power of red laser 1
    redPower1 = max([min([redPower1(1),1]), 0]) .* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
    redPower2 = unique(BasicVarVals(ismember(BasicVarNames,'RedPowerTwo'),cTrial)); %power of red laser 2
    redPower2 = max([min([redPower2(1),1]), 0]) .* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
    bluePower1 = unique(BasicVarVals(ismember(BasicVarNames,'BluePowerOne'),cTrial)); %power of blue laser 1
    bluePower1 = max([min([bluePower1(1),1]), 0]) .* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
    bluePower2 = unique(BasicVarVals(ismember(BasicVarNames,'BluePowerTwo'),cTrial)); %power of blue laser 2
    bluePower2 = max([min([bluePower2(1),1]), 0]) .* 3.3; %make sure its a value between 0 and 1 and mulitply with 3.3V
    
    % get current optocase and determine analog outputs
    RedSensoryPulses = unique(BasicVarVals(ismember(BasicVarNames,'RedSensoryPulses'),cTrial));
    cOptoVals = BasicVarVals(ismember(BasicVarNames,optoNames),cTrial)>0; %get active optogenetic cases
    optoOut = []; %optogenetic output profile
    if ~isempty(optoNames(cOptoVals))
        switch optoNames{cOptoVals}
            case 'RedOne'
                optoOut = 51;  %red laser on location 1
                W.loadWaveform(3,optoStim(1,:) * redPower1); %signal 3 is red 1 
            case 'BlueOne'
                optoOut = 52; %blue laser on location 1
                W.loadWaveform(5,optoStim(2,:) * bluePower1); %signal 5 is blue 1
            case 'RedTwo'
                optoOut = 53; %red laser on location 2
                W.loadWaveform(4,optoStim(1,:) * redPower2); %signal 4 is red 2
            case 'BlueTwo'
                optoOut = 54; %blue laser on location 2
                W.loadWaveform(6,optoStim(2,:) * bluePower2); %signal 6 is blue 2
            case 'RedBoth'
                optoOut = 55; %red laser on location 1 and 2
                W.loadWaveform(3,optoStim(1,:) * redPower1); %signal 3 is red 1 
                W.loadWaveform(4,optoStim(1,:) * redPower2); %signal 4 is red 2
            case 'BlueBoth'
                optoOut = 56; %blue laser on location 1 and 2
                W.loadWaveform(5,optoStim(2,:) * bluePower1); %signal 5 is blue 1
                W.loadWaveform(6,optoStim(2,:) * bluePower2); %signal 6 is blue 2
        end
        if RedSensoryPulses
            switch optoNames{cOptoVals}
                case 'RedOne'
                    optoOut = []; cStim = cStim + 10; %red laser on location 1 but with sensory stimuli
                    W.loadWaveform(redPulseSignals(1), pyraPulse * redPower1); %signal 7 is red pulse or pyramid 1
                case 'RedTwo'; optoOut = []; cStim = cStim + 20; %red laser on location 2  but with sensory stimuli
                    W.loadWaveform(redPulseSignals(2), pyraPulse * redPower2); %signal 8 is red pulse or pyramid 2
                case 'RedBoth'; optoOut = []; cStim = cStim + 30; %red laser on location 1 and 2 but with sensory stimuli
                    W.loadWaveform(redPulseSignals(1), pyraPulse * redPower1); %signal 7 is red pulse or pyramid 1
                    W.loadWaveform(redPulseSignals(2), pyraPulse * redPower2); %signal 8 is red pulse or pyramid 2
            end
        end
        
        if ~isempty(handles.Arduino)
            % check byte for enable signal
            if contains(optoNames{cOptoVals}, 'Red')
                enableByte = 1; %red lasers
            else
%                 enableByte = 2; %blue lasers
                enableByte = 3; %violet lasers
            end
            
            % switch enable lines through connected teensy (if present). Send one
            % value for each laser to identify enabled wavelength.
            % 1 = red, 2 = cyan, 3 = violet.
            fwrite(handles.Arduino, [150 enableByte enableByte]);
        end
    end
    
    % check byte for valve signal
    % switch lines for valve 1 and 2 through connected teensy (if present).
    if cStim == 9 %use valve 2 if this is a tactile control trial. Otherwise, use valve 1
        fwrite(handles.Arduino, [151 2]);
    else
        fwrite(handles.Arduino, [151 1]);
    end
            
    %identify case for grating texture
    temp = zeros(length(FlexNames),nrCases);
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
    TemporalFreq = BasicVarVals(ismember(BasicVarNames,'TemporalFreq'),cTrial);
    SpatialFreq = BasicVarVals(ismember(BasicVarNames,'SpatialFreq'),cTrial);
    pixPerFrame = SpatialFreq / round(1/ifi) * TemporalFreq; %how many pixels per frame should be shifted to get the bar speed

    useVisNoise = BasicVarVals(ismember(BasicVarNames,'UseVisNoise'),cTrial) == 1;

    if BasicVarVals(ismember(BasicVarNames,'UseAperture'),cTrial) == 1
        visibleSize = round(BasicVarVals(ismember(BasicVarNames,'ApertureSize'),cTrial)/2)*2+1;
    else
        visibleSize = sqrt(screenXpixels^2+screenYpixels^2)+1;
    end
    srcRect = [0 0 visibleSize visibleSize]; %rectangle from source texture
    [temp(1),temp(2)] = RectCenter([0 0 screenXpixels screenYpixels]);
    destRect = CenterRectOnPointd(srcRect,xPosition+temp(1),-yPosition+temp(2)); %rectangle on screen
    
    
    %assess timing and duration of current trial (duration of optogenetic stimulus is not taken into account)
    optoShift = BasicVarVals(ismember(BasicVarNames,'OptoShift'),cTrial); %shift of optogenetics relative to sensory stimulus
    optoShift = optoShift * any(BasicVarVals(ismember(BasicVarNames,optoNames),cTrial) > 0); %only shift in optogenetic trials
    pulseDur = BasicVarVals(ismember(BasicVarNames,'PulseDur'),cTrial); %duration of a sensory pulse
    pulseGap = BasicVarVals(ismember(BasicVarNames,'PulseGap'),cTrial); %gap between pulses
    stimDur = (pulseDur + pulseGap) * pulseCount;
    if optoShift < 0
        stimDur = stimDur - optoShift;
    end
    
    timeStamps = zeros(1,ceil(stimDur*round(1/ifi)));
    Cnt = 0; %counter for timeStamps
    if IsWin
        Priority(2); % Set to realtime priority to ensure high drawing speed
    else
        Priority(9); % Set to realtime priority to ensure high drawing speed
    end
    sTime = Screen('Flip', window); % Sync start time to the vertical retrace
    cTime = sTime; %current time equals start time
    absStimTime = sTime + stimDur; %total trial duration
    pulseStart = cTime - (pulseDur + pulseGap); %first pulse can start immediately
    pulseCnt = 0;
    trigerAnalog = false;
    trigerCamera = false;
    
    if optoShift < 0
        sensoryStart = cTime - optoShift; %shift start time of first sensory pulse
        optoStart = cTime; %optogenetics starts immediately
    else
        optoStart = cTime + optoShift; %shift optogenetics relative to stim onset
        sensoryStart = cTime;
    end
    
    noiseCnt = 0; %counter for noise Texture

    % start running stimulation
    while(cTime < absStimTime)
        
        xoffset = mod(Cnt*pixPerFrame,SpatialFreq);
        srcRect=[xoffset 0 xoffset + visibleSize visibleSize];
        Cnt = Cnt+1;
        
        pulseOn = ((cTime - pulseStart) >= pulseDur + pulseGap) && pulseCnt < pulseCount && cTime >= sensoryStart;
        if pulseOn
            pulseStart = cTime; %start next pulse
            pulseCnt = pulseCnt + 1;
            trigerAnalog = cStim ~= 1;
            trigerCamera = pulseCnt == 1; %send camera trigger on first pulse
        end
        
        if (cTime - pulseStart) < pulseDur && cTime >= sensoryStart && visualOn
            if useVisNoise
                
                noiseCnt = noiseCnt + 1;
                if noiseCnt > length(noiseTex)
                    noiseCnt = 1;
                end
                Screen('DrawTexture', window, noiseTex(noiseCnt), srcNoise, [0 0 screenXpixels screenYpixels], 0);
                
            else
                % Draw grating texture, rotated by "angle":
                Screen('DrawTexture', window, gratingtex(iCase), srcRect, destRect, BasicVarVals(ismember(BasicVarNames,'StimAngle'),cTrial),[],[],[],[],kPsychUseTextureMatrixForRotation);
                Screen('DrawTexture', window, masktex(mCase), [0 0 visibleSize visibleSize], destRect); %Draw aperture
            end
            Screen('FillRect', window, 255, [0 0 TrigSize TrigSize]); %draw indicator
            
        else
            Screen('FillRect', window, Background); %show background during gap
            Screen('FillRect', window, 0, [0 0 TrigSize TrigSize]); %make indicator black
        end
        
%         %draw frame by frame indicator - dont use this if indicator should show onset of visual stimulus
%         if rem(Cnt,2) == 1 %show white square on even frame count
%             % Screen('FillRect', window,255,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
%             Screen('FillRect', window,255,[0 0 TrigSize TrigSize])
%         else %show black square on uneven frame count
%             % Screen('FillRect', window,0,[screenXpixels-TrigSize screenYpixels-TrigSize screenXpixels screenYpixels])
%             Screen('FillRect', window,0,[0 0 TrigSize TrigSize])
%         end
        
        %show frame
        cTime = Screen('Flip', window, cTime + 0.5 * ifi);
        timeStamps(Cnt) = cTime;
        
        % trigger analog stimuli if needed
        if trigerAnalog
            handles.WavePlayer.play(cStim); %output sensory stimuli
            trigerAnalog = false;
        end
        
        if cTime >= optoStart && ~isempty(optoOut)
            handles.WavePlayer.play(optoOut);
            optoOut = []; %dont trigger again this trial
        end
        
        if trigerCamera && ~isempty(handles.Arduino)
            %send camera trigger from arduino
            fwrite(handles.Arduino, handles.camByte);
            trigerCamera = false;
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
    Screen('Flip', window);
    
    %stop analog output
    handles.WavePlayer.stop;
    
    % stop camera trigger after sequence is over
    fwrite(handles.Arduino, handles.stopCamByte)
    
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
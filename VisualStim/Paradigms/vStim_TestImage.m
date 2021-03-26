function vStim_TestImage(handles)
% Test image. Draws a frame on the screen edges and a cross through the
% center. Hopefully helps to get proper alignment.

%% Isolate all relevant variables from MainGUI
Cnt = 1;
for iVars = 1:size(handles.StaticVariableNames.String,1)
    temp = textscan(handles.StaticVariableNames.String{iVars},'%s%f');
    BasicVarNames{Cnt} = temp{1}{1};
    BasicVarVals(Cnt) = temp{2};
    Cnt = Cnt+1;
    clear temp
end
for iVars = 1:size(handles.FlexibleVariableNames.String,1)
    temp = textscan(handles.FlexibleVariableNames.String{iVars},'%s%f');
    BasicVarNames{Cnt} = temp{1}{1};
    BasicVarVals{Cnt} = str2num(handles.FlexibleVariableNames.String{iVars}(length(temp{1}{1})+1:end));
    Cnt = Cnt+1;
    clear temp
end

CrossWidth = BasicVarVals(ismember(BasicVarNames,'CrossWidth'));

% Open PTB screen window
PsychDefaultSetup(1);
Screen('Preference', 'SkipSyncTests', 1)
screenNumber = max(Screen('Screens')); % Draw to the external screen if avaliable

window = Screen('OpenWindow', screenNumber, 0); %open ptb window and save handle in pSettings
handles.Settings.rRate=Screen('GetFlipInterval', window); %refresh rate
[screenXpixels, screenYpixels] = Screen('WindowSize', window); % Get the size of the current window
handles.ScreenRes.String = [num2str(screenXpixels) ' x ' num2str(screenYpixels)];

% Set priority to ensure high drawing speed
Priority(MaxPriority(window));
[screenXpixels,screenYpixels] = Screen('WindowSize', window); % Get the size of the current window

% get cm coordinates
screenSize = textscan(handles.ScreenSize.String,'%f%c%f'); %get screen size in cm
screenWidth = screenSize{1};
screenHeight = screenSize{3};

%% Draw image
Screen('FrameRect', window, [], [], BasicVarVals(ismember(BasicVarNames,'FrameWdith')));
Screen('FillRect', window,255,[round(screenXpixels/2)-(CrossWidth/2) 0 round(screenXpixels/2)+(CrossWidth/2) screenYpixels]);
Screen('FillRect', window,255,[0 round(screenYpixels/2)-(CrossWidth/2) screenXpixels round(screenYpixels/2)+(CrossWidth/2)]);

%draw horizontal cm lines.
pxPerCM = screenYpixels / screenHeight;
for x = 1 : floor(screenHeight)-1
    Screen('FillRect', window,128,[0 round(pxPerCM*x) screenXpixels round(pxPerCM*x)+5]);
end

pxPerCM = screenXpixels / screenWidth;
for x = 1 : floor(screenWidth)-1
    Screen('FillRect', window,128,[round(pxPerCM*x) 0 round(pxPerCM*x)+5 screenYpixels]);
end

Screen('Flip', window); %
handles.Status.String = 'Testimage - Press any key to release'; 
drawnow();

% wait for paradigm to be released
KbStrokeWait;
handles.LoadParadigm.Value = false;
handles.LoadParadigm.String = 'Load selected paradigm';
handles.Status.String = 'No paradigm loaded';
sca
end
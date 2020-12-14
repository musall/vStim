function [stim, events] = vStim_getStimSequence(pulseDur, analogRate)
% Function to generate stimulus waveforms for the multisensory stimulation
% paradigm.

resFreq = 205; %resonant frequency of LRA stimulators (for tactile stimulus)
buzzLength = S.BuzzDuration;
clickLength = S.ClickDuration;
analogRate = S.sRate; %sampling rate of output module
events = cell(2,2);

%% Generate individual event waveforms
%auditory stimuli (multifrequency convolved click)
tt = -pi:2*pi*1000/(analogRate*clickLength):pi; tt(end) = [];
click = (1+cos(tt)).*(sin(2*tt)+sin(4*tt)+sin(6*tt)+sin(8*tt)+sin(16*tt));
click = click*S.ClickLoudness/max(click);

%stim event for somatosensory stimulator
tt = -pi:2*pi*1000/(analogRate*buzzLength):pi; tt(end) = [];
buzz = (cos(tt)*(S.BuzzStrength)+1)/2;
buzz = (buzz(1,:) - min(buzz)) / max((buzz(1,:) - min(buzz)));
for x = 1 : 2
    events{x} = 0 : 1/S.BuzzRate : pulseDur-1/analogRate; %produce regular sequence if remaining time is too short to allow ITI variability
end

%% create sequences
for x = 1 : 2
    events{x} = 0 : 1/S.ClickRate : pulseDur-1/analogRate; %produce regular sequence if remaining time is too short to allow ITI variability
end
for x = 1 : 2
    events{x+2} = 0 : 1/S.BuzzRate : pulseDur-1/analogRate; %produce regular sequence if remaining time is too short to allow ITI variability
end

stim = zeros(4,round(analogRate*pulseDur)); %preallocate stimulus sequence
for x = 1 : numel(events)
    for y = 1 : length(events{x})
        cEvent = round(events{x}(y)*analogRate);
        if x < 3 %auditory stimulus, x > 2 are visual
            stim(x,cEvent+1:cEvent+size(click,2)) = click;
        elseif x > 2 && x < 5
            stim(x,cEvent+1:cEvent+size(buzz,2)) = buzz;
        end
    end
end

% add resonant frequency to tactile stimulus
cycles = 2*pi*resFreq*(size(stim,2)/analogRate);
tt = 0:cycles/(round(analogRate*(size(stim,2)/analogRate))):cycles; tt(end-1) = [];
stim(3:4,:) = bsxfun(@times,(stim(3:4,:) +0.01),sin(tt));

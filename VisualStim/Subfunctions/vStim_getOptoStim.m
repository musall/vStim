function optoStim = vStim_getOptoStim(sRate, optoDur, optoRamp, optoFreq, onRamp)
% Generate optogenetic stimulus sequences. This creates stimulus sequences
% for blue and red laser light. One is a sine-wave sequence for
% excitation, the other is a ramping stimulus for inhibition.

if ~exist('onRamp','var') || isempty(onRamp)
    onRamp = 0.05; %ramp onset of square wave to avoid light artifacts
end

if optoRamp > optoDur - onRamp
    fprintf('!!! Optoramp(%fs) is longer as maximum pulse duration(%fs). Shortened accordingly !!!\n',optoRamp,optoDur - onRamp);
    optoRamp = optoDur - onRamp;
end

optoStim = zeros(2, round(sRate*optoDur)); %two channels. First is blue, second is red.

%stim sequence for red laser channels (first channel)
% % sine-wave pulses (can be problematic when using laser where output does not start at 0)
% redDur = floor(optoFreq * optoDur) * (1/optoFreq); %make sure this red signal only contains complete cycles
% tt = 0 : 1/sRate : redDur - 1/sRate;
% optoStim(1, 1:length(tt)) = (cos(2*pi*optoFreq*tt-pi)+1) ./2;

%ramping stimulus for blue laser channels (second channel)
rampLength = round(sRate * optoDur);
signal = ones(1, rampLength);
signal(1 : onRamp*sRate) = (1/sRate : 1/sRate : onRamp) ./ onRamp; %change first part to ramp
signal(round(sRate * (optoDur - optoRamp)) : end) = (optoRamp : -1/sRate : 0) ./ optoRamp; %change last part to ramp

% optoStim(2, :) = signal;
optoStim(1:2, :) = repmat(signal,2,1);
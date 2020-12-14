function optoStim = vStim_getOptoStim(sRate, optoDur, optoRamp, optoFreq, redPower, bluePower)
% Generate optogenetic stimulus sequences. This creates stimulus sequences
% for blue and red laser light. One is a square-wave sequence for
% excitation, the other is a ramping stimulus for inhibition.

optoStim = zeros(4, round(sRate*optoDur)); %four channels, first two rows are for laser 1, row 3-4 is for laser 2

%stim sequence for red laser channels (first channel)
tt = 0 : 1/sRate : optoDur - 1/sRate;
signal = (square(2*pi*optoFreq*tt)+1) ./2 .* redPower;
optoStim(1:2,:) = repmat(signal, 2, 1);

%ramping stimulus for blue laser channels (second channel)
signal = ones(1, sRate * optoDur);
signal(sRate * (optoDur - optoRamp) : end) = (optoRamp : -1/sRate : 0) ./ optoRamp; %change end part to ramp
optoStim(3:4,:) = repmat(signal, 2, 1) .* bluePower;
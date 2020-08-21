function [imStack, coords, rawIms] = sparseNoiseFor2PRig(stimP, dispP)
% [imStack, coords, rawIms] = sparseNoiseFor2PRig(stimP, dispP)
% 
% Produce a spherically-transformed sparse noise stimulus for use with
% Matt's 2P rig. Example inputs with explanation given below. Outputs are
% the transformed image stack, coordinates and colors of the squares (see
% sparseNoisePatterns for explanation), and the untransformed image stack
% (mostly useful for debugging).
% 
% Note: using a spherical coordinate tranformation, not a pincushion
% distortion. This means that the coordinates are spherical coordinates
% (azimuth and altitude, as for LabRigger), NOT that distances are
% preserved (as for the Allen).
% 
stimP.squareSize = 4;           % Size of squares, in degrees
stimP.minDist = 10;             % Distance between squares (center to corner) in degrees
stimP.nTrials = 100;            % Number of trials to generate
stimP.bright = [255 0 127];     % Brightness of [white black background]

dispP.resX = 1080;        % output resolution, X (for tilted monitor)
dispP.resY = 1920;        % output resolution, Y (for tilted monitor)
dispP.physWidth  = 26.5;  % width of screen, in cm (for tilted monitor)
dispP.physHeight = 47.5;  % height of screen, in cm (for tilted monitor)
dispP.eyeX = 15;          % eye X location from left, in cm (for tilted monitor)
dispP.eyeY = 29;          % eye Y location from top, in cm (for tilted monitor)
dispP.zDistBottom = 5;   % Distance to bottom of screen, along the horizontal eye line, in cm
dispP.zDistTop    = 12;   % Same thing for top of screen, in cm
dispP.degPerPixel = stimP.squareSize;    % Number of degrees per pixel of the undistorted image (real number)
dispP.bright = [255 0 127];     % Brightness of [white black background]


% Figure out size of screen in degrees and squares
minZ = min(dispP.zDistBottom, dispP.zDistTop);

wInDeg = atand(dispP.eyeX / minZ) + atand((dispP.physWidth - dispP.eyeX) / minZ);
wInSq = ceil(wInDeg / stimP.squareSize);
hInDeg = atand(dispP.eyeY / minZ) + atand((dispP.physHeight - dispP.eyeY) / minZ);
hInSq = ceil(hInDeg / stimP.squareSize);
screenSize = [hInDeg wInDeg];


% Compute where on the image to call the center
% Center this point in the undistorted image on the mouse's eye (degrees)
% Add +1 to account for padding in gridFromCoordsSingle
dispP.cartCX = wInSq * atand(dispP.eyeX / minZ) / wInDeg + 1;
dispP.cartCY = hInSq * atand(dispP.eyeY / minZ) / hInDeg + 1;

% Produce sparse noise pattern
coords = sparseNoisePatterns(screenSize, stimP.squareSize, stimP.minDist, stimP.nTrials);

% Warning: padding of image by a 1-pixel border occurs here (had to do this
% or interpolation gets messed up when we do the spherical distortion)
rawIms = gridFromCoordsSingle(ceil(screenSize ./ stimP.squareSize), coords, stimP.bright);

% Run spherical distortion
imStack = sphereDistortImage(rawIms, dispP);

% Account for monitor tilt
imStack = rot90(imStack, 1);

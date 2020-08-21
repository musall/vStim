function zImStack = sphereDistortImage(imStack, P)
% zImStack = sphereDistortImage(imStack, P)
% 
% Spherically distort an image for display in the imaging rigs. Based on
% code from Labrigger (with some substantial improvements, modifications,
% and bug fixes): http://labrigger.com/blog/2012/03/06/mouse-visual-stim/
% 
% imStack should be a height x width x nImages matrix. It does not have to
% be at the desired full resolution, but it should have the correct aspect
% ratio. nImages may be 1.
% 
% P should have all of the following fields (filled in here as an example):
% 
% P.resX = 1920;        % output resolution, X
% P.resY = 1080;        % output resolution, Y
% P.physWidth  = 53;    % width of screen, in cm
% P.physHeight = 30;    % height of screen, in cm
% P.eyeX = 26.50;       % eye X location, in cm
% P.eyeY = 11.42;       % eye Y location, in cm
% P.zDistBottom = 24;   % Distance to bottom of screen, along the horizontal eye line, in cm
% P.zDistTop    = 14;   % Same thing for top of screen, in cm
% P.degPerPixel = 4;    % Number of degrees per pixel of the undistorted image (real number)
% P.cartCX = 10;        % Center this point in the undistorted image on the mouse's eye (X, degrees)
% P.cartCY = 5;         % Same for Y
% 
% zImStack will be a uint8 array of size:
% P.resY x P.resY x size(imStack, 3)
% 
% For a test run:
% checkSize = 4; % pixels per side of each check
% w = 192; % width, in pixels
% h = 108; % height, in pixels
% im = double(checkerboard(checkSize,round(h/checkSize), round(w/checkSize)) > 0.5);
% zIm = sphereDistortImageLabrigger(im, P);
% figure; hold on;
% subplot(2, 1, 1); imshow(im);
% subplot(2, 1, 2); imshow(zIm);


% Background for interpolation. Using 127 will make any unfilled areas
% blend in; using other values will make it obvious if you're not using the
% whole screen area
% gray = P.bright(3);
gray = 100;


%% Compute useful coordinates

top = P.physHeight - P.eyeY;
% bottom = -dispP.eyeY;
% right = dispP.eyeX;
% left = dispP.eyeX - dispP.physWidth;
left = -P.eyeX;
% right = P.physWidth - P.eyeX;


%% Project monitor pixels into spherical coordinates
% In image space, x and y are width and height of monitor and z is the
% distance from the eye. I want Theta to correspond to azimuth and Phi to
% correspond to elevation, but these are measured from the x-axis and x-y
% plane, respectively. So I need to exchange the axes this way, prior to
% converting to spherical coordinates:
% orig (image) -> for conversion to spherical coords
% Z -> X
% X -> Y
% Y -> Z

% Theta is reversed from the convention we want (left low) so we flip xi

[xi, yi] = meshgrid(linspace(1, 0, P.resX), linspace(0, 1, P.resY));
cartPointsXFull = left + P.physWidth .* xi;
cartPointsYFull = top - P.physHeight .* yi;
cartPointsZFull = P.zDistTop + (P.zDistBottom - P.zDistTop) .* yi;

[sphrPointsTh, sphrPointsPh] = cart2sph(cartPointsZFull, cartPointsXFull, cartPointsYFull);


%% Recenter image pixels in spherical coordinates
[initResY, initResX, nIm] = size(imStack);

[xi, yi] = meshgrid(0:initResX-1, 0:initResY-1);
imX = pi / 180 * P.degPerPixel * (xi - P.cartCX + 0.5);
imY = pi / 180 * P.degPerPixel * (yi - P.cartCY + 0.5);

%% Use interp2 to find nearest image pixel for each monitor pixel

zImStack = repmat(uint8(gray), P.resY, P.resX, nIm);
for fr = 1:nIm
  zImStack(:, :, fr) = interp2(imX, imY, imStack(:, :, fr), sphrPointsTh, sphrPointsPh, 'nearest', gray);
end


% Here's the line of code to use for the reverse transformation.
% To see what the visual stimulus would look like from the mouse's point-of-view (MPOV) if it were not corrected.
% (From LabRigger, needs to be fixed to work with this code)
% 
% ZI_origMPOV = griddata(sphrPointsTh,sphrPointsPh,im,cartPointsX.*fx,cartPointsY.*fy);

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

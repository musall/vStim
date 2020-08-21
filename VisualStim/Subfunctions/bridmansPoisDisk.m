function coords = bridmansPoisDisk(siz, minDist, discrete)
% coords = bridmansPoisDisk(siz, minDist [, discrete])
% 
% Implement Bridman's Poisson disk sampling algorithm (also called
% blue-noise sampling):
% http://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
% 
% This algorithm lets you produce random samples of points with nearly
% uniform density, within the rectangle of size siz, that are fairly
% densely packed and are no closer than minDist to one another.
% 
% siz is the size of the bounding rectangle, minDist is the minimum
% distance between points, and coords is an nPts x 2 array of coordinates.
% If discrete is true, all coordinates will be integers. Default false.
% 
% The algorithm basically works like this. Initialize with a random point
% in the rectangle. Add this point to the "active" list. For every
% iteration, first pick a point from the active list. Now choose a point
% randomly from the annulus from minDist to 2*minDist from the active
% point. Check if it's too close to any other points. If it's ok, keep it.
% If it's too close, try a different point in the annulus. If we fail 30
% times in a row, remove this point from the active list. Stop when the
% active list is empty.
% 
% To make the distance check fast, we use a background grid (which can hold
% at most one point per square if we obey the minimum distance rule). The
% resulting algorithm runs in linear time and memory with search area.

% Parameters
maxAttempts = 30;

% Optional argument
if ~exist('discrete', 'var')
  discrete = 0;
end

% Need to make the background grid boxes this size because two points could
% be at opposite corners of the box
gridSize = minDist / (sqrt(2) + 1e-8);

minDist2 = minDist ^ 2;

% Background grid -- each square in the grid can contain at most one point.
% This will hold the index of the point in that square (if a point is
% present in the square)
grd = zeros(ceil(siz / gridSize));

% Pre-allocate a coordinates list. We'll trim it at the end.
coords = NaN(numel(grd), 2);

% Pre-allocate the active list
active = NaN(numel(grd), 1);

% Grid boxes to test
% Note that it's 5 wide because the width of a grid box is sqrt(2) smaller
% than minDist
[gAn1, gAn2] = meshgrid(1:5, 1:5);
gAnCoords = [gAn1(:), gAn2(:)] - 3;  % center is now (0, 0)

% Choose initial point using uniform distribution
nCoords = 1;
nActive = 1;
if discrete
  coords(1, :) = [randi(siz(1), 1), randi(siz(2), 1)];
else
  coords(1, :) = [siz(1) * rand, siz(2) * rand];
end
active(1, :) = 1;
grd(ceil(coords(1, 1) / gridSize), ceil(coords(1, 2) / gridSize)) = 1;

% Main loop of algorithm
while nActive > 0
  % Choose an active point
  acti = randi(nActive, 1);
  
  % Get the active point's coordinates
  ac = coords(active(acti), :);
  
  % Try points in the annulus until success or failure
  found = 0;
  for i = 1:maxAttempts
    % Choose a point to try
    % Note: the distance must be chosen taking into account that the
    % perimeter of a circle increases with radius. sqrt is the right
    % correction to get a flat density. The "3" is because we want
    % r2^2-r1^2, r1 is minDist, and r2 is 2*minDist.
    theta = 2 * pi * rand;
    dist = sqrt(rand * 3 * minDist2 + minDist2);
    if discrete
      c = ac + round(dist * [cos(theta) sin(theta)]);
    else
      c = ac + dist * [cos(theta) sin(theta)];
    end
    
    % Can't try points outside the bounds
    if any(c < 1) || any(c > siz)
      continue;
    end
    
    % Coordinates in the grid
    g = ceil(c / gridSize);
    cg = bsxfun(@plus, g, gAnCoords);
    cg = cg(all(cg > 0 & bsxfun(@le, cg, size(grd)), 2), :); % don't run off the edge
    cg = sub2ind(size(grd), cg(:, 1), cg(:, 2));
    
    % Test against possible conflicting points, using the background grid
    % to dramatically reduce the number of points we need to test in the
    % worst case, large numbers of points (to 25 points max)
    inds = grd(cg);
    inds = inds(inds > 0);
    if ~isempty(inds)
      dists2 = sum(bsxfun(@minus, c, coords(inds, :)) .^ 2, 2);
      if any(dists2 < minDist2)
        continue;
      end
    end
    
    % If you're debugging, can uncomment the below to perform an exhaustive
    % search
%     dists = sqrt(sum((coords(1:nCoords, :) - c) .^ 2, 2));
%     if min(dists) < minDist
%       error('Found points too close!');
%     end
    
    % If we get here, the candidate point is good
    nCoords = nCoords + 1;
    coords(nCoords, :) = c;
    nActive = nActive + 1;
    active(nActive) = nCoords;
    grd(g(1), g(2)) = nCoords;
    found = 1;
    break;
  end
  
  % If we failed this cycle, remove the active point from the active list
  % (It's fine to alter the order of the active list)
  if ~found
    active(acti) = active(nActive);
    nActive = nActive - 1;
  end
end

% Trim off unused part of coords array
coords = coords(1:nCoords, :);

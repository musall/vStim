function coords = sparseNoisePatterns(screenSize, squareSize, minDist, nTrials)
% Produces sparse noise (for visual stimuli) using a discrete version of
% Bridman's algorithm.
% 
% Result is a cell array with one element of the cell array per trial. Each
% cell contains an nSquares x 3 matrix, where each row is Y, X, white|black.
% white|black is 1 or 0.
% 
% screenSize is resolution in degrees as [Y X] (e.g., 40 x 60).
% squareSize is the width of the squares in degrees (e.g., 4).
% minDist is the minimum distance between squares (square center to nearest
%   square corner) as a multiple of square width (e.g., if a square spans 4
%   degrees and you want a minimum distance of 20 degrees, use 5).
% nTrials is the number of trials to generate.


%% Parameters

% Whether to balance white and black for each pixel.
balance = 1;

% Whether to use degrees or squares as pixels
degreeCoords = 0;

% How much we'll need to add to the grid coordinates in order to get a
% value in pixels
marg = [(mod(screenSize(1), squareSize) + squareSize) / 2, ...
  (mod(screenSize(2), squareSize) + squareSize) / 2];


%% Get coordinates in grid units

coords = cell(1, nTrials);
for tr = 1:nTrials
  % Need to add sqrt(2)/2 to minDist, because the squares are not points and
  % we measure from the center of a square to the corner of the nearest
  % square
  coords{tr} = bridmansPoisDisk(floor(screenSize / squareSize), ...
    floor(minDist / squareSize) + sqrt(2) / 2, 1);
end


%% Add black (0) or white (1) to each coordinate pair as a third column

if balance
  % Compute number of samples in each square
  grid = zeros(floor(screenSize / squareSize));
  for tr = 1:nTrials
    for p = 1:size(coords{tr}, 1)
      grid(coords{tr}(p, 1), coords{tr}(p, 2)) = grid(coords{tr}(p, 1), coords{tr}(p, 2)) + 1;
    end
    % Pre-allocate for third column
    coords{tr} = [coords{tr}, zeros(size(coords{tr}, 1), 1)];
  end
  
  % Create a cell array where each element holds the randomized array of
  % 1's and 0's for the corresponding grid element.
  colors = cell(size(grid));
  for g = 1:numel(grid)
    nW = ceil(grid(g) / 2);
    if rand > 0.5
      colors{g} = [ones(1, nW), zeros(1, grid(g) - nW)];
    else
      colors{g} = [zeros(1, nW), ones(1, grid(g) - nW)];
    end
    colors{g} = colors{g}(randperm(length(colors{g})));
  end
  
  % Dole out 1's and 0's
  onEl = ones(size(grid));
  for tr = 1:nTrials
    for p = 1:size(coords{tr}, 1)
      el = onEl(coords{tr}(p, 1), coords{tr}(p, 2));
      coords{tr}(p, 3) = colors{coords{tr}(p, 1), coords{tr}(p, 2)}(el);
      onEl(coords{tr}(p, 1), coords{tr}(p, 2)) = onEl(coords{tr}(p, 1), coords{tr}(p, 2)) + 1;
    end
  end
  
else
  for tr = 1:nTrials
    coords{tr} = [coords{tr}, randi(2, size(coords{tr}, 1), 1) - 1];
  end
end


%% Convert grid coordinates into degree coordinates

if degreeCoords
  for tr = 1:nTrials
    coords{tr} = bsxfun(@minus, coords{tr}, [1 1 0]);
    coords{tr} = bsxfun(@times, coords{tr}, [squareSize squareSize 1]);
    coords{tr} = floor(bsxfun(@plus, coords{tr}, [marg 0]));
  end
end


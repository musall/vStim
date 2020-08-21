function grd = gridFromCoordsSingle(siz, coords, bright)
% grd = gridFromCoordsSingle(siz, coords)
% 
% Take the output of sparseNoisePatterns and convert it into an image.
% bright should specify levels for white, black, and background (e.g., [255
% 0 127]). Values will be type single.

black = bright(2);
white = bright(1) - black;
back = bright(3);
grd = repmat(single(back), [siz+2 length(coords)]);

for tr = 1:length(coords)
  for c = 1:size(coords{tr}, 1)
    grd(coords{tr}(c, 1)+1, coords{tr}(c, 2)+1, tr) = coords{tr}(c, 3) * white + black;
  end
end

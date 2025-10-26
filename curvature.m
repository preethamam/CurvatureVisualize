function [shape, Icurv] = curvature(image, boundaryPoint, curvatureThresh, ...
                                    bp_tangent, interp_resolution, loopclose)

%%***********************************************************************%
%*                         Curvature measure                            *%
%*                  Measure shape properties of loops.                  *%
%*                                                                      *%
%* Original author: Dr. Meghan Driscoll                                 *%
%* Modified by: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam
%* Date: 08/02/2021                                                     *%
%************************************************************************%
%
%************************************************************************%
%
% Usage: [shape, Icurv] = curvature(image, boundaryPoint, curvatureThresh, ...
%                                     bp_tangent, interp_resolution, loopclose)
% Inputs: % Inputs:
%           Image                - input image
%           boundaryPoint        - number of boundary points over which curvature is found 
%           curvatureThresh      - the largest curvature magnitude allowed in the cutoff curvature
%           bp_tangent           - the number of boundary points over which the boundary tangent 
%                                  direction is measured
%           interp_resolution    - interpolation resolution -- the minimum number of pixels seperating 
%                                  boundary points after interpolation
%           loopclose            - 0 - if open boundaries | 1 - if closed boundaries
% 
% Outputs: 
%           shape                
%           .curvature          - the boundary curvature at each boundary point (uses snakeNum) 
%                                 Curvatures above or below a cutoff are given the magnitude of the cutoff
%           .meanNegCurvature   - the mean negative curvature
%           .numIndents         - the number of boundary regions over which the curvature is negative
%           .tangentAngle       - the angle of the tangent to the boundary at each boundary point
%           .tortuosity         - the boundary tortuousity (a measure of how bendy the boundary is)
%           .uncutCurvature     - the uncut boundary curvature at each boundary point (uses snakeNum)

%           Icurv               - Output image (padded image to make 3 channel. This fixes the plot)
%--------------------------------------------------------------------------
% 
%%%%%%%%%%%%%%%% MEASURE SHAPE %%%%%%%%%%%%%%%%
% Original author: Dr. Meghan Driscoll
% Modified and compacted/concised a complicated codebase by: Preetham Manjunatha
% 
% Thanks to Dr. Meghan Driscoll who kindly shared her code for academic purpose.
% If you use this code for visualization and other academic/research/any purposes. 
% 
% Please cite:
% 
% Reference:
% Driscoll MK, McCann C, Kopace R, Homan T, Fourkas JT, Parent C, et al. (2012) 
% Cell Shape Dynamics: From Waves to Migration. 
% PLoS Comput Biol 8(3): e1002392. 
% https://doi.org/10.1371/journal.pcbi.1002392
% 
% Important note: This code initially used a parfor to speed up the things. 
% Now it has been vectorized for fast calculations.

% RGB to binarization
if(size(image,3) == 3)
    Igray = rgb2gray(image);
    binaryimage = imbinarize(Igray);
    Icurv = image;
elseif (islogical(image))
    binaryimage = image;      
    Icurv = im2uint8(repmat(image,1,1,3));
else
    binaryimage = imbinarize(I);
    Icurv = im2uint8(repmat(image,1,1,3));
end

% Find boundaries X and Y coordinates
cc_index = 1;
boundaries = bwboundaries(binaryimage,8);
x = boundaries{cc_index}(:, 2);
y = boundaries{cc_index}(:, 1);

% Perimeter of the binary component
stats = regionprops(binaryimage,'perimeter');
perimeter = cat(1,stats(cc_index).Perimeter);

% Interpolate for more points
% Every coordinates
xn = x+rand(size(x))*1e-8;
yn = y+rand(size(x))*1e-8;
[xi,yi] = snakeinterp1(xn,yn,interp_resolution);

% Every nth coordinates
switch loopclose
    case 0
        xn = xi(1:boundaryPoint:end);
        yn = yi(1:boundaryPoint:end);
    case 1
        xn = [xi(1:boundaryPoint:end), xi(1)];
        yn = [yi(1:boundaryPoint:end), yi(1)];
end
shape_XY = [xn;yn]';
M = numel(xn);

% initialize variables    
shape_curvature        = zeros(1, M);
shape_uncutCurvature   = zeros(1, M);
shape_meanNegCurvature  = NaN(1,1);
shape_numIndents = NaN(1,1);
shape_tortuosity = NaN(1,1);
shape_tangentAngle = NaN(1,M);

% calculate the curvature (by finding the radius of the osculating circle using three neaby boundary points)
bp_positions = [shape_XY(M-1-boundaryPoint:M-1,:); shape_XY(1:M-1,:); shape_XY(1:boundaryPoint+1,:)];

%% Vectorized curvature calculation start
%--------------------------------------------------------------------------
% Indices for triplets
idx = (1:M).';
i1  = idx;
i2  = idx + boundaryPoint;
i3  = idx + 2*boundaryPoint;

% Gather points
P1 = bp_positions(i1, :);   % Mx2
P2 = bp_positions(i2, :);
P3 = bp_positions(i3, :);

% Initial slopes
s12 = (P1(:,2) - P2(:,2)) ./ (P1(:,1) - P2(:,1));
s23 = (P2(:,2) - P3(:,2)) ./ (P2(:,1) - P3(:,1));

% ---- First swap: if slope12 is Inf/-Inf or 0, swap point2 <-> point3
mask1 = ~isfinite(s12) | (s12 == 0);
if any(mask1)
    tmp           = P2(mask1, :);
    P2(mask1, :)  = P3(mask1, :);
    P3(mask1, :)  = tmp;
end

% Recompute slopes after first swap
s12 = (P1(:,2) - P2(:,2)) ./ (P1(:,1) - P2(:,1));
s23 = (P2(:,2) - P3(:,2)) ./ (P2(:,1) - P3(:,1));

% ---- Second swap: if slope23 is Inf/-Inf, swap point1 <-> point2
mask2 = ~isfinite(s23);
if any(mask2)
    tmp           = P1(mask2, :);
    P1(mask2, :)  = P2(mask2, :);
    P2(mask2, :)  = tmp;
end

% Recompute slopes after second swap
s12 = (P1(:,2) - P2(:,2)) ./ (P1(:,1) - P2(:,1));
s23 = (P2(:,2) - P3(:,2)) ./ (P2(:,1) - P3(:,1));

% Flat boundary where slope12 == slope23  --> curvature = 0 (already zeroed)
isFlat = (s12 == s23);

% Curved cases to compute
curved = ~isFlat;
if any(curved)
    % Slice arrays for curved rows
    P1c  = P1(curved, :);
    P2c  = P2(curved, :);
    P3c  = P3(curved, :);
    s12c = s12(curved);
    s23c = s23(curved);

    % Centers and midpoints
    % x_center
    xc = ( s12c .* s23c .* (P1c(:,2) - P3c(:,2)) ...
         + s23c .* (P1c(:,1) + P2c(:,1)) ...
         - s12c .* (P2c(:,1) + P3c(:,1)) ) ...
         ./ (2 .* (s23c - s12c));

    mid12 = 0.5 .* (P1c + P2c);
    mid13 = 0.5 .* (P1c + P3c);

    % y_center = (-1/slope12) * (x_center - mid12_x) + mid12_y
    yc = (-1 ./ s12c) .* (xc - mid12(:,1)) + mid12(:,2);

    % Radius and uncut curvature
    r   = sqrt( (P1c(:,1) - xc).^2 + (P1c(:,2) - yc).^2 );
    uncut = 1 ./ r;

    % Start with uncut, then clamp for visualization
    curv = uncut;
    big  = curv > curvatureThresh;
    curv(big) = curvatureThresh;

    % Sign via inpolygon(midpoint13)
    [In, On] = inpolygon(mid13(:,1), mid13(:,2), shape_XY(:,1), shape_XY(:,2));

    % If not inside, negate
    notIn = ~In;
    if any(notIn)
        curv(notIn)  = -curv(notIn);
        uncut(notIn) = -uncut(notIn);
    end

    % If On or uncut is non-finite, set to zero
    bad = On | ~isfinite(uncut);
    if any(bad)
        curv(bad)  = 0;
        uncut(bad) = 0;
    end

    % Write back to outputs
    shape_curvature(1, curved)      = curv.';
    shape_uncutCurvature(1, curved) = uncut.';
end

% Vectorized version end
%--------------------------------------------------------------------------

% find the mean negative curvature (really this should use a constant dist snake)
listCurve = shape_uncutCurvature(1,1:M-1);
listNegCurve = abs(listCurve(listCurve < 0));
if ~isempty(listNegCurve) 
    shape_meanNegCurvature(1,1) = sum(listNegCurve)/(M-1);
else
    shape_meanNegCurvature(1,1) = 0;
end

% find the number of negative boundary curvature regions
curveMask = (listCurve < 0);
curveMaskLabeled = bwlabel(curveMask);
numIndents = max(curveMaskLabeled);
if curveMask(1) && curveMask(end)
    numIndents  = numIndents-1;
end
shape_numIndents(1,1) = numIndents;

% find the tortuosity (should correct units)
shape_tortuosity(1,1) = sum(gradient(shape_uncutCurvature(1,1:M-1)).^2)/perimeter;

% calculate the direction of the tangent to the boundary 
bp_positions_tangent=[shape_XY(M-1-bp_tangent:M-1,:); shape_XY(1:M-1,:); shape_XY(1:bp_tangent+1,:)];
for j=1:M
    point1 = bp_positions_tangent(j,:);
    point2 = bp_positions_tangent(j+2*bp_tangent,:); 
    shape_tangentAngle(1,j) = mod(atan2(point1(1,2)-point2(1,2), point1(1,1)-point2(1,1)), pi);
end

shape.XY = shape_XY;
shape.curvature         = shape_curvature;
shape.uncutCurvature    = shape_uncutCurvature;
shape.meanNegCurvature  = shape_meanNegCurvature;
shape.numIndents = shape_numIndents;
shape.tortuosity = shape_tortuosity;
shape.tangentAngle = shape_tangentAngle;

end


%% Auxillary fucntions
%--------------------------------------------------------------------------
function [xi,yi] = snakeinterp1(x,y,RES)
% SNAKEINTERP1  Interpolate the snake to have equal distance RES
%     [xi,yi] = snakeinterp(x,y,RES)
%
%     RES: resolution desired

%     update on snakeinterp after finding a bug

%      Chenyang Xu and Jerry L. Prince, 3/8/96, 6/17/97
%      Copyright (c) 1996-97 by Chenyang Xu and Jerry L. Prince
%      image Analysis and Communications Lab, Johns Hopkins University

% convert to column vector
x = x(:); y = y(:);

N = length(x);  

% make it a circular list since we are dealing with closed contour
x = [x;x(1)];
y = [y;y(1)];

dx = x(2:N+1)- x(1:N);
dy = y(2:N+1)- y(1:N);
d = sqrt(dx.*dx+dy.*dy);  % compute the distance from previous node for point 2:N+1

d = [0;d];   % point 1 to point 1 is 0 

% now compute the arc length of all the points to point 1
% we use matrix multiply to achieve summing 
M = length(d);
d = (d'*uppertri(M,M))';

% now ready to reparametrize the closed curve in terms of arc length
maxd = d(M);

if (maxd/RES<3)
   error('RES too big compare to the length of original curve');
end

di = 0:RES:maxd;

xi = interp1(d,x,di);
yi = interp1(d,y,di);

N = length(xi);

if (maxd - di(length(di)) <RES/2)  % deal with end boundary condition
   xi = xi(1:N-1);
   yi = yi(1:N-1);
end
end


function q = uppertri(M,N)                      %added by Ilya
% UPPERTRI   Upper triagonal matrix 
%            UPPER(M,N) is a M-by-N triagonal matrix

[J,I] = meshgrid(1:M,1:N);
q = (J>=I);
end
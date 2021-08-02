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
%           bp_tangent           - the number of boundary points over which the boundary tangent direction is measured
%           savePath             - the directory where the variable shape is saved
%           loopclose            - 0 - if open boundaries | 1 - if closed boundaries
% 
% Outputs: 
%           shape                
%           .curvature          - the boundary curvature at each boundary point (uses snakeNum) Curvatures above or below a cutoff are given the magnitude of the cutoff
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
% Important note: This code uses parfor to speed up the things. If you do not have 
% the Matlab parallel computing toolbox. Please make 'parfor' as 'for' in this
% function.
%
% %%%%%%%% This code is way too slow! (curvature should not be in a for loop) %%%%%%%%%
% If I have time I will try to improve this. If anyone improves it, please
% share the modified version code with me.

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
shape_curvature         = NaN(1,M);
shape_uncutCurvature    = NaN(1,M);
shape_meanNegCurvature  = NaN(1,1);
shape_numIndents = NaN(1,1);
shape_tortuosity = NaN(1,1);
shape_tangentAngle = NaN(1,M);

% calculate the curvature (by finding the radius of the osculating circle using three neaby boundary points)
bp_positions = [shape_XY(M-1-boundaryPoint:M-1,:); shape_XY(1:M-1,:); shape_XY(1:boundaryPoint+1,:)];
parfor j = 1:M

    % assign the three points that the circle will be fit to such that the slopes are not infinite 
    point1 = bp_positions(j,:);
    point2 = bp_positions(j+boundaryPoint,:);
    point3 = bp_positions(j+2*boundaryPoint,:); 
    slope12 = (point1(1,2)-point2(1,2))/(point1(1,1)-point2(1,1));
    slope23 = (point2(1,2)-point3(1,2))/(point2(1,1)-point3(1,1));

    if slope12==Inf || slope12==-Inf || slope12 == 0
        point0 = point2; point2 = point3; point3 = point0;
        slope12 = (point1(1,2)-point2(1,2))/(point1(1,1)-point2(1,1));
        slope23 = (point2(1,2)-point3(1,2))/(point2(1,1)-point3(1,1));    
    end

    if slope23==Inf || slope23==-Inf
        point0 = point1; point1 = point2; point2 = point0;
        slope12 = (point1(1,2)-point2(1,2))/(point1(1,1)-point2(1,1));
        slope23 = (point2(1,2)-point3(1,2))/(point2(1,1)-point3(1,1));    
    end

    % if the boundary is flat
    if slope12==slope23  
        shape_curvature(1,j) = 0;

    % if the boundary is curved
    else

        % calculate the curvature
        x_center = (slope12*slope23*(point1(1,2)-point3(1,2))+slope23*(point1(1,1)+point2(1,1))...
                   -slope12*(point2(1,1)+point3(1,1)))/(2*(slope23-slope12));
        midpoint12 = (point1+point2)/2;
        midpoint13 = (point1+point3)/2;
        y_center = (-1/slope12)*(x_center-midpoint12(1,1))+midpoint12(1,2);
        shape_uncutCurvature(1,j) = 1/sqrt((point1(1,1)-x_center)^2+(point1(1,2)-y_center)^2);

        % cutoff the curvature (for visualization)
        shape_curvature(1,j) = shape_uncutCurvature(1,j);
        if shape_curvature(1,j) > curvatureThresh
            shape_curvature(1,j) = curvatureThresh;
        end

        % determine if the curvature is positive or negative
        [In, On] = inpolygon(midpoint13(1,1),midpoint13(1,2),shape_XY(:,1),shape_XY(:,2)); 

        if ~In              
            shape_curvature(1,j) = -1*shape_curvature(1,j);
            shape_uncutCurvature(1,j) = -1*shape_uncutCurvature(1,j);
        end

        if On || ~isfinite(shape_uncutCurvature(1,j))
            shape_curvature(1,j) = 0;
            shape_uncutCurvature(1,j) = 0;
        end

    end 

end

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

dx = x([2:N+1])- x(1:N);
dy = y([2:N+1])- y(1:N);
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



clear; close all; clc;
Start = tic;

%% Inputs
%--------------------------------------------------------------------------
% Image files
I = imread('T.png');    %'2019_04_19_19_33_26_456_Ang_0.bmp');

% Measure shape parameters
boundaryPoint = 10;         %12, 10 number of boundary points curvature is found over 
curvatureThresh = 0.25;     %0.06 the maximum allowed value of the curvature measure
bp_tangent = 10;            % number of boundary points the tangent angle is found over 
interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
loopclose = 1;              % 0 - if open boundaries | 1 - if closed boundaries


%% Find the curvature
[shape_details, Icurv] = curvature(I, boundaryPoint, curvatureThresh, bp_tangent, ...
                                 interpdmin, loopclose);
        
%% Plot curvature
X = shape_details.XY(:,1);
Y = shape_details.XY(:,2);
Z = zeros(size(X));
C = shape_details.curvature'*1;

% Plot
figure;
s1 = subplot(1,2,1); imshow(Icurv)
s2 = subplot(1,2,2); 
imshow(Icurv)
hold on
surf([X(:) X(:)], [Y(:) Y(:)], [Z(:) Z(:)], [C C], ...  % Reshape and replicate data
 'FaceColor', 'none', ...    % Don't bother filling faces with color
 'EdgeColor', 'interp', ...  % Use interpolated color for edges
 'LineWidth', 3);            % Make a thicker line
hold off

cmap = jet;
colormap(cmap);
cb = colorbar;  % Add a colorbar
cb.Label.String = 'Curvature';

s1Pos = get(s1,'position');
s2Pos = get(s2,'position');
s2Pos(3:4) = [s1Pos(3:4)];
set(s2,'position',s2Pos);




%% End parameters
%--------------------------------------------------------------------------
Runtime = toc(Start);

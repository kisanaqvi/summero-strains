function [x_rotated, y_rotated] = drawEllipse(p, majorAxes, minorAxes, centroid_X, centroid_Y, angles, conversionFactor)

% adapted from drawEllipse.m, which was made for videos. this script is slightly
% edited to work with single frames.

% last edited: jen, 2021 April 2
% last commit: first commit


% 0. isolate data from particles in current image
% 1. load major axes, min axes, centroids, angles, conversion factor
% 2. for all particles,
%           3. calculate points on ellipse centered at centroid
%           4. rotate ellipse origin, then bring back to centroid
%                   i. define center of rotation
%                  ii. create matrix by which to translate ellipse to and from origin
%                 iii. rotate ellipse around origin with object centered at origin
%                  iv. translate ellipse back to original centroid
% 5. overlay onto current image

a = majorAxes(p)/2; % horizontal radius
b = minorAxes(p)/2; % vertical radius

t = -pi:0.01:pi;

outline1 = (centroid_X+a*cos(t))/conversionFactor;
outline2 = (centroid_Y+b*sin(t))/conversionFactor;


% define the x- and y-data for the original line we would like to rotate
x = outline1;
y = outline2;

% create a matrix of these points, which will be useful in future calculations
v = [x;y];

% choose a point which will be the center of rotation
x_center = centroid_X/conversionFactor;
y_center = centroid_Y/conversionFactor;

% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(x));

% define a 60 degree counter-clockwise rotation matrix
theta = -pi*angles(p)/(180);       % pi/3 radians = 60 degrees
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

%
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin

so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation

% this can be done in one line as:
% vo = R*(v - center) + center

% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);


end
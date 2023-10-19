function R = rotamatz(angle, units)
% ROTAMATZ - create rotation matrix for a rotation of 'angle' about the z axis.	
% R = rotamatz(angle, units)
%   
%   INPUTS:
%   angle  = angle in degrees or radians. 
%   units  = 'deg' or 'rad'.
%   OUTPUTS:
%  	R      = 4x4 rotation matrix.

%   Created by: Matthew R. Walker
%   <mrwalker@shrinenet.org>
%   Last modified: April 18, 2003


if ~exist('units', 'var'), 
    disp('You must define angle units.')
    return
end

if isequal(units, 'deg'),
    angle = (angle*pi)/180; % degs to rads.
end

cosa = cos(angle);
sina = sin(angle);
R = [cosa -sina 0 0; ...
    sina cosa 0 0; ...
    0 0 1 0; ...
    0 0 0 1];	
%--------------------------------------------------------------------------
function out = unit(in)
% UNIT - create a unit vector aligned with input vector.

% 	Created by: Matthew Walker
%	<mrwalker@shrinenet.org>
%   Last modified: June 11, 2004


out = in/sqrt(sum(diag(in'*in)));
%--------------------------------------------------------------------------
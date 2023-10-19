function ACS=FemurACS(pathname,fullivfile,shaftivfile,cylinder_fit_axis,cylinder_fit_base_pt,cylinder_fit_height)
% ACS=FEMURACS(PATHNAME,FULLIVFILE,SHAFTIVFILE,CYLINDER_FIT_AXIS,CYLINDER_FIT_BASE_PT,CYLINDER_FIT_HEIGHT)
%   PATHNAME: path name
%   FULLIVFILE: full femur iv file
%   SHAFTIVFILE: cropped shaft iv file
%   CYLINDER_FIT_AXIS: axis of cylinder fit to cropped condyles
%   CYLINDER_FIT_BASE_PT: base of cylinder fit to cropped condyles
%   CYLINDER_FIT_HEIGHT: height of cylinder fit to cropped condyles

%% calculate mass properties of iv femur model
[pts_CT conn_CT] = read_vrml_fast(fullfile(pathname,fullivfile));
conn_CT(:,4) = [];
conn_CT(:) = conn_CT(:)+1;

% determine full model center of mass as well as the full model inertial axes
[full_centroid,full_surface_area,Volume,full_eigenvalues,full_eigenvectors,full_I1,full_I2,full_I_CoM,full_I_origin,full_patches] = mass_properties(pts_CT,conn_CT);

%% calculate mass properties of iv shaft model
[pts_shaft conn_shaft] = read_vrml_fast(fullfile(pathname,shaftivfile));
conn_shaft(:,4) = [];
conn_shaft(:) = conn_shaft(:)+1;

% determine shaft model center of mass as well as the shaft model inertial axes
[shaft_centroid,shaft_surface_area,Volume,shaft_eigenvalues,shaft_eigenvectors,shaft_I1,shaft_I2,shaft_I_CoM,shaft_I_origin,shaft_patches] = mass_properties(pts_shaft,conn_shaft);

%% determine center of ACS from the midpoint of the cylinder fit to the condyles
T=cylinder_fit_base_pt+(cylinder_fit_height/2)*unit(cylinder_fit_axis);

%% assign shaft vector as the smallest inertial axis of the shaft (long axis)
shaft_vector=shaft_eigenvectors(:,1);

% make sure shaft vector is pointing toward the proximal femur
correct_direction=unit(shaft_centroid-full_centroid); % determine vector between centroid of shaft and centroid of entire femur bone to compare to shaft vector
angular_difference1=AngleDiff(correct_direction,shaft_vector); % comparison of shaft vector and correct direction vector

% if the shaft vector and direction are facing in opposite directions, then invert the shaft vector
if angular_difference1 > 90
    shaft_vector= -shaft_vector;
end

%% create rotation matrix from cylinder fit vector and shaft vector

% make sure that the medial lateral axis is pointing in a direction that
% will allow long axis to be pointing proximal and anterior posterior axis
% pointing posterior.  if this is not the case negate the cylinder fit axis
angular_difference2=AngleDiff(correct_direction,unit(cross(unit(shaft_vector),unit(cylinder_fit_axis))));
if angular_difference2>90
    cylinder_fit_axis=-cylinder_fit_axis;
end

% assign x-axis as the medial lateral axis based on the cylinder fit vector
% assign y-axis as the anterior posterior axis based on the cross between
% the shaft vector and the cylinder fit vector
% assign z-axis as the long axis based on the cross between the medial
% lateral axis and the anterior posterior axis
R=eye(3);
R(:,1)=unit(unit(cylinder_fit_axis));
R(:,2)=unit(cross(unit(shaft_vector),R(:,1)));
R(:,3)=unit(cross(R(:,1),R(:,2)));

%% compile ACS as an RT from the R and T
ACS=[R,T';[0 0 0 1]];
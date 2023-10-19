function ACS=TibiaACS(pathname,fullivfile,plateauivfile,anteriorpt)
% ACS=TIBIAACS(PATHNAME,FULLIVFILE,PLATEAUIVFILE,ANTERIORPT)
%   PATHNAME: path name
%   FULLIVFILE: full tibia iv file
%   PLATEAUIVFILE: cropped plateau iv file
%   ANTERIORPT: 3D point toward anterior side of tibia

%% calculate mass properties of iv tibia model
[pts_full conn_full] = read_vrml_fast(fullfile(pathname,fullivfile));
conn_full(:,4) = [];
conn_full(:) = conn_full(:)+1;

% determine full model center of mass as well as the full model inertial axes
[full_centroid,full_surface_area,Volume,full_eigenvalues,full_eigenvectors,full_I1,full_I2,full_I_CoM,full_I_origin,full_patches] = mass_properties(pts_full,conn_full);

%% calculate mass properties of iv shaft model
[pts_plateau conn_plateau] = read_vrml_fast(fullfile(pathname,plateauivfile));
conn_plateau(:,4) = [];
conn_plateau(:) = conn_plateau(:)+1;

% determine shaft model center of mass as well as the shaft model inertial axes
[plateau_centroid,plateau_surface_area,Volume,plateau_eigenvalues,plateau_eigenvectors,plateau_I1,plateau_I2,plateau_I_CoM,plateau_I_origin,plateau_patches] = mass_properties(pts_plateau,conn_plateau);

%% determine center of ACS from the plateau crop center of mass
T=plateau_centroid;

%% assign shaft vector as the largest inertial axis of the cropped plateau
shaft_vector=plateau_eigenvectors(:,3);

% check to make sure long axis is pointing proximal by comparing its
% direction to the direction from the center of mass of the full tibia to
% the center of mass of the cropped tibial plateau
if AngleDiff(unit(plateau_centroid-full_centroid),shaft_vector) > 90
    % if the two vectors are pointing in opposite directions negate the
    % shaft vector to make it point upwards
    shaft_vector=-shaft_vector;
end

%% make sure anterior posterior axis is pointing forward
% check to make sure anterior posterior axis is poinging anterior by
% comparing its direction to the direction from the center of mass of the
% tibial plateau crop to the specified anterior point
anterior_direction=plateau_eigenvectors(:,2);
if AngleDiff(unit(anteriorpt-plateau_centroid),anterior_direction) > 90
    % if the two vectors are pointing in opposite directions negate the 2nd
    % inertial axis and make it point anterior
    anterior_direction=-anterior_direction;
end

%% creat rotation matrix from the shaft vector and remaining inertial axes
R=eye(3);
R(:,3)=shaft_vector;
R(:,2)=anterior_direction;
R(:,1)=unit(cross(anterior_direction,shaft_vector));

%% compile ACS as an RT from the R and T
ACS=[R,T';[0 0 0 1]];
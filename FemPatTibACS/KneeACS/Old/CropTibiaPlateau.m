function CropTibiaPlateau(pathname,fullivfile)
% CROPTIBIAPLATEAU
%   PATHNAME: Path where files should be loaded and saved
%   FULLIVFILE: IV tibia model filename

%% hard coded numbers
slice_thickness=0.625; % slice thickness value

%% determine slice properties of femur model, this includes the area of each slice (bounding box) and the centroid of each slice (bounding box)
[slice_properties]=SliceProperties(fullfile(pathname,fullivfile),slice_thickness);

%% load full tibia iv file
[pts_CT conn_CT] = read_vrml_fast(fullfile(pathname,fullivfile));
conn_CT(:,4) = [];
conn_CT(:) = conn_CT(:)+1;

% determine centroid and inertial axes of full femur model
[full_centroid,full_surface_area,full_volume,full_eigenvalues,full_eigenvectors,full_I1,full_I2,full_I_CoM,full_I_origin,full_patches] = mass_properties(pts_CT,conn_CT);

%% create transformation matrix from the inertial axes and center of mass
%  and then orient it in the positive z direction in order to crop upwards
RT_inertia = [full_eigenvectors; full_centroid];
pts_inertia = transformShell(pts_CT,RT_inertia,-1,1);

[max_area widest_slice_index] = max(slice_properties.area);
widest_pt = slice_properties.centroid(widest_slice_index,:); % this is the point at the widest cross sectional area

% lets figure out which way along X of the inertial axes points me towards
% the tibial platau. The centroid (0,0,0) now should be closer to the
% plautau since its larger and has more mass.
if (max(pts_inertia(:,1)) > abs(min(pts_inertia(:,1))))
    % if the max value is greater, then we are pointed the wrong way, flip
    % X & Y to keep us straight
    RT_inertia(1:3,1) = -RT_inertia(1:3,1);
    RT_inertia(1:3,2) = -RT_inertia(1:3,2);
end

% we now want to change the coordinate system, so that z points in the
% positive z direction. To do so, make z the new x, and x the negated z, we
% are basically rotating around the y axis by 90°
RT_positive_z=eye(3);
RT_positive_z(:,3)=RT_inertia(1:3,1);
RT_positive_z(:,1)=-RT_inertia(1:3,3);
RT_positive_z(:,2)=RT_inertia(1:3,2);

%% crop tibial plateau
% create a 4x4 transformation matrix from rotation matrix and bottom crop
% pt to be used for cropping
RT_crop = RT_to_fX4(RT_positive_z(1:3,1:3),widest_pt);
RT_crop_inverted = RT_crop^-1;
% crop plane at RT_crop_inverted
[pts_plateau conn_plateau]=cropIVFileToPlane(pts_CT,conn_CT,RT_crop_inverted);

%% determine centroid and inertial axes of cropped tibial plateau
[plateau_centroid,plateau_surface_area,plateau_volume,plateau_eigenvalues,plateau_eigenvectors,plateau_I1,plateau_I2,plateau_I_CoM,plateau_I_origin,plateau_patches] = mass_properties(pts_plateau,conn_plateau);

%% create transformation matrix from the inertial axes and center of mass
%  and then orient it in the positive z direction in order to make second
%  crop upwards
RT_plateau=[plateau_eigenvectors; plateau_centroid];
pts_plateau_inertia=transformShell(pts_plateau,RT_plateau,-1,1);

% lets figure out which way along Z of the inertial axes points me towards
% the tibial platau. the full centroid should be below the tibial plateau
% centroid
correct_direction=unit(plateau_centroid-full_centroid);
if AngleDiff(correct_direction,RT_plateau(1:3,3))>90
    % if the the z-axis is pointing distal, then invert it and the y-axis
    % (to keep right handed coordinate system)
    RT_plateau(1:3,2:3)=-RT_plateau(1:3,2:3);
end

%% crop tibial plateau again using the inertial axes of the tibial plateau
% create a 4x4 transformation matrix from rotation matrix and bottom crop
% pt to be used for cropping
RT_plateau_crop = RT_to_fX4(RT_plateau(1:3,1:3),widest_pt);
RT_plateau_crop_inverted = RT_plateau_crop^-1;
% crop plane at RT_plateau_crop_inverted
cropIVFileToPlane(pts_CT,conn_CT,RT_plateau_crop_inverted,fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_plateau_crop.iv']));

%% create file to visualize cropped plateau, crop points, and crop planes
fid = fopen(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_plateau_crop_visualize.iv']),'w');
fprintf(fid,createInventorHeader());
fprintf(fid,createInventorLink(fullfile(pathname,fullivfile),eye(3),[0 0 0],[1 1 1],0.5)); % full tibia model
fprintf(fid,createInventorLink(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_plateau_crop.iv']),eye(3),[0 0 0],[1 0 0],.4)); % cropped tibia model
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT_inertia(1:3,1:3),RT_inertia(4,:)));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT_plateau(1:3,1:3),RT_plateau(4,:)));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT_positive_z(1:3,1:3),widest_pt));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT_plateau(1:3,1:3),widest_pt));
fprintf(fid,createInventorSphere(full_centroid,2));
fprintf(fid,createInventorSphere(widest_pt,2));
fprintf(fid,createInventorSphere(plateau_centroid,2));
fclose('all');
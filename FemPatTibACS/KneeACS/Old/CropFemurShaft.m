function CropFemurShaft(pathname,fullivfile)
% CROPFEMURSHAFT(PATHNAME,FULLIVFILE)
%   PATHNAME: Path where files should be loaded and saved
%   FULLIVFILE: IV femur model filename

%% hard coded numbers
slice_thickness=0.625; % slice thickness value

%% determine properties of femur cross section slices
[slice_properties]=SliceProperties(fullfile(pathname,fullivfile),slice_thickness);

%% calculate mass properties of iv femur model
[pts_CT conn_CT] = read_vrml_fast(fullfile(pathname,fullivfile)); % load femur iv file
conn_CT(:,4) = [];
conn_CT(:) = conn_CT(:)+1;
% determine centroid and inertial axes of full femur model
[full_centroid,full_surface_area,full_volume,full_eigenvalues,full_eigenvectors,full_I1,full_I2,full_I_CoM,full_I_origin,full_patches] = mass_properties(pts_CT,conn_CT);

%% assign tranformation matrix based on inirtial axes and centroid
RT_full_inertia = [full_eigenvectors; full_centroid]; % full inertial based transformation matrix
pts_inertia = transformShell(pts_CT,RT_full_inertia,-1,1); % model moved to coordinate system based on inertial axes of full femur model

%% determine point where shaft begins
[area_max area_max_index]=max(slice_properties.area); % largest cross sectional area
r=range(slice_properties.area)/2;  % half range of the min and max cross sectional area
for i=1:length(slice_properties.area),d(i)=dist(slice_properties.area(i),r);end
[condyle_end condyle_end_index]=min(d(area_max_index:length(slice_properties.area))); condyle_end_index=condyle_end_index+area_max_index; % index where condyles end
shaft_start_index=round(1.3*condyle_end_index); % index where shaft begins

[min_distance min_distance_index]=min(abs(slice_properties.index-slice_properties.index(shaft_start_index)));
bottom_crop_pt=slice_properties.centroid(min_distance_index,:); % point at which to crop condyles
bottom_crop_pt(:,1:2)=full_centroid(:,1:2); % x and y pts can be same as the center of mass' x and y pts

%% crop the bottom of the femur shaft based on full inertial axes
% lets figure out which way along X of the inertial axes points me towards
% the femur condyles. The centroid (0,0,0) now should be closer to the
% condyles since its larger and has more mass.
if (max(pts_inertia(:,1)) < abs(min(pts_inertia(:,1))))
    % if the max value is smaller, then we are pointed the wrong way, flip
    % x & y to keep us straight
    RT_full_inertia(1:3,1) = -RT_full_inertia(1:3,1);
    RT_full_inertia(1:3,2) = -RT_full_inertia(1:3,2);
end;

% we now want to change the coordinate system, so that z points in the
% positive z direction. To do so, make z the new x, and x the negated z, we
% are basically rotating around the y axis by 90°
RT_positive_z=RT_full_inertia;
X_vec = RT_positive_z(1:3,1);
RT_positive_z(1:3,1) = -RT_positive_z(1:3,3);
RT_positive_z(1:3,3) = X_vec;

% create a 4x4 transformation matrix from rotation matrix and bottom crop
% pt to be used for cropping
RT4_bottom = RT_to_fX4(RT_positive_z(1:3,1:3),bottom_crop_pt);
RT4_bottom_inverted = RT4_bottom^-1;

% crop bottom of shaft
cropIVFileToPlane(pts_CT,conn_CT,RT4_bottom_inverted,fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_shaft_crop.iv']));

%% create file to visualize cropped shaft, crop points, and crop planes
fid = fopen(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_shaft_crop_visualize.iv']),'w');
fprintf(fid,createInventorHeader());

fprintf(fid,createInventorLink(fullfile(pathname,fullivfile),eye(3),[0 0 0],[1 1 1],0.5)); % full femur model
fprintf(fid,createInventorLink(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_shaft_crop.iv']),eye(3),[0 0 0],[1 0 0],.4)); % cropped shaft model

fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT_full_inertia(1:3,1:3),RT_full_inertia(4,:))); % inertial axes of full femur model
fprintf(fid,createInventorSphere(full_centroid,2)); % sphere at center of mass

fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT4_bottom(1:3,1:3),RT4_bottom(1:3,4))); % 1st bottom crop point and crop plane
fprintf(fid,createInventorSphere(bottom_crop_pt,2)); % sphere at bottom crop point

fclose('all');







%% OLD CODE THAT WAS USED TO CROP TOP OF SHAFT --> REFERENCE FOR FUTURE
% %% 2nd crop of femur shaft at the top based on inertial axes of first shaft crop
% % determine centroid and inertial axes of shaft model
% [bottom_cropped_centroid,bottom_cropped_surface_area,bottom_cropped_volume,bottom_cropped_eigenvalues,bottom_cropped_eigenvectors,bottom_cropped_I1,bottom_cropped_I2,bottom_cropped_I_CoM,bottom_cropped_I_origin,bottom_cropped_patches] = mass_properties(pts_bottom_cropped,conn_bottom_cropped);
% RT_bottom_cropped_inertia = [bottom_cropped_eigenvectors; bottom_cropped_centroid]; % shaft based transformation matrix
% pts_bottom_cropped_inertia = transformShell(pts_bottom_cropped,RT_bottom_cropped_inertia,-1,1); % model moved to coordinate system based on inertial axes of femur shaft
% 
% %lets figure out which way along X of the inertial axes points me towards
% %the femur condyles. The centroid (0,0,0) now should be closer to the
% %condyles since its larger and has more mass.
% if (max(pts_inertia(:,1)) < abs(min(pts_inertia(:,1))))
%     % if the max value is smaller, then we are pointed the wrong way, flip
%     % X & Y to keep us straight
%     RT_bottom_cropped_inertia(1:3,1) = -RT_bottom_cropped_inertia(1:3,1);
%     RT_bottom_cropped_inertia(1:3,2) = -RT_bottom_cropped_inertia(1:3,2);
% end;
% 
% %we now want to change the coordinate system, so that z points in the
% %positive z direction. To do so, make z the new x, and x the negated z, we
% %are basically rotating around the y axis by 90°
% RT_bottom_cropped_positive_z=RT_bottom_cropped_inertia;
% X_bottom_cropped_vec = RT_bottom_cropped_positive_z(1:3,1);
% RT_bottom_cropped_positive_z(1:3,1) = -RT_bottom_cropped_positive_z(1:3,3);
% RT_bottom_cropped_positive_z(1:3,3) = X_bottom_cropped_vec;
% 
% % find point at proximal end of femur model to crop in order to make crop
% % at both proximal and distal ends flush
% top_pts_index=pts_bottom_cropped(:,3)>(max(pts_bottom_cropped(:,3))-0.1);
% top_pts=pts_bottom_cropped_inertia(top_pts_index,:);
% [junk shaft_end_index]=min(pts_bottom_cropped_inertia(top_pts_index,1));
% top_crop_pt=transformShell(top_pts(shaft_end_index,:),RT_bottom_cropped_inertia,1,1);
% 
% % create a 4x4 transformation matrix from rotation matrix and top crop pt
% RT4_top=RT_to_fX4(RT_bottom_cropped_positive_z(1:3,1:3),top_crop_pt);
% RT4_top(1:3,2:3)=-RT4_top(1:3,2:3);
% RT4_top_inverted=RT4_top^-1;
% 
% % crop top of shaft
% [pts_bottom_top_cropped, conn_bottom_top_cropped]=cropIVFileToPlane(pts_bottom_cropped,conn_bottom_cropped,RT4_top_inverted);
% 
% %% 3rd crop of femur shaft at the bottom and top based on inertial axes of first two crops
% 
% % determine centroid and inertial axes of shaft model
% [bottom_top_cropped_centroid,bottom_top_cropped_surface_area,bottom_top_cropped_volume,bottom_top_cropped_eigenvalues,bottom_top_cropped_eigenvectors,bottom_top_cropped_I1,bottom_top_cropped_I2,bottom_top_cropped_I_CoM,bottom_top_cropped_I_origin,bottom_top_cropped_patches] = mass_properties(pts_bottom_top_cropped,conn_bottom_top_cropped);
% RT_bottom_top_cropped_inertia = [bottom_top_cropped_eigenvectors; bottom_top_cropped_centroid]; % shaft based transformation matrix
% pts_bottom_top_cropped_inertia = transformShell(pts_bottom_top_cropped,RT_bottom_top_cropped_inertia,-1,1); % model moved to coordinate system based on inertial axes of femur shaft
% 
% %lets figure out which way along X of the inertial axes points me towards
% %the femur condyles. The centroid (0,0,0) now should be closer to the
% %condyles since its larger and has more mass.
% if (max(pts_inertia(:,1)) < abs(min(pts_inertia(:,1))))
%     % if the max value is smaller, then we are pointed the wrong way, flip
%     % X & Y to keep us straight
%     RT_bottom_top_cropped_inertia(1:3,1) = -RT_bottom_top_cropped_inertia(1:3,1);
%     RT_bottom_top_cropped_inertia(1:3,2) = -RT_bottom_top_cropped_inertia(1:3,2);
% end;
% 
% %we now want to change the coordinate system, so that z points in the
% %positive z direction. To do so, make z the new x, and x the negated z, we
% %are basically rotating around the y axis by 90°
% RT_bottom_top_cropped_positive_z=RT_bottom_top_cropped_inertia;
% X_bottom_top_cropped_vec = RT_bottom_top_cropped_positive_z(1:3,1);
% RT_bottom_top_cropped_positive_z(1:3,1) = -RT_bottom_top_cropped_positive_z(1:3,3);
% RT_bottom_top_cropped_positive_z(1:3,3) = X_bottom_top_cropped_vec;
% % create 4x4 transformation matrix from shaft rotation matrix and bottom
% 
% % crop pt
% RT4_bottom_final = RT_to_fX4(RT_bottom_top_cropped_positive_z(1:3,1:3),bottom_crop_pt);
% RT4_bottom_inverted_final = RT4_bottom_final^-1;
% 
% % create 4x4 transformation matrix from shaft rotation matrix and top crop
% % pt
% RT4_top_final = RT_to_fX4(RT_bottom_top_cropped_positive_z(1:3,1:3),top_crop_pt);
% RT4_top_final(1:3,2:3)=-RT4_top_final(1:3,2:3);
% RT4_top_inverted_final=RT4_top_final^-1;
% 
% % crop shaft last time
% [pts_cropped_final, conn_cropped_final] = cropIVFileToPlane(pts_CT,conn_CT,RT4_bottom_inverted_final);
% cropIVFileToPlane(pts_cropped_final,conn_cropped_final,RT4_top_inverted_final,fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_shaft_crop.iv']));
% 
% %% Create file to visualize cropped shaft, crop points, and crop planes
% fid = fopen(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_shaft_crop_visualize.iv']),'w');
% fprintf(fid,createInventorHeader());
% 
% fprintf(fid,createInventorLink(fullfile(pathname,fullivfile),eye(3),[0 0 0],[1 1 1],0.5)); % full femur model
% fprintf(fid,createInventorLink(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_shaft_crop.iv']),eye(3),[0 0 0],[1 0 0],.4)); % cropped shaft model
% 
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT_full_inertia(1:3,1:3),RT_full_inertia(4,:))); % inertial axes of full femur model
% fprintf(fid,createInventorSphere(full_centroid,2)); % sphere at center of mass
% 
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT4_bottom(1:3,1:3),RT4_bottom(1:3,4))); % 1st bottom crop point and crop plane
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT4_bottom_final(1:3,1:3),RT4_bottom_final(1:3,4))); % 2nd bottom crop point and crop plane
% fprintf(fid,createInventorSphere(bottom_crop_pt,2)); % sphere at bottom crop point
% 
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT4_top(1:3,1:3),RT4_top(1:3,4))); % 1st top crop point and crop plane
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT4_top_final(1:3,1:3),RT4_top_final(1:3,4))); % 2nd bottom crop point and crop plane
% fprintf(fid,createInventorSphere(top_crop_pt,2)); % sphere at top crop point
% 
% fclose(fid);
% clear all;close all;fclose('all');
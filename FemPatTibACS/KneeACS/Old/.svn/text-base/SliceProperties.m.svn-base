function [output]=SliceProperties(ivfile,slice_thickness)

[pts_CT conn_CT] = read_vrml_fast(ivfile); % load iv file
conn_CT(:,4) = [];
conn_CT(:) = conn_CT(:)+1;

[CT_centroid,CT_sa,CT_v,CT_ev123,CT_eigenvectors,CT_I1,CT_I2,CT_I_CoM,CT_I_origin,CT_patches] = mass_properties(pts_CT,conn_CT);

RT_inertia = [CT_eigenvectors; CT_centroid]; % full inertial based transformation matrix
pts_inertia = transformShell(pts_CT,RT_inertia,-1,1); % model moved to coordinate system based on inertial axes of full femur model

% lets figure out which way along X of the inertial axes points me towards
% the femur condyles. The centroid (0,0,0) now should be closer to the
% condyles since its larger and has more mass.
if (max(pts_inertia(:,1)) < abs(min(pts_inertia(:,1))))
    
    % if the max value is smaller, then we are pointed the wrong way, flip
    % X & Y to keep us straight
    RT_inertia(1:3,1) = -RT_inertia(1:3,1);
    RT_inertia(1:3,2) = -RT_inertia(1:3,2);
    pts_inertia=transformShell(pts_CT,RT_inertia,-1,1);
    
end

% max x-coordinate
[max_x max_x_index] = max(pts_inertia(:,1));
% min x-coordinate
[min_x min_x_index] = min(pts_inertia(:,1));

for i = 1:ceil(abs(min_x-max_x)/slice_thickness)
    
%     close all
    
    % Find slice points
    poly_pts_index = find((pts_inertia(:,1) >= (min_x + (i-1)*slice_thickness) & pts_inertia(:,1) < (min_x + i*slice_thickness)));
    
    % Find length and width of bounding box
    r_y(i,1) = range(pts_inertia(poly_pts_index,2));
    r_z(i,1) = range(pts_inertia(poly_pts_index,3));
    
    area(i,1) = r_y(i,1)*r_z(i,1); % calculate area of bounding box
    centroid(i,1)=mean(pts_inertia(poly_pts_index,1)); % calculate x coordinate centroid
    centroid(i,2)=min(pts_inertia(poly_pts_index,2)) + range(pts_inertia(poly_pts_index,2))/2; % calculate y coordinate centroid
    centroid(i,3)=min(pts_inertia(poly_pts_index,3)) + range(pts_inertia(poly_pts_index,3))/2; % calculate z coordinate centroid
    
    [min_y_pt(i,2) min_index]=min(pts_inertia(poly_pts_index,2));
    min_y_pt(i,1)=centroid(i,1);
    min_y_pt(i,3)=centroid(i,3);
    
    [max_y_pt(i,2) max_index]=max(pts_inertia(poly_pts_index,2));
    max_y_pt(i,1)=centroid(i,1);
    max_y_pt(i,3)=centroid(i,3);
    
%     scatter(pts_inertia(poly_pts_index,2),pts_inertia(poly_pts_index,3));
%     hold on
%     scatter(centroid(i,2),centroid(i,3),'r');
%     close all 
end

min_y_pt_TF=transformShell(min_y_pt,RT_inertia,1,1);
max_y_pt_TF=transformShell(max_y_pt,RT_inertia,1,1);
[max_r_y max_r_y_index]=max(r_y);

% move centroid pts back to CT space
output.area=area;
output.centroid=transformShell(centroid,RT_inertia,1,1);
output.index=1:length(area);
output.min_ML_pt=min_y_pt_TF(max_r_y_index,:);
output.max_ML_pt=max_y_pt_TF(max_r_y_index,:);
output.ML_vector=unit(output.max_ML_pt-output.min_ML_pt);
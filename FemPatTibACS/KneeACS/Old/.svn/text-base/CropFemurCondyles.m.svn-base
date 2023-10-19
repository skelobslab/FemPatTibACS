function CropFemurCondyles(pathname,fullivfile,shaftivfile)
% CROPFEMURCONDYLES(PATHNAME,FULLIVFILE,SHAFTIVFILE)
%   PATHNAME: Path where files should be loaded and saved
%   FULLIVFILE: IV femur model filename
%   SHAFTIVFILE: Cropped IV femur shaft model filename

% Danny Miranda

%% hard coded numbers
maxDeviation=25; % max angle of deviation
slice_thickness=0.625; % slice thickness value
pt_multiplication_factor=300;
sample=5;

%% determine slice properties of femur model, this includes the area of each slice (bounding box) and the centroid of each slice (bounding box)
[slice_properties]=SliceProperties(fullfile(pathname,fullivfile),slice_thickness);

%% calculate mass properties of iv femur model
[pts_CT conn_CT] = read_vrml_fast(fullfile(pathname,fullivfile));
conn_CT(:,4) = [];
conn_CT(:) = conn_CT(:)+1;

% determine full model center of mass as well as the full model inertial axes
[full_centroid,full_surface_area,Volume,full_eigenvalues,full_eigenvectors,full_I1,full_I2,full_I_CoM,full_I_origin,full_patches] = mass_properties(pts_CT,conn_CT);

%% assign tranformation matrix based on inirtial axes and centroid
RT_full_inertia = [full_eigenvectors; full_centroid];

%% determine point where condyles begin

[area_max area_max_index]=max(slice_properties.area); % determine largest cross sectional area, slice where condyles is largest
r=range(slice_properties.area)/2;  % half range of the min and max cross sectional area, to determine where the condyles end

% determine the distance from each cross sectional area and the approximation of where the condyles end
for i=1:length(slice_properties.area)
    d(i)=dist(slice_properties.area(i),r);
end

% the cross sectional area that is closest to where the condyles end (r) isused to approximate where the condyles end
[condyle_end condyle_end_index]=min(d(area_max_index:length(slice_properties.area))); condyle_end_index=condyle_end_index+area_max_index;

%% determine long axis of femur shaft from its inertial axes
% determine centroid and inertial axes of femur shaft model
[shaft_centroid,shaft_surface_area,shaft_volume,shaft_eigenvalues,shaft_eigenvectors,shaft_I1,shaft_I2,shaft_I_CoM,shaft_I_origin,shaft_patches] = mass_properties(pathname,shaftivfile);
shaft_vector=shaft_eigenvectors(:,1)'; % set shaft vector as the smallest inertial axis because we want it to be the 'long' axis

% make sure shaft vector is pointing toward the distal femur
correct_direction=unit(full_centroid-shaft_centroid); % determine vector between centroid of shaft and centroid of entire femur bone to compare to shaft vector
angular_difference_shaft=AngleDiff(correct_direction,shaft_vector); % comparison of shaft vector and correct direction vector

% if the shaft vector and direction are facing in opposite directions, then invert the shaft vector
if angular_difference_shaft > 90
    shaft_vector= -shaft_vector;
end

%% determine where vector through shaft intersects bottom of condyles
distal_pt = calculateIntersectionLineIVFile(shaft_centroid,shaft_centroid + shaft_vector*300,fullfile(pathname,fullivfile));

%% create rotation matrix using the y- and z-axes from the full inertial coordinate system, and the x-axis as the shaft vector

% set x-axis as the shaft vector pointing distally. z-axis is determined 
% first because it is required to point posterior and the y-axis can 
% point in any direction
R_inertia_with_shaft=eye(3);
R_inertia_with_shaft(:,1)=shaft_vector;
R_inertia_with_shaft(:,3)=unit(cross(shaft_vector,RT_full_inertia(1:3,2)));
R_inertia_with_shaft(:,2)=unit(cross(R_inertia_with_shaft(:,3),shaft_vector));

% make sure that the z-axis determined above is pointing posterior by
% comparing its direction to the vector going from the shaft centroid to
% the full model centroid
z_is_posterior=0;
angular_difference_ap=AngleDiff(correct_direction,R_inertia_with_shaft(:,3));
if angular_difference_ap < 90
    z_is_posterior=1;
end

% negate y- and z-axes if the original z-axis is pointing anterior
if ~z_is_posterior
    R_inertia_with_shaft(:,2:3)= -R_inertia_with_shaft(:,2:3);
end


%% determine most proximal point where condyles should be cropped (on surface) the point on the surface in the direction of the z-axis
proximal_pt = calculateIntersectionLineIVFile(slice_properties.centroid(condyle_end_index,:),slice_properties.centroid(condyle_end_index,:)+R_inertia_with_shaft(:,3)'*pt_multiplication_factor,pts_CT,conn_CT);

% create a crop rotation matrix with the x-axis being the vector connecting
% the proximal point to the distal point. again, the z-axis is determined
% first because it is required to point posterior and the y-axis can point
% in any direction.
R_crop_inertia_=eye(3);
R_crop_inertia(:,1)=unit(distal_pt-proximal_pt);
R_crop_inertia(:,3)=unit(cross(R_crop_inertia(:,1),R_inertia_with_shaft(:,2)));
R_crop_inertia(:,2)=unit(cross(R_crop_inertia(:,3),R_crop_inertia(:,1)));
RT_crop_inertia_proximal_pt = [R_crop_inertia; proximal_pt];

% transform ct pts into proximal point coordinate system to perform first condyle crop
pts_proximal_inertia = transformShell(pts_CT,RT_crop_inertia_proximal_pt,-1,1);

%% 1st condyles crop based on inertial RT at proximal point
condyles_inertia_indices = pts_proximal_inertia(:,3)>0;  % indices of all condyle points
pts_condyles_inertia = pts_CT(condyles_inertia_indices,:); % find all condyle points
translation_inertia = ones(length(pts_CT),1);
translation_inertia(:) = -1;

translation_inertia(condyles_inertia_indices,:) = 1:length(pts_condyles_inertia);

conn_num_condyles_inertia = size(conn_CT,1);
conn_condyles_inertia = zeros(size(conn_CT,1),4);

numNewConn_inertia=0;
for i=1:conn_num_condyles_inertia,
    if (translation_inertia(conn_CT(i,1))~=-1 && translation_inertia(conn_CT(i,2))~=-1 && translation_inertia(conn_CT(i,3))~=-1),
        numNewConn_inertia = numNewConn_inertia+1;
        conn_condyles_inertia(numNewConn_inertia,:) = [translation_inertia(conn_CT(i,1)) translation_inertia(conn_CT(i,2)) translation_inertia(conn_CT(i,3)) -1];
    end    
end;
conn_condyles_inertia(conn_condyles_inertia(:,1)==0,:)=[];

%% cylinder fit to 1st condyles crop
% determine condyle bounding box dims
dim(1)=max(pts_condyles_inertia(:,1))-min(pts_condyles_inertia(:,1));
dim(2)=max(pts_condyles_inertia(:,2))-min(pts_condyles_inertia(:,2));
dim(3)=max(pts_condyles_inertia(:,3))-min(pts_condyles_inertia(:,3));
% determine index's for dims
dim_max=find(dim==max(dim));
dim_mid=find(dim==median(dim));
dim_min=find(dim==min(dim));
% determine index of pts_condyles that refer to the most medial pt and most lateral pt
a0pt1_=find(pts_condyles_inertia(:,dim_max)==max(pts_condyles_inertia(:,dim_max)));
a0pt2_=find(pts_condyles_inertia(:,dim_max)==min(pts_condyles_inertia(:,dim_max)));

x0 = mean(pts_condyles_inertia)'; % estimate pt on axis
a0 = unit(pts_condyles_inertia(a0pt2_,:)-pts_condyles_inertia(a0pt1_,:))'; % estimate axis direction
r0 = (((max(pts_condyles_inertia(:,dim_mid))-min(pts_condyles_inertia(:,dim_mid)))/2)+((max(pts_condyles_inertia(:,dim_min))-min(pts_condyles_inertia(:,dim_min)))/2))/2; % estimate radius
tolp = 0.1;
tolg = 0.1;

% cylinder fit
[x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog,a, R0, R] = lscylinder(pts_condyles_inertia(1:sample:end,:), x0, a0, r0, tolp, tolg);

%% repeat finding the crop planes based on the axis through the cylinder fit of the original condyle cropping

% make sure both inertia medial lateral axis and cylinder medial lateral axis are pointing in the same direction
if AngleDiff(unit(an),unit(R_crop_inertia(:,2)))>90
    an=-an; % if directions are not the same negate cylinder axis
end

% create a new crop rotation matrix with the x-axis being the vector
% connecting the proximal point to the distal point. again, the z-axis is 
% determined first because it is required to point posterior and the y-axis
% can point in any direction. instead of using the inertial axis to 
% determine the z-axis, the vector through the cylinder fit is used.
R_crop_cylinder=eye(3);
R_crop_cylinder(:,1)=unit(distal_pt-proximal_pt);
R_crop_cylinder(:,3)=unit(cross(R_crop_cylinder(:,1),an));
R_crop_cylinder(:,2)=unit(cross(R_crop_cylinder(:,3),R_crop_cylinder(:,1)));
RT_crop_cylinder_proximal_pt = [R_crop_cylinder; proximal_pt];

% transform ct pts into proximal point coordinate system to perform second condyle crop
pts_proximal_cylinder = transformShell(pts_CT,RT_crop_cylinder_proximal_pt,-1,1);

%% 2nd condyles crop based on cylinder fit or original condyles crop
condyles_cylinder_indices = pts_proximal_cylinder(:,3)>0;  % indices of all condyle points
pts_condyles_cylinder = pts_CT(condyles_cylinder_indices,:); % find all condyle points
translation_cylinder = ones(length(pts_CT),1);
translation_cylinder(:) = -1;

translation_cylinder(condyles_cylinder_indices,:) = 1:length(pts_condyles_cylinder);

conn_num_condyles_cylinder = size(conn_CT,1);
conn_condyles_cylinder = zeros(size(conn_CT,1),4);

numNewConn_cylinder=0;
for i=1:conn_num_condyles_cylinder,
    if (translation_cylinder(conn_CT(i,1))~=-1 && translation_cylinder(conn_CT(i,2))~=-1 && translation_cylinder(conn_CT(i,3))~=-1),
        numNewConn_cylinder = numNewConn_cylinder+1;
        conn_condyles_cylinder(numNewConn_cylinder,:) = [translation_cylinder(conn_CT(i,1)) translation_cylinder(conn_CT(i,2)) translation_cylinder(conn_CT(i,3)) -1];
    end    
end;
conn_condyles_cylinder(conn_condyles_cylinder(:,1)==0,:)=[];

% write iv file with only condyles
patch2iv(pts_condyles_cylinder,conn_condyles_cylinder(:,1:3),fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_condyles_crop.iv']));

%% Create axis to compare triangle normals
% % cylinder fit parameters
% axis_pt_est = mean(pts_condyle)';
% axis_dir_est = [0 0 1]'; % orient in Z
% radius_est = (max(pts_condyle(:,1))-min(pts_condyle(:,1)))/2;
% tolp = 0.1;
% tolg = 0.1;
% % cylinder fit
% [x0n, condyle_axis, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog, a, R0, R] = lscylinder(pts_condyle(1:10:end,:), axis_pt_est, axis_dir_est, radius_est, tolp, tolg);
% 
% %% Find all the surfaces within 25 degrees of the cylinder normals
% 
% % lets create the normal for each patch, and then check if we want to keep
% % it or not
% conn_reduced_condyles = conn_condyles;
% condyles_reduced_indices = ones(length(conn_reduced_condyles),1);
% for i=1:length(conn_reduced_condyles),
%     v1_2 = unit(pts_condyle(conn_reduced_condyles(i,1),:) - pts_condyle(conn_reduced_condyles(i,2),:));
%     v3_2 = unit(pts_condyle(conn_reduced_condyles(i,3),:) - pts_condyle(conn_reduced_condyles(i,2),:));
%     normal_axis = cross(v1_2, v3_2);
%     angular_difference_shaft = acosd(dot(normal_axis,condyle_axis));
%     if (abs(angular_difference_shaft-90)>maxDeviation)
%         condyles_reduced_indices(i)=0;
%     end;
% end;
% 
% %% Discard all triangles that are alone in 3D space (unconnected)
% [pts_reduced_condyles conn_reduced_condyles] = selectSubregionSurfaceByPatches(pts_condyle,conn_reduced_condyles, condyles_reduced_indices==1);
% [condyle_pieces] = splitSurface(pts_reduced_condyles,conn_reduced_condyles);
% num_condyle_pieces = zeros(length(condyle_pieces),1);
% for i=1:length(condyle_pieces), num_condyle_pieces(i)=size(condyle_pieces(i).conn,1); end;
% [max_piece max_piece_index] = max(num_condyle_pieces);
% patch2iv(condyle_pieces(max_piece_index).pts,condyle_pieces(max_piece_index).conn(:,1:3),fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_final_condyles_crop.iv']));

%% create file to visualize cropped condyles
fid = fopen(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_condyles_crop_visualize.iv']),'w');
fprintf(fid,createInventorHeader());
fprintf(fid,createInventorLink(fullfile(pathname,fullivfile),eye(3),[0 0 0],[1 1 1],0.4)); % full femur model
fprintf(fid,createInventorLink(fullfile(pathname,[fullivfile(1:length(fullivfile)-3),'_condyles_crop.iv']),eye(3),[0 0 0],[1 0 0],0)); % cropped shaft model
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',full_eigenvectors,full_centroid)); % inertial axes of full femur model
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',shaft_eigenvectors,shaft_centroid)); % inertial axes of shaft
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',R_inertia_with_shaft,distal_pt)); % full and shaft inertial based axes at distal point
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',R_crop_inertia,proximal_pt)); % crop axes based on inertial axis
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',R_crop_cylinder,proximal_pt)); % crop axes based on cylinder fit axis
fprintf(fid,createInventorSphere(distal_pt,2)); % sphere at distal crop point
fprintf(fid,createInventorSphere(proximal_pt,2)); % sphere at slice crop point
fprintf(fid,createInventorSphere(full_centroid,2)); % sphere at center of mass
fprintf(fid,createInventorSphere(shaft_centroid,2)); % sphere at center of mass of shaft
fclose('all');
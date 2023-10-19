function fACS = buildfACS(pname, fm)
% fACS = buildfACS(pname, fm)
%
%   fACS:   4x4 matrix containing medial-lateral axis, anterior-posterior
%           axis, long axis, and origin of the femoral anatomical
%           coordinate system
%
%   pname:  path where 3-D femur iv model is located
%   fm:     file name of 3-D femur iv model (e.g. femur0123.iv)
%
%   This script builds an ACS from a 3-D model of the distal femur using
%   it's diaphysis and condyles
%Miranda DL, Rainbow MJ, Leventhal EL, Crisco JJ, Fleming BC. 
%Automatic determination of anatomical coordinate systems for 
%three-dimensional bone models of the isolated human knee. 
%J Biomech. 2010 May 28;43(8):1623–6. 

%% Hard coded numbers
slice_thickness = 0.625;    % slice thickness value
pad = 2;    % pt padding for cylinder fit

%% Initial steps
%   load points and connections of 3-D femur model 
[pts.m conn.m] = read_vrml_fast(fullfile(pname, fm)); % load points and connections of 3-D femur model
conn.m(:,4) = [];
conn.m(:) = conn.m(:)+1;

%   Determine inertial properties
[centroid.m,sa.m,v.m,evals.m,inertial_axes.m,I1.m,I2.m,CoM.m,I_origin.m,patches.m] = mass_properties(pts.m,conn.m);

%   Create transformation matrix using the 3-D femur models inertial axes
%   and centroid
rt.i = [inertial_axes.m; centroid.m];

%   Register femur model to it's inertial axes and centroid
pts.i = transformShell(pts.m, rt.i, -1, 1);

%% Determine axial slice properties of the 3-D femur model
slice_properties = sliceProperties(pts.m, pts.i, rt.i, slice_thickness);

%% Isolate femoral diaphysis
% The first step in building the femoral ACS is isolating only the points
% comprising the femoral diaphysis
[pts.d conn.d] = fDiaph(pts.m, conn.m, pts.i, rt.i, slice_properties, centroid.m, pname, fm);

% determine centroid and inertial axes of femur diaphysis model
[centroid.d,sa.d,v.d,evals.d,inertial_axes.d,I1.d,I2.d,CoM.d,I_origin.d,patches.d] = mass_properties(pts.d, conn.d);
diaphysis_vector=inertial_axes.d(:,1)'; % set diaphysis vector as the smallest inertial axis because we want it to be the 'long' axis

%% Isolate femoral condyles
% The second step in building the femoral ACS is isolating only the points
% comprising the femoral condyles
[pts.c conn.c] = fCondyles(pts.m, conn.m, rt.i, slice_properties, centroid.m, centroid.d, diaphysis_vector, pname, fm);

% determine centroid and inertial axes of condyles model
[centroid.c,sa.c,v.c,evals.c,inertial_axes.c,I1.c,I2.c,CoM.c,I_origin.c,patches.c] = mass_properties(pts.c,conn.c);

%% Determine initial conditions and fit cylinder to condyles
% determine condyle bounding box dims
dim(1)=max(pts.c(:,1))-min(pts.c(:,1));
dim(2)=max(pts.c(:,2))-min(pts.c(:,2));
dim(3)=max(pts.c(:,3))-min(pts.c(:,3));
% determine index's for dims
dim_max=find(dim==max(dim));
dim_mid=find(dim==median(dim));
dim_min=find(dim==min(dim));
% determine index of pts.c that refer to the most medial pt and most lateral pt
a0pt1_ = pts.c(:,dim_max)==max(pts.c(:,dim_max));
a0pt2_ = pts.c(:,dim_max)==min(pts.c(:,dim_max));

x0 = mean(pts.c)'; % estimate pt on axis
a0 = unit(pts.c(a0pt2_,:)-pts.c(a0pt1_,:))'; % estimate axis direction
r0 = (((max(pts.c(:,dim_mid))-min(pts.c(:,dim_mid)))/2)+((max(pts.c(:,dim_min))-min(pts.c(:,dim_min)))/2))/2; % estimate radius
tolp = 0.1;
tolg = 0.1;

% cylinder fit
[x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog,a, R0, R] = lscylinder(pts.c(1:pad:end,:), x0, a0, r0, tolp, tolg);

% determine the length of the condyles, from medial to lateral
condyles_length=abs(max(pts.c(:,1))-min(pts.c(:,2)));

%% Calculate femoral ACS
fACS = femurACS(centroid.m, centroid.d, diaphysis_vector, an', [x0n-(0.5*condyles_length)*an]', condyles_length);

% save femur ACS files
dlmwrite(fullfile(pname,[fm(1:length(fm)-3),'_ACS.txt']),fACS,'delimiter','\t','precision','%.10f');
function [fACS tACS] = kneeACSpackage(pathname, femurivfile, tibiaivfile, anteriortibiapt)
% [fACS tACS] = kneeACSpackage(pathname, femurivfile, tibiaivfile, anteriortitbiapt)
%
%   INPUT
%   pathname: path where iv files are located
%   femurivfile: full femur model file name
%   tibiaivfile: full tibia model file name
%   anteriortibiapt: 3-D point located on the anterior half of tibia
%
%   OUTPUT
%   fACS: femur ACS output
%   tACS: tibia ACS output
tstart=tic;
% pt padding for cylinder fit
pad=3;

% obtain femur shaft crop from full femur model
CropFemurShaft(pathname,femurivfile);

% assign shaft iv file name
shaftivfile=[femurivfile(1:length(femurivfile)-3),'_shaft_crop.iv'];

% obtain femur condyles crop from full femur model
CropFemurCondyles(pathname,femurivfile,shaftivfile);

% assign condyles iv file name
condylesivfile=[femurivfile(1:length(femurivfile)-3),'_condyles_crop.iv'];

% read in condyles points and connections
[pts_condyles conn_condyles] = read_vrml_fast(fullfile(pathname,condylesivfile));
conn_condyles(:,4) = [];
conn_condyles(:) = conn_condyles(:)+1;
% determine centroid and inertial axes of condyles model
[condyles_centroid,condyles_surface_area,condyles_volume,condyles_eigenvalues,condyles_eigenvectors,condyles_I1,condyles_I2,condyles_I_CoM,condyles_I_origin,condyles_patches] = mass_properties(pts_condyles,conn_condyles);

% set cylinder fit initial conditions

% determine condyle bounding box dims
dim(1)=max(pts_condyles(:,1))-min(pts_condyles(:,1));
dim(2)=max(pts_condyles(:,2))-min(pts_condyles(:,2));
dim(3)=max(pts_condyles(:,3))-min(pts_condyles(:,3));
% determine index's for dims
dim_max=find(dim==max(dim));
dim_mid=find(dim==median(dim));
dim_min=find(dim==min(dim));
% determine index of pts_condyles that refer to the most medial pt and most lateral pt
a0pt1_=find(pts_condyles(:,dim_max)==max(pts_condyles(:,dim_max)));
a0pt2_=find(pts_condyles(:,dim_max)==min(pts_condyles(:,dim_max)));

x0 = mean(pts_condyles)'; % estimate pt on axis
a0 = unit(pts_condyles(a0pt2_,:)-pts_condyles(a0pt1_,:))'; % estimate axis direction
r0 = (((max(pts_condyles(:,dim_mid))-min(pts_condyles(:,dim_mid)))/2)+((max(pts_condyles(:,dim_min))-min(pts_condyles(:,dim_min)))/2))/2; % estimate radius
tolp = 0.1;
tolg = 0.1;

% cylinder fit
[x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog,a, R0, R] = lscylinder(pts_condyles(1:pad:end,:), x0, a0, r0, tolp, tolg);

% determine the length of the condyles, from medial to lateral
condyles_length=abs(max(pts_condyles(:,1))-min(pts_condyles(:,2)));

% calculate femur ACS
fACS=FemurACS(pathname,femurivfile,shaftivfile,an',[x0n-(0.5*condyles_length)*an]',condyles_length);

% obtain tibia plateau crop
CropTibiaPlateau(pathname,tibiaivfile)

% assign tibia platea iv file name
plateauivfile=[tibiaivfile(1:length(tibiaivfile)-3),'_plateau_crop.iv'];

% calculate tibia ACS
tACS=TibiaACS(pathname,tibiaivfile,plateauivfile,anteriortibiapt);

% save femur and tibia ACS files
dlmwrite(fullfile(pathname,[femurivfile(1:length(femurivfile)-3),'_ACS.txt']),fACS,'delimiter','\t','precision','%.3f');
dlmwrite(fullfile(pathname,[tibiaivfile(1:length(tibiaivfile)-3),'_ACS.txt']),tACS,'delimiter','\t','precision','%.3f');
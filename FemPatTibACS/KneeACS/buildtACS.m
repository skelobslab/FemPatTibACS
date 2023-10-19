function tACS = buildtACS(pname, tm, anteriorpt)
% tACS = buildtACS(pname, fm, anteriorpt)
%
%   tACS:       4x4 matrix containing medial-lateral axis, anterior-posterior
%               axis, long axis, and origin of the tibial anatomical
%               coordinate system
%
%   pname:      path where 3-D tibia iv model is located
%   tm:         file name of 3-D tibia iv model (e.g. tibia0123.iv)
%   anteriorpt: any point on the anterior half of the tibia
%
%   This script builds an ACS from a 3-D model of the proximal tibia using
%   it's plateau
%Miranda DL, Rainbow MJ, Leventhal EL, Crisco JJ, Fleming BC. 
%Automatic determination of anatomical coordinate systems for 
%three-dimensional bone models of the isolated human knee. 
%J Biomech. 2010 May 28;43(8):1623–6. 
%% Hard coded numbers
slice_thickness = 0.625;    % slice thickness value

%% Initial steps
%   load points and connections of 3-D tibia model
[pts.m conn.m] = read_vrml_fast(fullfile(pname, tm)); % load points and connections of 3-D tibia model
conn.m(:,4) = [];
conn.m(:) = conn.m(:)+1;

%   Determine inertial properties
[centroid.m,sa.m,v.m,evals.m,inertial_axes.m,I1.m,I2.m,CoM.m,I_origin.m,patches.m] = mass_properties(pts.m,conn.m);

%   Create transformation matrix using the 3-D tibia models inertial axes
%   and centroid
rt.i = [inertial_axes.m; centroid.m];

%   Register tibia model to it's inertial axes and centroid
pts.i = transformShell(pts.m, rt.i, -1, 1);

%% Determine axial slice properties of the 3-D tibia model
slice_properties = sliceProperties(pts.m, pts.i, rt.i, slice_thickness);

%% Isolate tibial plateau
% The first step in building the tibial ACS is isolating only the points
% comprising the plateau
[pts.p conn.p] = tPlateau(pts.m, conn.m, pts.i, rt.i, centroid.m, slice_properties, pname, tm);

% determine centroid and inertial axes of femur diaphysis model
[centroid.p,sa.p,v.p,evals.p,inertial_axes.p,I1.p,I2.p,CoM.p,I_origin.p,patches.p] = mass_properties(pts.p, conn.p);

%% Calculate tibial ACS
tACS = tibiaACS(centroid.m, centroid.p, inertial_axes.p, anteriorpt);

% save femur ACS files
dlmwrite(fullfile(pname,[tm(1:length(tm)-3),'_ACS.txt']),tACS,'delimiter','\t','precision','%.10f');
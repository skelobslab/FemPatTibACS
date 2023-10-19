function [] = DannyTestTibiaACS()
%%
grayfile = 'F:\LocalCopies\DannyKnee\226 Tibia_grayvalues.txt';
ivFile = 'F:\LocalCopies\DannyKnee\HUMAN KNEES_226 Tibia_001.iv';
saveFile = 'F:\LocalCopies\DannyKnee\cropTest.iv';
%%

centroidData = readGreyvalueFile(grayfile);

%%
centroid = cell2mat({centroidData(:).centroid}');
area = cell2mat({centroidData(:).area}');
locations = cell2mat({centroidData(:).location}');

num_slices = length(centroidData);

% TODO: need to determine the height of the transition from diaphysis to
% metaphysis
scatter(locations,area);

%%
[pts conn] = read_vrml_fast(ivFile);
conn(:,4) = [];
conn(:) = conn(:)+1;
%%
[Centroid,SurfaceArea,Volume,CoM_ev123,CoM_eigenvectors,I1,I2,I_CoM,I_origin,patches] = mass_properties(ivFile);
%%
RT = [CoM_eigenvectors; Centroid];
pts2 = transformShell(pts,RT,-1,1);

[maxarea widestSliceIndex] = max(area);
loc = centroid(widestSliceIndex,:); % this is the z-value we will want!

%lets figure out which way along X of the inertial axes points me towards
%the tibial platau. The centroid (0,0,0) now should be closer to the
%plautau since its larger and has more mass.
if (max(pts2(:,1)) > abs(min(pts2(:,1))))
    % if the max value is greater, then we are pointed the wrong way, flip
    % X & Y to keep us straight
    RT(1:3,1) = -RT(1:3,1);
    RT(1:3,2) = -RT(1:3,2);
end;

%we now want to change the coordinate system, so that z points in the
%positive z direction. To do so, make z the new x, and x the negated z, we
%are basically rotating around the y axis by 90°
X_vec = RT(1:3,1);
RT(1:3,1) = -RT(1:3,3);
RT(1:3,3) = X_vec;


%%
%lets crop the file now....
RT4 = RT_to_fX4(RT(1:3,1:3),loc);
RT4 = RT4^-1;
cropIVFileToPlane(pts, conn,RT4,saveFile);




%%
% ACS = RT;
fid = fopen('F:\LocalCopies\DannyKnee\test2.iv','w');
fprintf(fid,createInventorHeader());
fprintf(fid,createInventorLink(ivFile,eye(3),[0 0 0],[1 1 1],0.4));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',RT(1:3,1:3),RT(4,:)));
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',ACS2(1:3,1:3),ACS2(4,:)));
% fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',ACS3(1:3,1:3),ACS3(4,:)));
% makeCylIV(fid, 1, [1 0 0], ACS(4,:)', ACS(1:3,1), 50);
% makeCylIV(fid, 1, [0 1 0], ACS(4,:)', ACS(1:3,2), 20);
% makeCylIV(fid, 1, [0 0 1], ACS(4,:)', ACS(1:3,3), 20);
% makeCylIV(fid, 1, [0 0 1], shaft1, -x, 200);
% fprintf(fid,createInventorSphere(headJunction,2));
fprintf(fid,createInventorSphere(loc,2));
% fprintf(fid,createInventorSphere(distal,3));
% fprintf(fid,createInventorSphere(metCent,1));
% fprintf(fid,createInventorSphere(styCent,1));
fclose(fid);

% ! F:\LocalCopies\DannyKnee\test.iv
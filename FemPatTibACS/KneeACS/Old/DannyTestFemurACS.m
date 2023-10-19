function [] = DannyTestFemurACS(grayfile,ivFile)
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
%%
b = -30;
a = -110;

threshold1 = -175;

[minDistance index] = min(abs(locations-threshold1));
headJunction = centroid(index,:);

[minDistance index] = min(abs(locations-a));
shaft1 = centroid(index,:);
[minDistance index] = min(abs(locations-b));
shaft2 = centroid(index,:);

shaft_vector = shaft1-shaft2;

%%
intersects = calculateIntersectionLineIVFile(shaft1 + shaft_vector*300,shaft1 - shaft_vector*300,ivFile);
%%
origin = intersects(1,:);

AB = -shaft_vector; %need to flip it so it points down.
start_z = CoM_eigenvectors(:,3);
x = unit(AB);
y = unit(cross(x,start_z));
z = unit(cross(x,y));
R = [x',y',z'];
ACS2 = [R; origin];

%%
%now lets get the point on the surface.... in the direction of the z-axis
slicePoint = calculateIntersectionLineIVFile(headJunction,headJunction + ACS2(1:3,3)'*300,ivFile);



%%
newX = unit(slicePoint - origin);
newZ = unit(cross(newX,y));
R2 = [newX',y',newZ'];
ACS3 = [R2; origin];
pts2 = transformShell(pts,ACS3,-1,1);

% pts3 = pts2(pts2(:,3)>0,:);
%now lets find a way to crop this ting....

%%

indices = pts2(:,3)>0;

% new_pts = pts(find(angles>startAngle & angles<endAngle),:);
new_pts = pts(indices,:);
translation = ones(length(pts),1);
translation(:) = -1;

translation(indices,:) = 1:length(new_pts);

num_con = size(conn,1);
new_con = zeros(size(conn,1),4);
%%
numNewCon=0;
for i=1:num_con,
    if (translation(conn(i,1))~=-1 && translation(conn(i,2))~=-1 && translation(conn(i,3))~=-1),
        numNewCon = numNewCon+1;
        new_con(numNewCon,:) = [translation(conn(i,1)) translation(conn(i,2)) translation(conn(i,3)) -1];
    end    
end;
new_con(new_con(:,1)==0,:)=[];


patch2iv(new_pts,new_con(:,1:3),'F:\LocalCopies\DannyKnee\subregion.iv');
%%

%first cylinder fitting
X = new_pts(1:10:end,:);
x0 = mean(new_pts)';
a0 = [0 0 1]'; % orient in Z
r0 = (max(new_pts(:,1))-min(new_pts(:,1)))/2;
tolp = 0.1;
tolg = 0.1;
[x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog, ... 
          a, R0, R] = lscylinder(X, x0, a0, r0, tolp, tolg);

      
% Okay, now lets go and find all the surfaces within X° from the cylinder
% normals

maxDeviation = 20; 

%lets create the normal for each patch, and then check if we want to keep
%it or not
new_conn = new_con;
indices2 = ones(length(new_conn),1);
for i=1:length(new_conn),
    v1_2 = unit(new_pts(new_conn(i,1),:) - new_pts(new_conn(i,2),:));
    v3_2 = unit(new_pts(new_conn(i,3),:) - new_pts(new_conn(i,2),:));
    normal = cross(v1_2, v3_2);
    angle = acosd(dot(normal,an));
    if (abs(angle-90)>maxDeviation)
        indices2(i)=0;
    end;
end;

%%
[pts2 conn2] = selectSubregionSurfaceByPatches(new_pts,new_conn, indices2==1);


patch2iv(pts2,conn2(:,1:3),'F:\LocalCopies\DannyKnee\out2.iv');
%%

%now lets try the second cylinder fitting on the new points
% create second cylinder
X = pts2(1:end,:);
x0 = mean(pts2)';
a0 = [0 0 1]'; % orient in Z
r0 = (max(pts2(:,1))-min(pts2(:,1)))/2;
tolp = 0.1;
tolg = 0.1;
[x0n2, an2, rn2, d, sigmah, conv, Vx0n, Van, urn, GNlog, ... 
          a, R0, R] = lscylinder(X, x0, a0, r0, tolp, tolg);

%%
l = 100;
ACS = RT;
fid = fopen('F:\LocalCopies\DannyKnee\test.iv','w');
fprintf(fid,createInventorHeader());
fprintf(fid,createInventorLink(ivFile,eye(3),[0 0 0],[1 1 1],0.4));
fprintf(fid,createInventorLink('F:\LocalCopies\DannyKnee\subregion.iv',eye(3),[0 0 0]));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',ACS(1:3,1:3),ACS(4,:)));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',ACS2(1:3,1:3),ACS2(4,:)));
fprintf(fid,createInventorLink('P:\TheCollective\IV_Files\ACS.iv',ACS3(1:3,1:3),ACS3(4,:)));
fprintf(fid, createInventorCylinder(x0n,an,l,rn,[0 1 0],0.7));
fprintf(fid, createInventorCylinder(x0n2,an2,l,rn2,[1 0 0],0.7));
% makeCylIV(fid, 1, [1 0 0], ACS(4,:)', ACS(1:3,1), 50);
% makeCylIV(fid, 1, [0 1 0], ACS(4,:)', ACS(1:3,2), 20);
% makeCylIV(fid, 1, [0 0 1], ACS(4,:)', ACS(1:3,3), 20);
fprintf(fid,createInventorSphere(headJunction,2));
fprintf(fid,createInventorSphere(slicePoint,2));
% fprintf(fid,createInventorSphere(distal,3));
% fprintf(fid,createInventorSphere(metCent,1));
% fprintf(fid,createInventorSphere(styCent,1));
fclose(fid);

% ! F:\LocalCopies\DannyKnee\test.iv
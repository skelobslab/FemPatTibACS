function [patT] = patellaACS(ivPath,ivFile,side)
%function [patT] = patellaACS(ivPath,ivFile,side)
%Generate an Automated ACS of the human patella based on surface topography. 
%Based on:  Rainbow, M. J. et al. Automatic determination of an anatomical
%coordinate system for a three-dimensional model of the human patella. 
%J Biomech (2013). doi:10.1016/j.jbiomech.2013.05.024
%
%INPUTS
%ivPath = directory containing 3D patella model
%ivFile = iv or wrl file of 3D patella model
%side = 'R' (Right) or 'L' (Left)
%
%OUTPUTS
%patT = 4 x 4 matrix containing the anatomical cooridnate system where,
%X-axis (medial / lateral) = patT(1:3,1), Y-axis (anterior / posterior = patT(1:3,2),
%Z-axis (proximal / distal) = patT(1:3,3), origin = patT(1:3,4)


global objf_out;


Xr = [0 0 0];
Yr = [0 0 0];
Zr = [0 0 0];

[Centroid,SurfaceArea,Volume,CoM_ev123,CoM_eigenvectors,I1,I2,I_CoM,I_origin,patches] = mass_properties(fullfile(ivPath,ivFile));
[pts cnt] = read_vrml_fast(fullfile(ivPath,ivFile));

%set eig3 to z-axis (A/P axis)
ACS_L_P_x = CoM_eigenvectors(1:3,1);
ACS_L_P_y = CoM_eigenvectors(1:3,2);
ACS_L_P_z = CoM_eigenvectors(1:3,3);

%eigenvectors may be flipped 180 degree from intended orienation (from
%posterior to anterior)  Use patella shape to check axis orienation
[ACS_L_P_z ACS_L_P_x] = checkPatellaAP(ACS_L_P_x, ACS_L_P_y, ACS_L_P_z,Centroid,pts);

%correct z-axis (so it is oriented from posterior to anterior)
ACS_P_temp = [unit(ACS_L_P_x) unit(ACS_L_P_y) unit(ACS_L_P_z)];
ACS_steer = patellaIGuess(pts, ACS_L_P_x, ACS_L_P_y, ACS_L_P_z, Centroid);

lb = -89.5;
ub = 89.5;
gstd = [0];

options = optimset('Algorithm','interior-point','MaxFunEvals',100000,'TolX',1e-8,'Display','notify');
counter = 1;
[x,fval,exitflag]=fmincon(@evalRidgeDistances,5,[],[],[],[],lb,ub,[],options); 
Zrot = rotamatz(x,'deg');
patT = ACS_steer * Zrot;

temp = [patT(1:3,2) patT(1:3,3) patT(1:3,1)];
patT = RT_to_fX4(temp, patT(1:3,4)');

if strcmp(side,'L')
    patT(:,1) = -patT(:,1);
    patT(:,3) = -patT(:,3);
end

patT = patellaACSRefine(ivPath,ivFile,patT,side);

function objf = evalRidgeDistances(coordRot)
    samples = 200;
    zpercent = 0.25;
    zrot = rotamatz(coordRot,'deg');
    ACSrot = ACS_steer * zrot;
    steerPtsL = RT_transform(pts,ACSrot(1:3,1:3),ACSrot(1:3,4)',0);
    moveaxis_inc = linspace((0.9 * min(steerPtsL(:,1))),(0.8 * max(steerPtsL(:,1))),samples); 
    
    
    for i = 1:samples
        
         I_poi = find(steerPtsL(:,1) > moveaxis_inc(i) & steerPtsL(:,1) < (moveaxis_inc(i) + 2)   & steerPtsL(:,3) < (zpercent * min(steerPtsL(:,3))) );
        
        tempPts = steerPtsL(I_poi,:);
        

        
        [yy II] = min(tempPts(:,3));
        
        xdev(i,1) = tempPts(II,2);
    end %i
    
        
    objf = std(xdev);
%     objf = -abs(max(xdev)-min(xdev));
    gstd(counter,1) = objf;
    gstd(counter,2) = coordRot;
    counter = counter + 1;
%     objf_out(s,1) = objf;
    obj_out = objf;
    
end %objective fucntion


end %main function


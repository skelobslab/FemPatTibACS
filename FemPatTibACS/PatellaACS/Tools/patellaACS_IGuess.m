function [ACS_steer] = patellaACS_IGuess(ivPath,ivFile,side)



%check: may only be needed for side determination

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

ACS_steer = [ACS_steer(:,2) ACS_steer(:,3) ACS_steer(:,1) ACS_steer(:,4)]; 

if strcmp(side,'L')
    ACS_steer(:,1) = -ACS_steer(:,1);
    ACS_steer(:,3) = -ACS_steer(:,3);
end



end %main function


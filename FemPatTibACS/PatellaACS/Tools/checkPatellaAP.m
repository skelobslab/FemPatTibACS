function [ACS_L_P_z ACS_L_P_x] = checkPatellaAP(ACS_L_P_x, ACS_L_P_y, ACS_L_P_z,Centroid,pts)
%This function fits a 4th order polynomial to the front and back of the
%patella and choses the side with the best R^2 as the front.  Need to
%implement a better way to do this.

surface_res = 0.5;
fit_order = 4;

ACS_P_temp = [ACS_L_P_x ACS_L_P_y ACS_L_P_z];
templocalpts = RT_transform(pts,ACS_P_temp,Centroid,0);

%z < 0 fit
Itmp = find(templocalpts(:,3) < 0);
Mn = polyfitn(templocalpts(Itmp,1:2),templocalpts(Itmp,3),fit_order);
zfit = polyvaln(Mn,templocalpts(Itmp,1:2));
[Xm Ym] = meshgrid(min(templocalpts(Itmp,1)):surface_res:max(templocalpts(Itmp,1)),min(templocalpts(Itmp,2)):surface_res:max(templocalpts(Itmp,2)));
for i = 1:size(Xm,2)
    Zm(:,i) = polyvaln(Mn,[Xm(:,i) Ym(:,i)]); 
end %i
   
%z > 0 fit
Itmp = find(templocalpts(:,3) > 0);
Mp = polyfitn(templocalpts(Itmp,1:2),templocalpts(Itmp,3),fit_order);
zfit = polyvaln(Mp,templocalpts(Itmp,1:2));
[Xm Ym] = meshgrid(min(templocalpts(Itmp,1)):surface_res:max(templocalpts(Itmp,1)),min(templocalpts(Itmp,2)):surface_res:max(templocalpts(Itmp,2)));
clear Zm
for i = 1:size(Xm,2)
    Zm(:,i) = polyvaln(Mp,[Xm(:,i) Ym(:,i)]); 
end %i
    
%check direction of 3rd inertial axis    
    z_axis = [0 0 1];
    y_axis = [0 1 0];
    x_axis = [1 0 0];
    pos_fit = 1;
    if Mp.R2 < Mn.R2
        ACS_L_P_z = -ACS_L_P_z;
        ACS_L_P_x = -ACS_L_P_x;
    else
        ACS_L_P_z = ACS_L_P_z;
        ACS_L_P_x = ACS_L_P_x;        
    end

end %function
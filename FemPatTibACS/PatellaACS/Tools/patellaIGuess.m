function ACS_steer = patellaIGuess(pts, ACS_L_P_x, ACS_L_P_y, ACS_L_P_z, Centroid)
%This function uses curvature to determine the lateral aspect of the
%patella.

ACS_P_temp = [unit(ACS_L_P_x) unit(ACS_L_P_y) unit(ACS_L_P_z)];

articular_fit_order = 6;
surface_res = 0.5;

templocalpts = RT_transform(pts,ACS_P_temp,Centroid,0);
percentO = 0.1; %fit polynomial to posterior surface
Itmp = find(templocalpts(:,3) < min(templocalpts(:,3)) .* percentO);
M = polyfitn(templocalpts(Itmp,1:2),templocalpts(Itmp,3),articular_fit_order);
zfit = polyvaln(M,templocalpts(Itmp,1:2));
[Xm Ym] = meshgrid(min(templocalpts(Itmp,1)):surface_res:max(templocalpts(Itmp,1)),min(templocalpts(Itmp,2)):surface_res:max(templocalpts(Itmp,2)));
clear Zm Zm_exclude Pmin_exclude Pmax_exclude Krms_all_exclude pxy_exclude
for i = 1:size(Xm,2)
    Zm(:,i) = polyvaln(M,[Xm(:,i) Ym(:,i)]); 
end %i

 [K H Pmax Pmin] = surfature(Xm, Ym, Zm);
 Krms_all = ((Pmin.^2 + Pmax.^2)./2).^(1/2);


% figure(s) 
% A = mesh(Xm,Ym,Zm,Pmin); hold on;
[y I] = max(reshape(Pmin,numel(Pmin),1));
% plot3(Xm(I),Ym(I),Zm(I),'ko','MarkerFaceColor','g')

latDimpleL = [Xm(I) Ym(I), Zm(I)];

d = min(sum([(templocalpts(:,1) - latDimpleL(1)).^2 (templocalpts(:,2) - latDimpleL(2)).^2 ...
         (templocalpts(:,3) - latDimpleL(3)).^2],2).^(1/2));
    tempPmin = Pmin;
while d > 2
    tempPmin(I) = -1;
    [y I] = max(reshape(tempPmin,numel(tempPmin),1));
    latDimpleL = [Xm(I) Ym(I), Zm(I)];
    
    d = min(sum([(templocalpts(:,1) - latDimpleL(1)).^2 (templocalpts(:,2) - latDimpleL(2)).^2 ...
         (templocalpts(:,3) - latDimpleL(3)).^2],2).^(1/2));
        
end %check that it is not a false min

latDimpleG = RT_transform(latDimpleL,ACS_P_temp,Centroid,1);

tempVert = unit(cross(unit(latDimpleG - Centroid),ACS_L_P_z'));
tempLat = unit(cross(ACS_L_P_z',tempVert));

ACS_steer = [tempVert' tempLat' ACS_L_P_z Centroid'; [0 0 0 1]];

end %patellaIGuess
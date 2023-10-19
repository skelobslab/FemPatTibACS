function [patT] = patellaRefine(ivPath,ivFile,initACS,side)

global objf_out;



Xr = [0 0 0];
Yr = [0 0 0];
Zr = [0 0 0];

Centroid = initACS(1:3,4)';

[pts cnt] = read_vrml_fast(fullfile(ivPath,ivFile));


ACS_steer = RT_to_fX4([initACS(1:3,3) initACS(1:3,1) initACS(1:3,2)], initACS(1:3,4)');

lb = -15;
ub = 15;
gstd = [0];
options = optimset('Algorithm','interior-point','MaxFunEvals',100000,'TolX',1e-8,'Display','notify');
counter = 1;
[x,fval,exitflag]=fmincon(@evalRidgeDistances,5,[],[],[],[],lb,ub,[],options); 
Zrot = rotamatz(x,'deg');
patT = ACS_steer * Zrot;

temp = [patT(1:3,2) patT(1:3,3) patT(1:3,1)];
patT = RT_to_fX4(temp, patT(1:3,4)');



function objf = evalRidgeDistances(coordRot)
    samples = 100;
    zpercent = 0.20;
    zrot = rotamatz(coordRot,'deg');
    ACSrot = ACS_steer * zrot;
    steerPtsL = RT_transform(pts,ACSrot(1:3,1:3),ACSrot(1:3,4)',0);
 
    moveaxis_inc = linspace((0.55 * min(steerPtsL(:,1))),(0.88 * max(steerPtsL(:,1))),samples); 
    

    
    for i = 1:samples
        
         I_poi = find(steerPtsL(:,1) > moveaxis_inc(i) & steerPtsL(:,1) < (moveaxis_inc(i) + 2)   & steerPtsL(:,3) < (zpercent * min(steerPtsL(:,3))) );
        
        tempPts = steerPtsL(I_poi,:);
        [yy II] = min(tempPts(:,3));
        
        xdev(i,1) = tempPts(II,2);
    end %i
    
   
    
    objf = std(xdev);
%     objf = -abs(max(xdev)-min(xdev));
    gstd(counter,1) = objf;
    counter = counter + 1;
%     objf_out(s,1) = objf;
    obj_out = objf;
    
end %objective fucntion




end %main function


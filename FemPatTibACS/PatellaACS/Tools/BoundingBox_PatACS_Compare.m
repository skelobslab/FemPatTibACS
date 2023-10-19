%evaluate bounding box and patACS methods

subjects = {'3976','4147','5043','5519','6221','6580','7004','7497','8436','9446'};
sides = {'R','L','R','R','L','L','L','R','L','L'};

ivdir = 'P:\Patella\PatellaACS\models_processed';
pT = createInventorHeader();
pTR = createInventorHeader();
pTB = createInventorHeader();
subject = ['XJUC0' subjects{1}];
ivFile = [subject '_Pat.iv'];


Xr = [0 0 0]; Yr = [0 0 0]; Zr = [0 0 0];
XrR = [0 0 0]; YrR = [0 0 0]; ZrR = [0 0 0];
XrB = [0 0 0]; YrB = [0 0 0]; ZrB = [0 0 0];
XrI = [0 0 0]; YrI = [0 0 0]; ZrI = [0 0 0];

CR = [0 0 0];
CB = [0 0 0];


orange = [234/255 178/255 22/255];
blue = [67/255 67/255 255/255];

% s = 3;
for s = 1:length(subjects)

subject = ['XJUC0' subjects{s}];
S = sides{s};

 geoT = loadGeomagicTransform(fullfile(ivdir,[subject '.tfm']));
    geoT(1:3,1) = unit(geoT(1:3,1));
    geoT(1:3,2) = unit(geoT(1:3,2));
    geoT(1:3,3) = unit(geoT(1:3,3));
 geoTinv = transposefX4(geoT);
 
 [Centroid,SurfaceArea,Volume,CoM_ev123,CoM_eigenvectors,I1,I2,I_CoM,I_origin,patches] = mass_properties(fullfile(ivdir,ivFile));

%evaluate coordinate systems before and after refinement step
patT = patellaACS(ivdir,ivFile,'L');
patTR = patellaACSRefine(ivdir,ivFile,patT,'L');

patI = patellaACS_IGuess(ivdir,ivFile,'L');

[pts cnt] = read_vrml_fast(fullfile(ivdir,ivFile));

[box boxT] = boundingBox(pts);

patTg = geoT * patT;
patTRg = geoT * patTR;
boxTg = geoT * boxT;
patIg = geoT * patI;

Xr = unit(patTg(1:3,1)') + Xr;
Yr = unit(patTg(1:3,2)') + Yr;
Zr = unit(patTg(1:3,3)') + Zr;

XrR = unit(patTRg(1:3,1)') + XrR;
YrR = unit(patTRg(1:3,2)') + YrR;
ZrR = unit(patTRg(1:3,3)') + ZrR;

XrB = unit(patTRg(1:3,1)') + XrB;
YrB = unit(patTRg(1:3,2)') + YrB;
ZrB = unit(patTRg(1:3,3)') + ZrB;

XrI = unit(patIg(1:3,1)') + XrI;
YrI = unit(patIg(1:3,2)') + YrI;
ZrI = unit(patIg(1:3,3)') + ZrI;

XrAll(s,:) = patTg(1:3,1)';
YrAll(s,:) = patTg(1:3,2)';
ZrAll(s,:) = patTg(1:3,3)';

XrAllR(s,:) = patTRg(1:3,1)';
YrAllR(s,:) = patTRg(1:3,2)';
ZrAllR(s,:) = patTRg(1:3,3)';

XrAllB(s,:) = boxTg(1:3,1)';
YrAllB(s,:) = boxTg(1:3,2)';
ZrAllB(s,:) = boxTg(1:3,3)';

XrAllI(s,:) = patIg(1:3,1)';
YrAllI(s,:) = patIg(1:3,2)';
ZrAllI(s,:) = patIg(1:3,3)';

CR = [CR; patTRg(1:3,4)'];
CB = [CB: boxTg(1:3,4)'];

end %s

Xr = unit(Xr);
Yr = unit(Yr);
Zr = unit(Zr);
% 
XrR = unit(XrR);
YrR = unit(YrR);
ZrR = unit(ZrR);

XrB = unit(XrB);
YrB = unit(YrB);
ZrB = unit(ZrB);

XrI = unit(XrI);
YrI = unit(YrI);
ZrI = unit(ZrI);

for s = 1:length(subjects)
    ang_out(s,1) = rad2deg(acos(dot(Xr,XrAll(s,:))));
    ang_out(s,2) = rad2deg(acos(dot(Yr,YrAll(s,:))));
    ang_out(s,3) = rad2deg(acos(dot(Zr,ZrAll(s,:))));
    
    ang_outR(s,1) = rad2deg(acos(dot(XrR,XrAllR(s,:))));
    ang_outR(s,2) = rad2deg(acos(dot(YrR,YrAllR(s,:))));
    ang_outR(s,3) = rad2deg(acos(dot(ZrR,ZrAllR(s,:))));
    
    ang_outB(s,1) = rad2deg(acos(dot(XrB,XrAllB(s,:))));
    ang_outB(s,2) = rad2deg(acos(dot(YrB,YrAllB(s,:))));
    ang_outB(s,3) = rad2deg(acos(dot(ZrB,ZrAllB(s,:))));
    
    ang_outI(s,1) = rad2deg(acos(dot(XrI,XrAllI(s,:))));
    ang_outI(s,2) = rad2deg(acos(dot(YrI,YrAllI(s,:))));
    ang_outI(s,3) = rad2deg(acos(dot(ZrI,ZrAllI(s,:))));
    
end %s 

 

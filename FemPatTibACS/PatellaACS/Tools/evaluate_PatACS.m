%eval patella ACS
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

CR = [];
CB = [];

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
 
 
    

ivFile = [subject '_Pat.iv'];

if s == 1
    pTR = [pTR createInventorLink(fullfile(ivdir,ivFile),geoT(1:3,1:3),geoT(1:3,4)',[0.8 0.8 0.8], 0.55)];
    pTB = [pTB createInventorLink(fullfile(ivdir,ivFile),geoT(1:3,1:3),geoT(1:3,4)',[0.8 0.8 0.8], 0.55)];
end 

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

%end evaluation


% pT = [pT createInventorArrow(patTg(1:3,4)',patTg(1:3,1)',50,0.5,[1 0 0],0)];
% pT = [pT createInventorArrow(patTg(1:3,4)',patTg(1:3,2)',50,0.5,[0 1 0],0)];
% pT = [pT createInventorArrow(patTg(1:3,4)',patTg(1:3,3)',50,0.5,[0 0 1],0)];
pTRo = createInventorHeader();
pTRb = createInventorHeader();

%blue
pTRb = [pTRb createInventorLink(fullfile(ivdir,ivFile),geoT(1:3,1:3),geoT(1:3,4)',blue, 0.55)];
pTRb = [pTRb createInventorArrow(patTRg(1:3,4)',-patTRg(1:3,1)',40,1,[0.5 0.5 1],0)];
pTRb = [pTRb createInventorArrow(patTRg(1:3,4)',-patTRg(1:3,2)',40,1,[1 0.5 0.5],0)];
pTRb = [pTRb createInventorArrow(patTRg(1:3,4)',patTRg(1:3,3)',40,1,[0.5 1 0.5],0)];
pTRb = [pTRb createInventorSphere(patTRg(1:3,4)',2,[1 1 1],0)];

pTR = [pTR createInventorArrow(patTRg(1:3,4)',patTRg(1:3,1)',40,0.5,[1 0.5 0.5],0)];
pTR = [pTR createInventorArrow(patTRg(1:3,4)',patTRg(1:3,2)',40,0.5,[0.5 1 0.5],0)];
pTR = [pTR createInventorArrow(patTRg(1:3,4)',patTRg(1:3,3)',40,0.5,[0.5 0.5 1],0)];
pTR = [pTR createInventorSphere(patTRg(1:3,4)',2,[1 1 1],0)];

pTB = [pTB createInventorArrow(boxTg(1:3,4)',boxTg(1:3,1)',40,0.5,[1 0.5 0.5],0)];
pTB = [pTB createInventorArrow(boxTg(1:3,4)',boxTg(1:3,2)',40,0.5,[0.5 1 0.5],0)];
pTB = [pTB createInventorArrow(boxTg(1:3,4)',boxTg(1:3,3)',40,0.5,[0.5 0.5 1],0)];
pTB = [pTB createInventorSphere(boxTg(1:3,4)',2,[1 1 1],0)];


%orange
pTRo = [pTRo createInventorLink(fullfile(ivdir,ivFile),geoT(1:3,1:3),geoT(1:3,4)',orange, 0.55)];
pTRo = [pTRo createInventorArrow(patTRg(1:3,4)',-patTRg(1:3,1)',45,1,[0 0 1],0)];
pTRo = [pTRo createInventorArrow(patTRg(1:3,4)',-patTRg(1:3,2)',45,1,[1 0 0],0)];
pTRo = [pTRo createInventorArrow(patTRg(1:3,4)',patTRg(1:3,3)',45,1,[0 1 0],0)];
pTRo = [pTRo createInventorSphere(patTRg(1:3,4)',2,[0 0 0],0)];


fid = fopen(['K:\pTRgeo_orange_' num2str(s) '.iv'],'w');
fprintf(fid,pTRo);
fclose(fid);

fid = fopen(['K:\pTRgeo_blue_' num2str(s) '.iv'],'w');
fprintf(fid,pTRb);
fclose(fid);



ivstring = createInventorHeader();

% ivstring = [ivstring createInventorLink(fullfile(ivdir,ivFile),eye(3,3),zeros(1,3),orange, 0.55)];
% ivstring = [ivstring createInventorLink(fullfile(ivdir,ivFile),eye(3,3),zeros(1,3),blue, 0.55)];

% ivstring = [ivstring createInventorArrow(patT(1:3,4)',patT(1:3,1)',50,0.5,[1 0 0],0)];
% ivstring = [ivstring createInventorArrow(patT(1:3,4)',patT(1:3,2)',50,0.5,[0 1 0],0)];
% ivstring = [ivstring createInventorArrow(patT(1:3,4)',patT(1:3,3)',50,0.5,[0 0 1],0)];

% %orange
% ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-patTR(1:3,1)',40,1,[0 0 1],0.1)];
% ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,2)',40,1,[1 0 0],0.1)];
% ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,3)',40,1,[0 1 0],0.1)];
% ivstring = [ivstring createInventorSphere(patTR(1:3,4)',2,[0 0 0],0)];

% %blue
% ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-patTR(1:3,1)',40,1,[0.5 0.5 1],0.1)];
% ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,2)',40,1,[1 0.5 0.5],0.1)];
% ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,3)',40,1,[0.5 1 0.5],0.1)];
% ivstring = [ivstring createInventorSphere(patTR(1:3,4)',2,[1 1 1],0)];

ivstring = [ivstring createInventorLink(fullfile(ivdir,ivFile),eye(3,3),zeros(1,3),[0.8 0.8 0.8], 0.55)];
ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-patTR(1:3,1)',45,1,[0 0 1],0.0)];
ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-patTR(1:3,2)',45,1,[1 0 0],0.0)];
ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,3)',45,1,[0 1 0],0.0)];
ivstring = [ivstring createInventorSphere(patTR(1:3,4)',2,[0 0 0],0)];

ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-CoM_eigenvectors(1:3,2)',40,1.5,[0 0 1],0.0)];
ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-CoM_eigenvectors(1:3,1)',40,1.5,[1 0 0],0.0)];
ivstring = [ivstring createInventorArrow(patTR(1:3,4)',-CoM_eigenvectors(1:3,3)',40,1.5,[0 1 0],0.0)];

% fid = fopen(['K:\test_orange' num2str(s) '.iv'],'w');
fid = fopen(['K:\test_coord_inert' num2str(s) '.iv'],'w');
fprintf(fid,ivstring);
fclose(fid);

CR = [CR; patTRg(1:3,4)'];
CB = [CB; boxTg(1:3,4)'];

end %s subjects

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

pTR = [pTR createInventorArrow(patTRg(1:3,4)',XrR,50,1.2,[1 0 0],0)];
pTR = [pTR createInventorArrow(patTRg(1:3,4)',YrR,50,1.2,[0 1 0],0)];
pTR = [pTR createInventorArrow(patTRg(1:3,4)',ZrR,50,1.2,[0 0 1],0)];

pTB = [pTB createInventorArrow(boxTg(1:3,4)',XrB,50,1.2,[1 0 0],0)];
pTB = [pTB createInventorArrow(boxTg(1:3,4)',YrB,50,1.2,[0 1 0],0)];
pTB = [pTB createInventorArrow(boxTg(1:3,4)',ZrB,50,1.2,[0 0 1],0)];
pTB = [pTB createInventorSphere(boxTg(1:3,4)',2,[1 1 1],0)];

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

% fid = fopen('K:\pTgeo.iv','w');
% fprintf(fid,pT);
% fclose(fid);

fid = fopen('K:\pTRgeo_all.iv','w');
fprintf(fid,pTR);
fclose(fid);

fid = fopen('K:\pTBgeo_all.iv','w');
fprintf(fid,pTB);
fclose(fid);

% [mean(ang_out,1); std(ang_out,1)]
% [mean(ang_outR,1); std(ang_outR,1)]
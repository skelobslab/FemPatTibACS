%plot beak ridge

%eval patella ACS
subjects = {'3976','4147','5043','5519','6221','6580','7004','7497','8436','9446'};
sides = {'R','L','R','R','L','L','L','R','L','L'};

ivdir = 'K:\D_and_M_unlimited\PatellaACS\models_processed';
figure(1); 

for s = 1:length(subjects)
subject = ['XJUC0' subjects{s}];
S = sides{s};

ivFile = [subject '_Pat.iv'];

%evaluate coordinate systems before and after refinement step
patT = patellaACS(ivdir,ivFile);

[pts cnt] = read_vrml_fast(fullfile(ivdir,ivFile));

ptsL = RT_transform(pts,patT(1:3,1:3),patT(1:3,4)',0);
zpercent = 0.70;
I = find(abs(ptsL(:,2)) < 0.5 & ptsL(:,3) < zpercent * min(ptsL(:,3)));
poi = ptsL(I,:);



subplot(2,5,s); hold on;
plot(ptsL(I,1),ptsL(I,3),'k.');

end %s subjects


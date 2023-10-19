%quick contour check

subjects = {'3976','4147','5043','5519','6221','6580','7004','7497','8436','9446'};
sides = {'R','L','R','R','L','L','L','R','L','L'};

ivdir = 'P:\Patella\PatellaACS\models_processed';
for s = 1:10
subject = ['XJUC0' subjects{s}];
ivFile = [subject '_Pat.iv'];

patT = patellaACS(ivdir,ivFile,'L');
patTR = patellaACSRefine(ivdir,ivFile,patT,'L');

[pts cnt] = read_vrml_fast(fullfile(ivdir,ivFile));

ptsT = RT_transform(pts,patTR(1:3,1:3),patTR(1:3,4)',0);

ordered_contour = shape_anal_iv_051608_patella(ivdir, ivFile,'L',patT);

if ordered_contour(1,2) < ordered_contour(20,2)
    index = sort((1:1:length(ordered_contour))',1,'descend');
    ordered_contour = ordered_contour(index,:);
end

[Y I] = min(ordered_contour(:,3));
[Y Imx] = max(ordered_contour(:,3));

poi = ordered_contour(1:Imx,:);

poiF = buttbutt(200,3,poi,0);

I = find(abs(poiF(:,3)) <=7);
meanRidge = mean(poiF(I,2));
stdRidge = std(poiF(I,2));
I = find(poiF(:,2) > meanRidge + 1.8 * stdRidge);
oPts = poi(I,:);
oPts = oPts(find(oPts(:,3) < 0),:);
[y I] = max(oPts(:,3))
I = find(poi(:,2) == oPts(I,2) & poi(:,3) == oPts(I,3));


% camlight('right'); camlight('left'); % camlight('headlight');
% lighting gouraud;
% 
% cnt(:,1:3) = cnt(:,1:3) + 1;
% patch('faces',cnt(:,1:3),'vertices',ptsT,...
%         'facecolor',[200/255 250/255 255/255],...
%         'edgecolor','none','facealpha',0.5);
%     hold on; 
% 
% 
% 
% plot3(ordered_contour(:,1),ordered_contour(:,2),ordered_contour(:,3),'r.')

figure; subplot(2,1,1);  plot(ordered_contour(:,3),ordered_contour(:,2)); hold on;
plot(ordered_contour(1,3),ordered_contour(1,2),'ko','MarkerFaceColor','g')
plot(ordered_contour(15,3),ordered_contour(15,2),'ko','MarkerFaceColor','c')


plot(poi(:,3),poi(:,2),'g-','LineWidth',1.5);
plot(poi(I,3),poi(I,2),'ro');

plot(poiF(:,3),poiF(:,2),'m','LineWidth',1.3);
% plot(oc_neg(1,3),oc_neg(1,2),'k.')
% plot(oc_neg(end,3),oc_neg(end,2),'ro')
% plot(ordered_contour(I,3),ordered_contour(I,2),'ko','MarkerFaceColor','c')
title(subject)
% axis('equal')

% subplot(2,1,2); plot(poi(:,3), doc_neg);

end %s
tilefigs

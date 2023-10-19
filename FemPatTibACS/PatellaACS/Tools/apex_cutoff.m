function zcutoff = apex_cutoff(ivdir,ivFile,side, patTR)


[pts cnt] = read_vrml_fast(fullfile(ivdir,ivFile));

ptsT = RT_transform(pts,patTR(1:3,1:3),patTR(1:3,4)',0);

ordered_contour = shape_anal_iv_051608_patella(ivdir, ivFile,side,patTR);

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
[zcutoff I] = max(oPts(:,3))
% I = find(poi(:,2) == oPts(I,2) & poi(:,3) == oPts(I,3));




% camlight('right'); camlight('left'); % camlight('headlight');
% lighting gouraud;
% 
% cnt(:,1:3) = cnt(:,1:3) + 1;
% patch('faces',cnt(:,1:3),'vertices',ptsT,...
%         'facecolor',[200/255 250/255 255/255],...
%         'edgecolor','none','facealpha',0.5);
%     hold on; 

end %function


%generate, view, and save patella ACS
subjects = {'3976','4147','5043','5519','6221','6580','7004','7497','8436','9446'};
sides = {'R','L','R','R','L','L','L','R','L','L'};

ivdir = 'P:\Patella\PatellaACS\models_original';
% ivdir = 'P:\Patella\PatellaACS\models_processed';
savedir = 'P:\Patella\PatellaACS\ACS';

for s = 1:length(subjects)
% for s = 9

    subject = ['XJUC0' subjects{s}];
    ivFile = [subject '_Pat_' sides{s} '.iv'];
%     ivFile = [subject '_Pat.iv'];
    
    patT = patellaACS(ivdir,ivFile,sides{s});
    patTR = patellaACSRefine(ivdir,ivFile,patT,sides{s});
%     
%     patT = patellaACS(ivdir,ivFile,'L');
%     patTR = patellaACSRefine(ivdir,ivFile,patT,'L');
    
    dlmwrite(fullfile(savedir,[subject '_Pat_' sides{s} '_ACS.txt']),patTR,'\t');
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(ivdir,ivFile),eye(3,3),zeros(1,3),[0.8 0.8 0.8],0.3)];
    ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,1)',50,1,[1 0 0],0)];
    ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,2)',50,1,[0 1 0],0)];
    ivstring = [ivstring createInventorArrow(patTR(1:3,4)',patTR(1:3,3)',50,1,[0 0 1],0)];
    
    fid = fopen(['K:\TempViz\PatViz\' subject sides{s} '_patACSViz.iv'],'w');
    fprintf(fid,ivstring);
    fclose(fid);
end %s loop through subjects
function [box boxACS] = boundingBox(pts);
R = eye(3,3);
t = zeros(1,3);

pX = pts(:,1); pY = pts(:,2); pZ = pts(:,3);

box(1,:) = [max(pX) min(pY) max(pZ)];
box(2,:) = [max(pX) max(pY) max(pZ)];
box(3,:) = [min(pX) max(pY) max(pZ)];
box(4,:) = [min(pX) min(pY) max(pZ)];
box(5,:) = [max(pX) min(pY) min(pZ)];
box(6,:) = [max(pX) max(pY) min(pZ)];
box(7,:) = [min(pX) max(pY) min(pZ)];
box(8,:) = [min(pX) min(pY) min(pZ)];

R(1:3,1) = -(unit(box(5,:) - box(8,:))');
R(1:3,2) = -(unit(box(6,:) - box(5,:))');
R(1:3,3) = unit(box(1,:) - box(5,:))';

t = mean(box,1);
boxACS = RT_to_fX4(R,t);
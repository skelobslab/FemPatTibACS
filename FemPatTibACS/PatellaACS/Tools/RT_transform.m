function [d1] = RT_transform(d0,R,trans,direction);% Simple rigid body transformation of data.  Given R,T forward generates new% position, back assume a calulated R,T and allows superposition of data.% Direction back = 0, forward  = 1n_pts = size(d0,1);trans_mat = [trans(1)*ones(n_pts,1) ...             trans(2)*ones(n_pts,1) ...			 trans(3)*ones(n_pts,1)];			 if direction == 1, 	R_d = R*d0';	R_d = R_d';	d1 = R_d + trans_mat;end;if direction == 0, 	d0 = d0 - trans_mat;	d1 = R'*d0';	d1 = d1';end;
function pts_ordered = SortCrossSectionPts(pts)

count=1;
count2=1;
pts_ordered(1,:)=pts(count,:);
pt1=pts(1,:);
a=length(pts);
pts=pts(2:length(pts),:);

while length(pts)>1

    [row_pts col_pts]=size(pts);
    clear d
    for i=1:row_pts
        d(i,1)=dist(pt1,pts(i,:));
    end
    
    [min_d min_d_index]=min(d);
    
    while (min_d>5 & count2<a)
   
        pts=pts(find(pts(min_d_index,1)~=pts(:,1) & pts(min_d_index,2)~=pts(:,2)),:);
        
        [row_pts col_pts]=size(pts);
        clear d
        for i=1:row_pts
            d(i,1)=dist(pt1,pts(i,:));
        end
        
        if exist('d')
        [min_d min_d_index]=min(d);
        end
        count2=count2+1;
    end

    if count2<a
    count=count+1;count2=count2+1;
    pts_ordered(count,:)=pts(min_d_index,:);
    pt1=pts(min_d_index,:);
    pts=pts(find(pt1(:,1)~=pts(:,1) & pt1(:,2)~=pts(:,2)),:);
    end
end
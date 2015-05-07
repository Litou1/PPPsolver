function nodeCon=nodeNeighbors(mat,lv_mat,uniMesh)
[r ,c]=find(mat);% find nonzero nodes
nodeCon=zeros(size(r,1),20); % node connection matrix
for i=1:size(r,1)
    
    % search depth determined by the the level matrix
    % for uniform nodes, depth is 2^n, where n is the
    % sub-mesh level, otherwise it is twice the current level
    
    nodeLv=lv_mat(r(i),c(i));
    if nodeLv==0
        searchD=uniMesh;
    else
        searchD=2*uniMesh/(2^lv_mat(r(i),c(i)));
    end
    above_k=0;
    below_k=0;
    right_k=0;
    left_k=0;
    above_right_k=0;
    above_left_k=0;
    below_right_k=0;
    below_left_k=0;
    
    all_above=zeros(1,searchD);
    all_below=zeros(1,searchD);
    all_right=zeros(1,searchD);
    all_left=zeros(1,searchD);
    all_above_right=zeros(1,searchD);
    all_above_left=zeros(1,searchD);
    all_below_right=zeros(1,searchD);
    all_below_left=zeros(1,searchD);
    
    
    for j=1:searchD
        
        try
            all_above(j)=mat(r(i),c(i)+j);
        catch
            above=NaN;
            above_k=NaN;
        end
        
        try
            all_below(j)=mat(r(i),c(i)-j);
        catch
            below=NaN;
            below_k=NaN;
        end
        
        
        try
            all_left(j)=mat(r(i)-j,c(i));
        catch
            left=NaN;
            left_k=NaN;
        end
        
        try
            all_right(j)=mat(r(i)+j,c(i));
        catch
            right=NaN;
            right_k=NaN;
        end
        
        
        try
            all_above_right(j)=mat(r(i)+j,c(i)+j);
        catch
            above_right=NaN;
            above_right_k=NaN;
        end
        
        try
            all_above_left(j)=mat(r(i)-j,c(i)+j);
        catch
            above_left=NaN;
            above_left_k=NaN;
        end
        
        try
            all_below_right(j)=mat(r(i)+j,c(i)-j);
        catch
            below_right=NaN;
            below_right_k=NaN;
        end
        
        try
            all_below_left(j)=mat(r(i)-j,c(i)-j);
        catch
            below_left=NaN;
            below_left_k=NaN;
        end
        
    end
    
    % find the first non-zeros and "non-NaN" entry, that is the nearset node!
    % 0 means there is no node in that direction
    if ~isnan(above_k)
        above_k=find(all_above,1,'first');
        if ~isempty(above_k)
            above=all_above(above_k);
        else
            above_k=0;
            above=0;
        end
    end
    
    if ~isnan(below_k)
        below_k=find(all_below,1,'first');
        if ~isempty(below_k)
            below=all_below(below_k);
        else
            below_k=0;
            below=0;
        end
    end
    
    if ~isnan(left_k)
        left_k=find(all_left,1,'first');
        if ~isempty(left_k)
            left=all_left(left_k);
        else
            left_k=0;
            left=0;
        end
    end
    
    if ~isnan(right_k)
        right_k=find(all_right,1,'first');
        if ~isempty(right_k)
            right=all_right(right_k);
        else
            right_k=0;
            right=0;
        end
    end
    
    
    if ~isnan(above_right_k)
        above_right_k=find(all_above_right,1,'first');
        if ~isempty(above_right_k)
            above_right=all_above_right(above_right_k);
        else
            above_right_k=0;
            above_right=0;
        end
    end
    
    if ~isnan(above_left_k)
        above_left_k=find(all_above_left,1,'first');
        if ~isempty(above_left_k)
            above_left=all_above_left(above_left_k);
        else
            above_left_k=0;
            above_left=0;
        end
    end
    
    if ~isnan(below_right_k)
        below_right_k=find(all_below_right,1,'first');
        if ~isempty(below_right_k)
            below_right=all_below_right(below_right_k);
        else
            below_right_k=0;
            below_right=0;
        end
    end
    
    if ~isnan(below_left_k)
        below_left_k=find(all_below_left,1,'first');
        if ~isempty(below_left_k)
            below_left=all_below_left(below_left_k);
        else
            below_left_k=0;
            below_left=0;
        end
    end
    
    nodeCon(i,:)=[mat(r(i),c(i)) r(i) c(i) above below left right ...
        above_left above_right below_left below_right...
        above_k below_k left_k right_k ...
        above_left_k above_right_k below_left_k below_right_k...
        nodeLv];
end
nodeCon=sortrows(nodeCon);
end
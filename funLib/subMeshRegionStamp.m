function [lv_mat,mat]=subMeshRegionStamp...
    (num,k,subMeshLevel,factor,ctrNode,left,right,below,above,lv_mat,mat)
% This is a helper function for submesh region nodes stamping
for j=1:num
    % initialize the temporary level mattix
    tmp_lv_mat=zeros(size(mat));
    for i=1:subMeshLevel
        % number of the current meshes contained in each level in X and Y
        % direction
        numX=(k+left(j,i)+right(j,i))/(factor/2^i)+1;
        numY=(k+below(j,i)+above(j,i))/(factor/2^i)+1;
        % mark the positions of sub-mesh nodes for stamping
        fillMat=ones(numX,numY);
        %transform matrix T
        T=zeros(factor/2^i,factor/2^i);
        T(1,1)=1;
        fillMat=kron(fillMat, T);
        fillMat=fillMat(1:end-(factor/2^i-1),1:end-(factor/2^i-1));
        % the indices for the current level stamp
        hugeIdX=ctrNode(j,1)-left(j,i)-k/2:ctrNode(j,1)+right(j,i)+k/2;
        hugeIdY=ctrNode(j,2)-below(j,i)-k/2:ctrNode(j,2)+above(j,i)+k/2;
        mat(hugeIdX,hugeIdY)=mat(hugeIdX,hugeIdY) | fillMat;
        tmpIdx=zeros(size(mat));
        tmpIdx(hugeIdX,hugeIdY)=fillMat;
        tmp_lv_mat(logical(tmpIdx))=i;
    end
    lv_mat= max(lv_mat,tmp_lv_mat);
end
end
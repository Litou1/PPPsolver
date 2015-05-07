function [mat,lv_mat,shortNum,sourceNum,smallestMeshSize,uniMesh]=...
    meshGen(roundViaFLag,viewMesh,planeSizeX,planeSizeY,...
    short,source,viaSize,k,subMeshLevel,meshT)
% This function generates the mesh geomemtry matrix. Each voltage node has
% an assigned number. If there is no node, fill up with zeros.
%
% returns:
% mat, mesh geomemtry matrix
% lv_mat, submesh level matrix
% shortNum,sourceNum
% smallestMeshSize, in mm
% uniMesh, in mm
%
%  Note that the x direction and y direction are defined as
%  the following for the convienence of matching matrix indice
%  x-row y-column
%   .----------y
%   |
%   |
%   |
%   |
%   |
%   x

shortNum=size(short,1);
sourceNum=size(source,1);
% uniform mesh size is determined by the smallest via size and sub-mesh level
smallestMeshSize=viaSize/k
uniMeshSize=smallestMeshSize*2^subMeshLevel

% =========================================================================
%  assign uniform mesh nodes
% =========================================================================

%uniform node number in X and Y directions
uniNodeNumX=round(planeSizeX/uniMeshSize)+1;
uniNodeNumY=round(planeSizeY/uniMeshSize)+1;
uniMesh=2^subMeshLevel; % number of smallest mesh cell in one uniform cell
uniMesh=int64(uniMesh);
mat=ones(uniNodeNumX,uniNodeNumY);

% transform matrix T, a clever trick to fill with a lot zeros
T=zeros(uniMesh,uniMesh);
T(1,1)=1;
% make the uniform mesh sparse according the sub-mesh level
mat=kron(mat,T);
% submesh level matrix: same size as mat storing the value of submesh level for
% each nodes
lv_mat=zeros(size(mat));

% move shorts and sources to the nearest nodes
short = [ round(short(:,1)/smallestMeshSize)*(smallestMeshSize) ...
    round(short(:,2)/smallestMeshSize)*(smallestMeshSize)];
source = [ round(source(:,1)/smallestMeshSize)*(smallestMeshSize) ...
    round(source(:,2)/smallestMeshSize)*(smallestMeshSize)];

% find the center of the uniform cell at the via in the unit of smallest mesh
shortCtr  = int64(short/smallestMeshSize);
sourceCtr = int64(source/smallestMeshSize);
shortCtrNode  = int64(short/smallestMeshSize+1);
sourceCtrNode = int64(source/smallestMeshSize+1);

% throw the rightmost and bottommost edges
% they are extra columns and rows produced by the kronecker product function
mat=mat(1:end-(uniMesh-1),1:end-(uniMesh-1));
lv_mat=lv_mat(1:end-(uniMesh-1),1:end-(uniMesh-1));

if subMeshLevel~=0
    % =========================================================================
    % find the number of smallest mesh cells to get to the nearest uniform cell
    % "walls" in all directions
    % =========================================================================
    
    shortLeftK  = int64(mod(shortCtr(:,1)-k/2,uniMesh));
    shortRightK= int64(uniMesh-mod(shortCtr(:,1)+k/2,uniMesh));
    shortBelowK = int64(mod(shortCtr(:,2)-k/2,uniMesh));
    shortAboveK = int64(uniMesh-mod(shortCtr(:,2)+k/2,uniMesh));
    
    sourceLeftK  = int64(mod(sourceCtr(:,1)-k/2,uniMesh));
    sourceRightK = int64(uniMesh-mod(sourceCtr(:,1)+k/2,uniMesh));
    sourceBelowK = int64(mod(sourceCtr(:,2)-k/2,uniMesh));
    sourceAboveK = int64(uniMesh-mod(sourceCtr(:,2)+k/2,uniMesh));
    
%     Check if there are near the edges of uniform cells. Add an additional 
%     uniform mesh cell if K's are smaller than the tolerance, which means  
%     that it too close to the edge so more small mesh cells are needed
    edgeTol=k*4; % edge tolerance
    shortLeftK=int64(shortLeftK<=edgeTol)*uniMesh+shortLeftK;
    shortRightK=int64(shortRightK<=edgeTol)*uniMesh+shortRightK;
    shortBelowK=int64(shortBelowK<=edgeTol)*uniMesh+shortBelowK;
    shortAboveK=int64(shortAboveK<=edgeTol)*uniMesh+shortAboveK;
    
    sourceLeftK=int64(sourceLeftK<=edgeTol)*uniMesh+sourceLeftK;
    sourceRightK=int64(sourceRightK<=edgeTol)*uniMesh+sourceRightK;
    sourceBelowK=int64(sourceBelowK<=edgeTol)*uniMesh+sourceBelowK;
    sourceAboveK=int64(sourceAboveK<=edgeTol)*uniMesh+sourceAboveK;
    
    % =========================================================================
    %  Calculate the chain thickness of each level
    % =========================================================================

    % divided submesh into two kinds,
    % family I  is far away from the via
    % family II is near the via
    % family I  is sparser
    % family II is denser
    f1Lmt=1;% COULD BE 0 OR 1, 1 is denser
    % sub-mesh level  3 4 5 6 7
    % f2level         2 2 3 3 4
    f2level=floor((subMeshLevel+1)/2); % family II levels, typically the median 
    % of submesh levels, whose thickness can be no more than this famliy II limit
    f2Lmts=[4 4 8 12];% emprical limits, works good
    f2Lmt=f2Lmts(f2level);
    nthLmt=5; % decide the max of number of (n-1)th level meshes to move to
    % the smallest mesh if available
    % eg. previously, last level thickness was 1
    %     now, 1+3x2=7, which means it has to move 3 level n-1 submesh to level
    %     n

    % find the solution for thickness of each level -- a_i
    % d_to_uniform+meshT*2^n=a1*2^(n-1)+a2*2^(n-2)+....+a1*2+a0
    t_shortLeft=solThknssCoef(shortLeftK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    t_shortRight=solThknssCoef(shortRightK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    t_shortBelow=solThknssCoef(shortBelowK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    t_shortAbove=solThknssCoef(shortAboveK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    
    t_sourceLeft=solThknssCoef(sourceLeftK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    t_sourceRight=solThknssCoef(sourceRightK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    t_sourceBelow=solThknssCoef(sourceBelowK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    t_sourceAbove=solThknssCoef(sourceAboveK,subMeshLevel,meshT,f2level,f1Lmt,...
        f2Lmt,nthLmt);
    
    % convert thickness into reversed-order culmulative sum in the unit of the
    % smallest mesh
    
    % unit thickness array eg.32 16 8 4 2 1
    uT=repmat(2.^(subMeshLevel-1:-1:0),shortNum,1);
    short_offset_left  = fliplr(cumsum(double(fliplr(t_shortLeft.*uT)),2));
    short_offset_right = fliplr(cumsum(double(fliplr(t_shortRight.*uT)),2));
    short_offset_below = fliplr(cumsum(double(fliplr(t_shortBelow.*uT)),2));
    short_offset_above = fliplr(cumsum(double(fliplr(t_shortAbove.*uT)),2));
    
    uT=repmat(2.^(subMeshLevel-1:-1:0),sourceNum,1);
    source_offset_left  = fliplr(cumsum(double(fliplr(t_sourceLeft.*uT)),2));
    source_offset_right = fliplr(cumsum(double(fliplr(t_sourceRight.*uT)),2));
    source_offset_below = fliplr(cumsum(double(fliplr(t_sourceBelow.*uT)),2));
    source_offset_above = fliplr(cumsum(double(fliplr(t_sourceAbove.*uT)),2));
    
    % =========================================================================
    %  begin stamping sub-meshes for shorts and sources region
    % =========================================================================
    
    [lv_mat,mat]=subMeshRegionStamp(shortNum,k,subMeshLevel,uniMesh,...
        shortCtrNode,short_offset_left,short_offset_right,short_offset_below,...
        short_offset_above,lv_mat,mat);
    [lv_mat,mat]=subMeshRegionStamp(sourceNum,k,subMeshLevel,uniMesh,...
     sourceCtrNode,source_offset_left,source_offset_right,source_offset_below,...
     source_offset_above,lv_mat,mat);
end

% =========================================================================
%  assign shorts and source nodes
% =========================================================================
numVia=shortNum+sourceNum;
ctrViaNode=[shortCtrNode;sourceCtrNode];
mat=viaRegionStamp(roundViaFLag,numVia,k,ctrViaNode,mat);

% fill in, eg. 2 vias, 3 4 5 6... for logical indices (mat>0)
% all nodes expect for sources and shorts, which heave already taken care of
mat(mat>0)=numVia+1:numVia+sum(mat(:)>0);

% display the mesh geometry matrix
if viewMesh==1
    colormap(flipud(gray));
    figure(1)
    img1 = imagesc(rot90(mat));
    impixelregion(img1)
end
fprintf('Estimated unknowns: %d \n', (max(mat(:))-k^2*(shortNum+sourceNum))*3)
end



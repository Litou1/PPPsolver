function C=C_mat(option,nodes,h,smallestMeshSize,er,shortNum,sourceNum,k)
% This function generates C matrix, returns in pF
%
% viaBranchX and viaBranchY are branche logical indices contained in a via
% Examples:
% C=C_mat(nodes,h,smallestMeshSize,er,shortNum,sourceNum,k)
%
% $Author: Leihao Wei $    $Date:1/29/2014  $    $Revision: 0.1 $
% Copyright: Rose-Hulman Institute of Technology and Missouri S&T


eps=er*8.854e-12;% permittivity
sideLen=nodes.min_k;
C=zeros(size(nodes.thisNode,1),1);

if option==1 || option==2
    %type one caps
    type1=[nodes.oneTopLeft;
        nodes.oneBotLeft;
        nodes.oneTopRight;
        nodes.oneBotRight];
    type1Idx=ismember(nodes.thisNode,type1);
    C(type1Idx)=(1.5*sideLen(type1Idx)).^2;
    % type two caps
    type2=[nodes.twoTop;
        nodes.twoBot;
        nodes.twoLeft;
        nodes.twoRight;
        nodes.threeTop;
        nodes.threeBot;
        nodes.threeLeft;
        nodes.threeRight];
    type2Idx=ismember(nodes.thisNode,type2);
    C(type2Idx)=1.5*sideLen(type2Idx).^2;
    % type three caps
    type3=[nodes.fiveTopLeft;
        nodes.fiveTopRight;
        nodes.fiveBotLeft;
        nodes.fiveBotRight;
        nodes.fiveLeftAbove;
        nodes.fiveLeftBelow;
        nodes.fiveRightAbove;
        nodes.fiveRightBelow];
    type3Idx=ismember(nodes.thisNode,type3);
    C(type3Idx)=sideLen(type3Idx).^2+1.5*sideLen(type3Idx);
    % corner caps
    corner=[intersect(nodes.topBoard, nodes.leftBoard);
        intersect(nodes.topBoard, nodes.rightBoard);
        intersect(nodes.botBoard, nodes.leftBoard);
        intersect(nodes.botBoard, nodes.rightBoard)];
    cornerIdx=ismember(nodes.thisNode,corner);
    C(cornerIdx)=(0.5*sideLen(cornerIdx)).^2;
    % edge caps
    edge=setdiff([nodes.leftBoard;nodes.rightBoard;nodes.topBoard;nodes.botBoard],corner);
    edgeIdx=ismember(nodes.thisNode,edge);
    C(edgeIdx)=0.5*(sideLen(edgeIdx)).^2;
    % all others
    othersIdx=~(type1Idx | type2Idx | type3Idx | cornerIdx | edgeIdx);
    C(othersIdx)=(sideLen(othersIdx)).^2;
    
    % delete repeated via nodes
    [node_bk,ia,~]=unique(nodes.thisNode,'stable');
    C=C(ia);
    % shorting via cap
    shortVia=1:shortNum;
    shortViaIdx=ismember(node_bk,shortVia);
    C(shortViaIdx)=(2*k+1);
    % note that source via has 0 capacitance due to antipad
    sourceVia=shortNum+1:shortNum+sourceNum;
    sourceViaIdx=ismember(node_bk,sourceVia);
    C(sourceViaIdx)=0;
    
    % stamp the main diagonal terms
    C=1e12*sparse(node_bk,node_bk,C)*smallestMeshSize^2*1e-6/(h*1e-3)*eps;
else
    C=sparse(max(nodes.thisNode),max(nodes.thisNode));
end
end
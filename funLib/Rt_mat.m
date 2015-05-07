% generate R times t matrix, where t is the depth that current uniformly
% distributed in the cross section area, R is in the unit of kohms

function Rt=Rt_mat(option,branchX,branchY,viaBranchX,viaBranchY,sigma)
if option==2 || option==3
Rtx=branchX.size./(sigma*(abs(branchX.sy-branchX.ey)));
Rty=branchY.size./(sigma*(abs(branchY.sx-branchY.ex)));
else
    Rtx=zeros(length(branchX.size),1);
    Rty=zeros(length(branchY.size),1);
end

Rtx (viaBranchX) = [];
Rty (viaBranchY) = [];
numXbranch=size(Rtx,1);
numYbranch=size(Rty,1);
Rt=sparse(1:numXbranch+numYbranch,1:numXbranch+numYbranch,[Rtx;Rty]);

end
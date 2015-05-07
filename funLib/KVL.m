%*************************************************************************
%Generate the voltage submatrix part of MNA matrix based on KVL
%*************************************************************************

function KVL = KVL (nodes, branchX,branchY)

nodeCount=max(nodes.thisNode);
branchCount=max(branchX.thisBranch);

Vmx1 = sparse (branchX.thisBranch, branchX.startNode, -1,branchCount ,nodeCount);
Vmx2 = sparse (branchX.thisBranch, branchX.endNode, 1, branchCount,nodeCount);
Vmx = Vmx1 + Vmx2;


branchCount=max(branchY.thisBranch);

Vmy1 = sparse (branchY.thisBranch, branchY.startNode, -1,branchCount ,nodeCount);
Vmy2 = sparse (branchY.thisBranch, branchY.endNode, 1, branchCount,nodeCount);
Vmy = Vmy1 + Vmy2;

KVL = [Vmx;Vmy];

end
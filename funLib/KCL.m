%*************************************************************************
%Generate the current submatrix part of MNA matrix based on weight KCL
%*************************************************************************
function KCL = KCL (subMeshLevel, Vm,nodes,branchX,branchY)
Im=-Vm';
if subMeshLevel~=0
branchXCount=numel(branchX.thisBranch);
% 
% =========================================================================
% Type 1 nodes
% =========================================================================
% -------------------------------------------------------------------------
% 
% top left corner
% weight left branches to -0.75
[idx ,~]=find(bsxfun(@eq,branchX.endNode',nodes.oneTopLeft)');
tp1_topLeft_left=sub2ind(size(Im),nodes.oneTopLeft,branchX.thisBranch(idx));
Im(tp1_topLeft_left)=-0.75;
% weight above branches to 0.75
[idx ,~]=find(bsxfun(@eq,branchY.startNode',nodes.oneTopLeft)');
tp1_topLeft_above=sub2ind(size(Im),nodes.oneTopLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp1_topLeft_above)=0.75;

% ------------------------------------------------------------------------
% 
% top right corner
% 
% weight right branches to 0.75
[idx ,~]=find(bsxfun(@eq,branchX.startNode',nodes.oneTopRight)');
tp1_topRight_right=sub2ind(size(Im),nodes.oneTopRight,branchX.thisBranch(idx));
Im(tp1_topRight_right)=0.75;
% weight above branches to 0.75
[idx ,~]=find(bsxfun(@eq,branchY.startNode',nodes.oneTopRight)');
tp1_topRight_above=sub2ind(size(Im),nodes.oneTopRight,branchY.thisBranch(idx)+branchXCount);
Im(tp1_topRight_above)=0.75;


% ------------------------------------------------------------------------
% 
% bottom left corner
% 
% weight left branches to -0.75
[idx ,~]=find(bsxfun(@eq,branchX.endNode',nodes.oneBotLeft)');
tp1_botLeft_left=sub2ind(size(Im),nodes.oneBotLeft,branchX.thisBranch(idx));
Im(tp1_botLeft_left)=-0.75;
% weight below branches to -0.75
[idx ,~]=find(bsxfun(@eq,branchY.endNode',nodes.oneBotLeft)');
tp1_botLeft_below=sub2ind(size(Im),nodes.oneBotLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp1_botLeft_below)=-0.75;

% ------------------------------------------------------------------------
% 
% bottom right corner
% 
% weight right branches to 0.75
[idx ,~]=find(bsxfun(@eq,branchX.startNode',nodes.oneBotRight)');
tp1_botRight_right=sub2ind(size(Im),nodes.oneBotRight,branchX.thisBranch(idx));
Im(tp1_botRight_right)=0.75;
% weight below branches to -0.75
[idx ,~]=find(bsxfun(@eq,branchY.endNode',nodes.oneBotRight)');
tp1_botRight_below=sub2ind(size(Im),nodes.oneBotRight,branchY.thisBranch(idx)+branchXCount);
Im(tp1_botRight_below)=-0.75;

%%
% =========================================================================
% Type 2 nodes
% =========================================================================
% -------------------------------------------------------------------------
% left edges
% find its left above X branch and left below X branch
% weight left above X branches to -0.25
tp2_leftEdge_aboveNode=nodes.above(ismember(nodes.thisNode,nodes.twoLeft));
[idx, ~]=find(bsxfun(@eq,branchX.endNode', tp2_leftEdge_aboveNode)');
tp2_leftEdge_above=sub2ind(size(Im),nodes.twoLeft,branchX.thisBranch(idx));
Im(tp2_leftEdge_above)=-0.25;
% weight left below X branches to -0.25
tp2_leftEdge_belowNode=nodes.below(ismember(nodes.thisNode,nodes.twoLeft));
[idx, ~]=find(bsxfun(@eq,branchX.endNode', tp2_leftEdge_belowNode)');
tp2_leftEdge_below=sub2ind(size(Im),nodes.twoLeft,branchX.thisBranch(idx));
Im(tp2_leftEdge_below)=-0.25;

% -------------------------------------------------------------------------
% right edges
% find its right above X branch and right below X branch
% weight right above X branches to 0.25
tp2_rightEdge_aboveNode=nodes.above(ismember(nodes.thisNode,nodes.twoRight));
[idx, ~]=find(bsxfun(@eq,branchX.startNode', tp2_rightEdge_aboveNode)');
tp2_rightEdge_above=sub2ind(size(Im),nodes.twoRight,branchX.thisBranch(idx));
Im(tp2_rightEdge_above)=0.25;
% weight right below X branches to 0.25
tp2_rightEdge_belowNode=nodes.below(ismember(nodes.thisNode,nodes.twoRight));
[idx, ~]=find(bsxfun(@eq,branchX.startNode', tp2_rightEdge_belowNode)');
tp2_rightEdge_below=sub2ind(size(Im),nodes.twoRight,branchX.thisBranch(idx));
Im(tp2_rightEdge_below)=0.25;

% -------------------------------------------------------------------------
% top edges
% find its left above Y branch and right above Y branch
% weight left above Y branch to 0.25
tp2_topEdge_leftNode=nodes.left(ismember(nodes.thisNode,nodes.twoTop));
[idx, ~]=find(bsxfun(@eq,branchY.startNode', tp2_topEdge_leftNode)');
tp2_topEdge_left=sub2ind(size(Im),nodes.twoTop,branchY.thisBranch(idx)+branchXCount);
Im(tp2_topEdge_left)=0.25;
% weight right above Y branches to 0.25
tp2_topEdge_rightNode=nodes.right(ismember(nodes.thisNode,nodes.twoTop));
[idx, ~]=find(bsxfun(@eq,branchY.startNode', tp2_topEdge_rightNode)');
tp2_topEdge_right=sub2ind(size(Im),nodes.twoTop,branchY.thisBranch(idx)+branchXCount);
Im(tp2_topEdge_right)=0.25;

% -------------------------------------------------------------------------
% bottom edges
% find its left below Y branch and right below Y branch
% weight left below Y branch to -0.25
tp2_botEdge_leftNode=nodes.left(ismember(nodes.thisNode,nodes.twoBot));
[idx, ~]=find(bsxfun(@eq,branchY.endNode', tp2_botEdge_leftNode)');
tp2_botEdge_left=sub2ind(size(Im),nodes.twoBot,branchY.thisBranch(idx)+branchXCount);
Im(tp2_botEdge_left)=-0.25;
% weight right below Y branches to -0.25
tp2_botEdge_rightNode=nodes.right(ismember(nodes.thisNode,nodes.twoBot));
[idx, ~]=find(bsxfun(@eq,branchY.endNode', tp2_botEdge_rightNode)');
tp2_botEdge_right=sub2ind(size(Im),nodes.twoBot,branchY.thisBranch(idx)+branchXCount);
Im(tp2_botEdge_right)=-0.25;

%%
% =========================================================================
% Type 3 nodes
% =========================================================================
% -------------------------------------------------------------------------
% left edges
% weight left X branches to -0.5
[idx, ~]=find(bsxfun(@eq,branchX.endNode', nodes.threeLeft)');
tp3_leftEdge_left=sub2ind(size(Im),nodes.threeLeft,branchX.thisBranch(idx));
Im(tp3_leftEdge_left)=-0.5;
% -------------------------------------------------------------------------
% right edges
% weight right X branches to 0.5
[idx, ~]=find(bsxfun(@eq,branchX.startNode', nodes.threeRight)');
tp3_rightEdge_right=sub2ind(size(Im),nodes.threeRight,branchX.thisBranch(idx));
Im(tp3_rightEdge_right)=0.5;
% -------------------------------------------------------------------------
% top edges
% weight above Y branches to 0.5
[idx, ~]=find(bsxfun(@eq,branchY.startNode', nodes.threeTop)');
tp3_topEdge_above=sub2ind(size(Im),nodes.threeTop,branchY.thisBranch(idx)+branchXCount);
Im(tp3_topEdge_above)=0.5;
% -------------------------------------------------------------------------
% bottom edges
% weight below Y branches to -0.5
[idx, ~]=find(bsxfun(@eq,branchY.endNode', nodes.threeBot)');
tp3_botEdge_below=sub2ind(size(Im),nodes.threeBot,branchY.thisBranch(idx)+branchXCount);
Im(tp3_botEdge_below)=-0.5;


%%
% 
% =========================================================================
% Type 4 nodes
% =========================================================================
% -------------------------------------------------------------------------
% 
% top left corner
% weight right branches to 2/3
[idx ,~]=find(bsxfun(@eq,branchX.startNode',nodes.fourTopLeft)');
tp4_topLeft_right=sub2ind(size(Im),nodes.fourTopLeft,branchX.thisBranch(idx));
Im(tp4_topLeft_right)=2/3;
% weight below branches to -2/3
[idx ,~]=find(bsxfun(@eq,branchY.endNode',nodes.fourTopLeft)');
tp4_topLeft_below=sub2ind(size(Im),nodes.fourTopLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp4_topLeft_below)=-2/3;

% ------------------------------------------------------------------------
% 
% top right corner
% 
% weight left branches to -2/3
[idx ,~]=find(bsxfun(@eq,branchX.endNode',nodes.fourTopRight)');
tp4_topRight_left=sub2ind(size(Im),nodes.fourTopRight,branchX.thisBranch(idx));
Im(tp4_topRight_left)=-2/3;
% weight below branches to -2/3
[idx ,~]=find(bsxfun(@eq,branchY.endNode',nodes.fourTopRight)');
tp4_topRight_below=sub2ind(size(Im),nodes.fourTopRight,branchY.thisBranch(idx)+branchXCount);
Im(tp4_topRight_below)=-2/3;


% ------------------------------------------------------------------------
% 
% bottom left corner
% 
% weight right branches to 2/3
[idx ,~]=find(bsxfun(@eq,branchX.startNode',nodes.fourBotLeft)');
tp4_botLeft_right=sub2ind(size(Im),nodes.fourBotLeft,branchX.thisBranch(idx));
Im(tp4_botLeft_right)=2/3;
% weight above branches to 2/3
[idx ,~]=find(bsxfun(@eq,branchY.startNode',nodes.fourBotLeft)');
tp4_botLeft_above=sub2ind(size(Im),nodes.fourBotLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp4_botLeft_above)=2/3;

% ------------------------------------------------------------------------
% 
% bottom right corner
% 
% weight left branches to -2/3
[idx ,~]=find(bsxfun(@eq,branchX.endNode',nodes.fourBotRight)');
tp4_botRight_left=sub2ind(size(Im),nodes.fourBotRight,branchX.thisBranch(idx));
Im(tp4_botRight_left)=-2/3;
% weight above branches to 2/3
[idx ,~]=find(bsxfun(@eq,branchY.startNode',nodes.fourBotRight)');
tp4_botRight_above=sub2ind(size(Im),nodes.fourBotRight,branchY.thisBranch(idx)+branchXCount);
Im(tp4_botRight_above)=2/3;




%%
% =========================================================================
% Type 5 nodes
% =========================================================================
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% left above 
% find its right above X branch and left below X branch
% weight right above X branches to 2/3*0.5=1/3
tp5_leftAboveEdge_aboveNode=nodes.above(ismember(nodes.thisNode,nodes.fiveLeftAbove));
[idx, ~]=find(bsxfun(@eq,branchX.startNode', tp5_leftAboveEdge_aboveNode)');
tp5_leftAboveEdge_above=sub2ind(size(Im),nodes.fiveLeftAbove,branchX.thisBranch(idx));
Im(tp5_leftAboveEdge_above)=1/3;
% weight right below X branches to 0.25
tp5_leftAboveEdge_belowNode=nodes.below(ismember(nodes.thisNode,nodes.fiveLeftAbove));
[idx, ~]=find(bsxfun(@eq,branchX.startNode', tp5_leftAboveEdge_belowNode)');
tp5_leftAboveEdge_below=sub2ind(size(Im),nodes.fiveLeftAbove,branchX.thisBranch(idx));
Im(tp5_leftAboveEdge_below)=0.25;
% -------------------------------------------------------------------------
% left below 
% find its right above X branch and left below X branch
% weight right above X branches to 0.25
tp5_leftBelowEdge_aboveNode=nodes.above(ismember(nodes.thisNode,nodes.fiveLeftBelow));
[idx, ~]=find(bsxfun(@eq,branchX.startNode', tp5_leftBelowEdge_aboveNode)');
tp5_leftBelowEdge_above=sub2ind(size(Im),nodes.fiveLeftBelow,branchX.thisBranch(idx));
Im(tp5_leftBelowEdge_above)=0.25;
% weight right below X branches to 2/3*0.5=1/3
tp5_leftBelowEdge_belowNode=nodes.below(ismember(nodes.thisNode,nodes.fiveLeftBelow));
[idx, ~]=find(bsxfun(@eq,branchX.startNode', tp5_leftBelowEdge_belowNode)');
tp5_leftBelowEdge_below=sub2ind(size(Im),nodes.fiveLeftBelow,branchX.thisBranch(idx));
Im(tp5_leftBelowEdge_below)=1/3;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% right above 
% find its left above X branch and left below X branch
% weight left above X branches to -2/3*0.5=-1/3
tp5_rightAboveEdge_aboveNode=nodes.above(ismember(nodes.thisNode,nodes.fiveRightAbove));
[idx, ~]=find(bsxfun(@eq,branchX.endNode', tp5_rightAboveEdge_aboveNode)');
tp5_rightAboveEdge_above=sub2ind(size(Im),nodes.fiveRightAbove,branchX.thisBranch(idx));
Im(tp5_rightAboveEdge_above)=-1/3;
% weight left below X branches to -0.25
tp5_rightAboveEdge_belowNode=nodes.below(ismember(nodes.thisNode,nodes.fiveRightAbove));
[idx, ~]=find(bsxfun(@eq,branchX.endNode', tp5_rightAboveEdge_belowNode)');
tp5_rightAboveEdge_below=sub2ind(size(Im),nodes.fiveRightAbove,branchX.thisBranch(idx));
Im(tp5_rightAboveEdge_below)=-0.25;
% -------------------------------------------------------------------------
% right below 
% find its left above X branch and left below X branch
% weight left above X branches to -0.25
tp5_rightBelowEdge_aboveNode=nodes.above(ismember(nodes.thisNode,nodes.fiveRightBelow));
[idx, ~]=find(bsxfun(@eq,branchX.endNode', tp5_rightBelowEdge_aboveNode)');
tp5_rightBelowEdge_above=sub2ind(size(Im),nodes.fiveRightBelow,branchX.thisBranch(idx));
Im(tp5_rightBelowEdge_above)=-0.25;
% weight left below X branches to -2/3*0.5=1/3
tp5_rightBelowEdge_belowNode=nodes.below(ismember(nodes.thisNode,nodes.fiveRightBelow));
[idx, ~]=find(bsxfun(@eq,branchX.endNode', tp5_rightBelowEdge_belowNode)');
tp5_rightBelowEdge_below=sub2ind(size(Im),nodes.fiveRightBelow,branchX.thisBranch(idx));
Im(tp5_rightBelowEdge_below)=-1/3;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% top left
% find its left below Y branch and right below Y branch
% weight left below Y branch to -1/3
tp5_TopLeftEdge_leftNode=nodes.left(ismember(nodes.thisNode,nodes.fiveTopLeft));
[idx, ~]=find(bsxfun(@eq,branchY.endNode', tp5_TopLeftEdge_leftNode)');
tp5_TopLeftEdge_left=sub2ind(size(Im),nodes.fiveTopLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp5_TopLeftEdge_left)=-1/3;
% weight right below Y branches to -0.25
tp5_TopLeftEdge_rightNode=nodes.right(ismember(nodes.thisNode,nodes.fiveTopLeft));
[idx, ~]=find(bsxfun(@eq,branchY.endNode', tp5_TopLeftEdge_rightNode)');
tp5_TopLeftEdge_right=sub2ind(size(Im),nodes.fiveTopLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp5_TopLeftEdge_right)=-0.25;
% -------------------------------------------------------------------------
% top right
% find its left below Y branch and right below Y branch
% weight left below Y branch to -0.25
tp5_TopRightEdge_leftNode=nodes.left(ismember(nodes.thisNode,nodes.fiveTopRight));
[idx, ~]=find(bsxfun(@eq,branchY.endNode', tp5_TopRightEdge_leftNode)');
tp5_TopRightEdge_left=sub2ind(size(Im),nodes.fiveTopRight,branchY.thisBranch(idx)+branchXCount);
Im(tp5_TopRightEdge_left)=-0.25;
% weight right below Y branches to -1/3
tp5_TopRightEdge_rightNode=nodes.right(ismember(nodes.thisNode,nodes.fiveTopRight));
[idx, ~]=find(bsxfun(@eq,branchY.endNode', tp5_TopRightEdge_rightNode)');
tp5_TopRightEdge_right=sub2ind(size(Im),nodes.fiveTopRight,branchY.thisBranch(idx)+branchXCount);
Im(tp5_TopRightEdge_right)=-1/3;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% bot left
% find its left above Y branch and right above Y branch
% weight left above Y branch to 1/3
tp5_BotLeftEdge_leftNode=nodes.left(ismember(nodes.thisNode,nodes.fiveBotLeft));
[idx, ~]=find(bsxfun(@eq,branchY.startNode', tp5_BotLeftEdge_leftNode)');
tp5_BotLeftEdge_left=sub2ind(size(Im),nodes.fiveBotLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp5_BotLeftEdge_left)=1/3;
% weight right above Y branches to 0.25
tp5_BotLeftEdge_rightNode=nodes.right(ismember(nodes.thisNode,nodes.fiveBotLeft));
[idx, ~]=find(bsxfun(@eq,branchY.startNode', tp5_BotLeftEdge_rightNode)');
tp5_BotLeftEdge_right=sub2ind(size(Im),nodes.fiveBotLeft,branchY.thisBranch(idx)+branchXCount);
Im(tp5_BotLeftEdge_right)=0.25;
% -------------------------------------------------------------------------
% bot right
% find its left above Y branch and right above Y branch
% weight left above Y branch to 0.25
tp5_BotRightEdge_leftNode=nodes.left(ismember(nodes.thisNode,nodes.fiveBotRight));
[idx, ~]=find(bsxfun(@eq,branchY.startNode', tp5_BotRightEdge_leftNode)');
tp5_BotRightEdge_left=sub2ind(size(Im),nodes.fiveBotRight,branchY.thisBranch(idx)+branchXCount);
Im(tp5_BotRightEdge_left)=0.25;
% weight right above Y branches to 1/3
tp5_BotRightEdge_rightNode=nodes.right(ismember(nodes.thisNode,nodes.fiveBotRight));
[idx, ~]=find(bsxfun(@eq,branchY.startNode', tp5_BotRightEdge_rightNode)');
tp5_BotRightEdge_right=sub2ind(size(Im),nodes.fiveBotRight,branchY.thisBranch(idx)+branchXCount);
Im(tp5_BotRightEdge_right)=1/3;
end
%%
KCL=Im;
end
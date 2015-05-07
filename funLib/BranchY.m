   classdef BranchY
   % Branch is an object linking two nodes, much like a link list with
   % only two nodes!
   %  Note that the x direction and y direction are defined as 
   %  the following for the convenience of matching matrix indice
   %   x-row y-column
   %   .----------y
   %   |
   %   |
   %   |
   %   |
   %   |
   %   x

       properties
       % define the properties of the class here, (like fields of a struct)
           thisBranch
           startNode;
           endNode; %
           sx  % these are coordinate, [x y] 1 by 2 vector
           ex
           sy
           ey
           center
           node1 % coordinates of gaussian quadrature
           node2 % 
           size
           lv
       end
       methods
       % methods, including the constructor are defined in this block
           function obj = BranchY(nodes,smallestMeshSize)
           % class constructor
               if(nargin > 0)
                % delete starting nodes which represent the rightmost edge 
                 startIdx=~isnan(nodes.above) & nodes.above~=0;
                 obj.startNode    = nodes.thisNode(startIdx);
                 obj.thisBranch=(1:length(obj.startNode))';

                 obj.endNode      = nodes.above(startIdx);
                 obj.size=nodes.above_k(startIdx)*smallestMeshSize;
                 
                 obj.sy           = nodes.y(startIdx) ;
                 [endIdx ,~]=find(bsxfun(@eq,nodes.thisNode',obj.endNode)');
                 obj.ey           = nodes.y(endIdx);
                 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %new update 1.26.14
                 startLevel=nodes.nodeLv(startIdx);
                 endLevel=nodes.nodeLv(endIdx);
                 obj.lv=max([startLevel endLevel],[],2);
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 x= nodes.x(startIdx);
                 len=numel(obj.sy);
                 obj.sx=zeros(len,1);
                 obj.ex=zeros(len,1);
                 % pick out the left edge branches in transition region
                 % new update 12.14.13  
                 leftEdgePool=[nodes.oneBotLeft;nodes.twoLeft;nodes.threeLeft;...
                     nodes.fourTopRight;nodes.fourBotRight;nodes.fiveRightAbove;nodes.fiveRightBelow];
                 leftBranchIdx=ismember(obj.startNode,leftEdgePool);
                 
                 obj.sx(leftBranchIdx)           =  x(leftBranchIdx)-obj.size(leftBranchIdx);
                 obj.ex(leftBranchIdx)           =  x(leftBranchIdx)+obj.size(leftBranchIdx)/2;
                 
                 
                 % pick out the right edge branches in transition region
                 % new update 12.14.13
                 rightEdgePool=[nodes.oneBotRight;nodes.twoRight;nodes.threeRight;...
                     nodes.fourTopLeft;nodes.fourBotLeft;nodes.fiveLeftAbove;nodes.fiveLeftBelow];
                 rightBranchIdx=ismember(obj.startNode,rightEdgePool);

                 obj.sx(rightBranchIdx)           =  x(rightBranchIdx)-obj.size(rightBranchIdx)/2;
                 obj.ex(rightBranchIdx)           =  x(rightBranchIdx)+obj.size(rightBranchIdx);
                 
                 % eventually, pick out the left and right board edge branches
                 leftBoardBranchIdx=ismember(obj.startNode,nodes.leftBoard);
                 obj.sx(leftBoardBranchIdx)           =  x(leftBoardBranchIdx)-0;
                 obj.ex(leftBoardBranchIdx)           =  x(leftBoardBranchIdx)+obj.size(leftBoardBranchIdx)/2;                 
            
                 rightBoardBranchIdx=ismember(obj.startNode,nodes.rightBoard);
                 obj.sx(rightBoardBranchIdx)           =  x(rightBoardBranchIdx)-obj.size(rightBoardBranchIdx)/2;
                 obj.ex(rightBoardBranchIdx)           =  x(rightBoardBranchIdx)+0;
                 
                 % all other branches;
                 othersIdx=~(leftBranchIdx | rightBranchIdx | leftBoardBranchIdx | rightBoardBranchIdx);
                 obj.sx(othersIdx)           =  x(othersIdx)-obj.size(othersIdx)/2;
                 obj.ex(othersIdx)           =  x(othersIdx)+obj.size(othersIdx)/2;
                 
                 % Gaussian nodes position
                 W1=0.211;
                 W2=0.789;
                 obj.center=[(obj.sx+obj.ex)/2 (obj.sy+obj.ey)/2 ];
                 obj.node1 =[(obj.sx+obj.ex)/2-(W2-W1)/2*(obj.ey-obj.sy) (obj.sy+obj.ey)/2 ];
                 obj.node2 =[(obj.sx+obj.ex)/2+(W2-W1)/2*(obj.ey-obj.sy) (obj.sy+obj.ey)/2 ];
               end
           end
           
        function obj=assignViaBranchLocNum(obj,roundViaFLag,viaCount,k)
            
            for i=1:viaCount
                if k==4 && roundViaFLag==1
                    b = ((-i+1)*(k+1)^2-1:-1:-i*(k+1)^2+(k/2-1)*4)';
                else
                    b = ((-i+1)*(k+1)^2-1:-1:-i*(k+1)^2)'; % element that needs to be replaced
                end
                [idx, ~]=find(bsxfun(@eq,obj.startNode',b)');  
                obj.startNode (idx)=i;
                [idx, ~]=find(bsxfun(@eq,obj.endNode',b)');
                obj.endNode (idx)=i;                
            end
            
        end
       end
   end
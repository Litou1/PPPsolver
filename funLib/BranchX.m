   classdef BranchX
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
           size % section length
           lv   %branch level
       end
       methods
       % methods, including the constructor are defined in this block
           function obj = BranchX(nodes,smallestMeshSize)
           % class constructor
               if(nargin > 0)
                % delete starting nodes which represent the rightmost edge 
                 startIdx=~isnan(nodes.right) & nodes.right~=0;
                 obj.startNode    = nodes.thisNode(startIdx);
                 obj.thisBranch=(1:length(obj.startNode))';
                  
                 obj.endNode      = nodes.right(startIdx);
                 obj.size=nodes.right_k(startIdx)*smallestMeshSize;
                 
                 obj.sx           = nodes.x(startIdx) ;
                 [endIdx ,~]=find(bsxfun(@eq,nodes.thisNode',obj.endNode)');
                 obj.ex           =nodes.x(endIdx);
                                 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %new updates 1.26.14
                 startLevel=nodes.nodeLv(startIdx);
                 endLevel=nodes.nodeLv(endIdx);
                 obj.lv=max([startLevel endLevel],[],2);
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  
                 y= nodes.y(startIdx);
                 len=numel(obj.sx);
                 obj.sy=zeros(len,1);
                 obj.ey=zeros(len,1);
                 % pick out the top edge branches in transition region
                 % new updates 12.14.13
                 topEdgePool=[nodes.oneTopLeft;nodes.twoTop;nodes.threeTop;...
                     nodes.fourBotLeft;nodes.fourBotRight;nodes.fiveBotRight;nodes.fiveBotLeft];
                 topBranchIdx=ismember(obj.startNode,topEdgePool);
                 
                 obj.sy(topBranchIdx)           =  y(topBranchIdx)-obj.size(topBranchIdx)/2;
                 obj.ey(topBranchIdx)           =  y(topBranchIdx)+obj.size(topBranchIdx);
                 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % pick out the bottom edge branches in transition region
                 % new updates 12.14.13
                 bottomEdgePool=[nodes.oneBotLeft;nodes.twoBot;nodes.threeBot;...
                     nodes.fourTopLeft;nodes.fourTopRight;nodes.fiveTopRight;nodes.fiveTopLeft];
                 botBranchIdx=ismember(obj.startNode,bottomEdgePool);
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 obj.sy(botBranchIdx)           =  y(botBranchIdx)-obj.size(botBranchIdx);
                 obj.ey(botBranchIdx)           =  y(botBranchIdx)+obj.size(botBranchIdx)/2;
                 
                 % eventually, pick out the top and bottom board edge branches
                 topBoardBranchIdx=ismember(obj.startNode,nodes.topBoard);
                 obj.sy(topBoardBranchIdx)           =  y(topBoardBranchIdx)-obj.size(topBoardBranchIdx)/2;
                 obj.ey(topBoardBranchIdx)           =  y(topBoardBranchIdx)+0;                 
            
                 botBoardBranchIdx=ismember(obj.startNode,nodes.botBoard);
                 obj.sy(botBoardBranchIdx)           =  y(botBoardBranchIdx)-0;
                 obj.ey(botBoardBranchIdx)           =  y(botBoardBranchIdx)+obj.size(botBoardBranchIdx)/2;
                 
                 % all other branches;
                 othersIdx=~(topBranchIdx | botBranchIdx | topBoardBranchIdx | botBoardBranchIdx);
                 obj.sy(othersIdx)           =  y(othersIdx)-obj.size(othersIdx)/2;
                 obj.ey(othersIdx)           =  y(othersIdx)+obj.size(othersIdx)/2;
                 
                 
                 % Gaussian nodes position
                 W1=0.211;
                 W2=0.789;
                 obj.center=[(obj.sx+obj.ex)/2 (obj.sy+obj.ey)/2 ];
                 obj.node1 =[(obj.sx+obj.ex)/2 (obj.sy+obj.ey)/2+(W2-W1)/2*(obj.ey-obj.sy) ];
                 obj.node2 =[(obj.sx+obj.ex)/2 (obj.sy+obj.ey)/2-(W2-W1)/2*(obj.ey-obj.sy) ];
               end
           end
           
        % assign via branch local number   
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
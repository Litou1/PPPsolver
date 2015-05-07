classdef Node
    % write a description of the class here.
    properties
        % define the properties of the class here, ...
        % like fields Node(nodeCon,smallestMeshSize)of a struct)
        
        thisNode; % current node
        x;        
        y;
        
        % neighboring nodes 
        above; 
        below;
        left;
        right;
        above_left;
        above_right;
        below_left;
        below_right;
        
        % neighboring nodes distance -- k
        above_k;
        below_k;
        left_k;
        right_k;
        above_left_k;
        above_right_k;
        below_left_k;
        below_right_k;
        
        % identify types of nodes and its sub-category
        % e.g. Top left -- on top edge, left side
        
        oneTopLeft;
        oneBotLeft;
        oneTopRight;
        oneBotRight;
        
        twoTop;
        twoBot;
        twoLeft;
        twoRight;
        
        threeTop;
        threeBot;
        threeLeft;
        threeRight;
        
        fourTopLeft;
        fourBotLeft;
        fourTopRight;
        fourBotRight;
        
        fiveTopLeft;
        fiveTopRight;
        fiveBotLeft;
        fiveBotRight;
        fiveLeftAbove;
        fiveLeftBelow;
        fiveRightAbove;
        fiveRightBelow;
        
        % nodes at edges
        topBoard;
        botBoard;
        leftBoard;
        rightBoard;
        
        
        min_k;% for determining partial capacitor sizes
        
        % type 2 or type 5 nodes
        t2t3t5Top;
        t2t3t5Bot;
        t2t3t5Left;
        t2t3t5Right;
        
        nodeLv;% store info. of level of nodes
    end
    
    methods
        % methods, including the constructor are defined in this block
        function obj = Node(nodeCon,smallestMeshSize)
            % class constructor
            if(nargin > 0)
                obj.x     = (nodeCon(:,2)-1)*smallestMeshSize;
                obj.y     = (nodeCon(:,3)-1)*smallestMeshSize;
                
                obj.thisNode   = nodeCon(:,1);
                obj.above = nodeCon(:,4);
                obj.below = nodeCon(:,5);
                obj.left  = nodeCon(:,6);
                obj.right = nodeCon(:,7);
                obj.above_left = nodeCon(:,8);
                obj.above_right = nodeCon(:,9);
                obj.below_left  = nodeCon(:,10);
                obj.below_right = nodeCon(:,11);
                obj.above_k = nodeCon(:,12);
                obj.below_k = nodeCon(:,13);
                obj.left_k  = nodeCon(:,14);
                obj.right_k = nodeCon(:,15);
                obj.above_left_k = nodeCon(:,16);
                obj.above_right_k = nodeCon(:,17);
                obj.below_left_k  = nodeCon(:,18);
                obj.below_right_k = nodeCon(:,19);
                obj.nodeLv=nodeCon(:,20);
                %find type1 nodes at top left corner
                oneTopLeftIdx=(obj.above_left_k==2*obj.below_right_k)&(obj.above_right_k==2*obj.below_right_k)&...
                    (obj.below_left_k==2*obj.below_right_k);
                obj.oneTopLeft=obj.thisNode(oneTopLeftIdx);
                
                %find type1 nodes at bottom left corner
                oneBotLeftIdx=(obj.above_left_k==2*obj.above_right_k)&(obj.below_left_k==2*obj.above_right_k)&...
                    (obj.below_right_k==2*obj.above_right_k);
                obj.oneBotLeft=obj.thisNode(oneBotLeftIdx);
                
                %find type1 nodes at top right corner
                oneTopRightIdx=(obj.above_left_k==2*obj.below_left_k)&(obj.above_right_k==2*obj.below_left_k)&...
                    (obj.below_right_k==2*obj.below_left_k);
                obj.oneTopRight=obj.thisNode(oneTopRightIdx);
                
                %find type1 nodes at bottom right corner
                oneBotRightIdx=(obj.above_right_k==2*obj.above_left_k)&(obj.below_left_k==2*obj.above_left_k)&...
                    (obj.below_right_k==2*obj.above_left_k);
                obj.oneBotRight=obj.thisNode(oneBotRightIdx);
                
                
                
                %find type4 nodes at top left corner
                %                 fourTopLeftIdx=(2*obj.above_left_k==obj.below_right_k)&(2*obj.above_right_k==obj.below_right_k)&...
                %                     (2*obj.below_left_k==obj.below_right_k);
                tar=obj.below_right_k/2;
                fourTopLeftIdx=(obj.above_k==tar)&(obj.below_k==tar)&(obj.left_k==tar)&(obj.right_k==tar)&(obj.above_left_k==tar)&(obj.above_right_k==tar)&(obj.below_left_k==tar);
                obj.fourTopLeft=obj.thisNode(fourTopLeftIdx);
                
                %find type4 nodes at bottom left corner
                %                 fourBotLeftIdx=(2*obj.above_left_k==obj.above_right_k)&(2*obj.below_left_k==obj.above_right_k)&...
                %                     (2*obj.below_right_k==obj.above_right_k);
                tar=obj.above_right_k/2;
                fourBotLeftIdx=(obj.above_k==tar)&(obj.below_k==tar)&(obj.left_k==tar)&(obj.right_k==tar)&(obj.above_left_k==tar)&(obj.below_left_k==tar)&(obj.below_right_k==tar);
                obj.fourBotLeft=obj.thisNode(fourBotLeftIdx);
                
                %find type4 nodes at top right corner
                %                 fourTopRightIdx=(2*obj.above_left_k==obj.below_left_k)&(2*obj.above_right_k==obj.below_left_k)&...
                %                     (2*obj.below_right_k==obj.below_left_k);
                tar=obj.below_left_k/2;
                fourTopRightIdx=(obj.above_k==tar)&(obj.below_k==tar)&(obj.left_k==tar)&(obj.right_k==tar)&(obj.above_left_k==tar)&(obj.above_right_k==tar)&(obj.below_right_k==tar);
                obj.fourTopRight=obj.thisNode(fourTopRightIdx);
                
                %find type4 nodes at bottom right corner
                %                 fourBotRightIdx=(2*obj.above_right_k==obj.above_left_k)&(2*obj.below_left_k==obj.above_left_k)&...
                %                     (2*obj.below_right_k==obj.above_left_k);
                tar=obj.above_left_k/2;
                fourBotRightIdx=(obj.above_k==tar)&(obj.below_k==tar)&(obj.left_k==tar)&(obj.right_k==tar)&(obj.above_right_k==tar)&(obj.below_left_k==tar)&(obj.below_right_k==tar);
                obj.fourBotRight=obj.thisNode(fourBotRightIdx);
                
                
                % find type5 nodes on each edge. Each edge contains 2
                % types of nodes
                
                % top edge "left" type, "left" means it is near the left
                % corner
                %                 fiveTopLeftIdx=(obj.below_k==0)&(obj.below_right_k==0)&(obj.below_left_k==obj.above_left_k)&(obj.above_left_k==obj.above_right_k);
                fiveTopLeftIdx=(obj.below_k==0)&(obj.below_left_k==obj.above_left_k)&(obj.above_left_k==obj.above_right_k);
                obj.fiveTopLeft=obj.thisNode(fiveTopLeftIdx);
                
                %                 fiveTopRightIdx=(obj.below_k==0)&(obj.below_left_k==0)&(obj.below_right_k==obj.above_left_k)&(obj.above_left_k==obj.above_right_k);
                fiveTopRightIdx=(obj.below_k==0)&(obj.below_right_k==obj.above_left_k)&(obj.above_left_k==obj.above_right_k);
                obj.fiveTopRight=obj.thisNode(fiveTopRightIdx);
                
                % bottom edge "left" type, "left" means it is near the left
                % corner
                %                 fiveBotLeftIdx=(obj.above_k==0)&(obj.above_right_k==0)&(obj.above_left_k==obj.below_left_k)&(obj.below_left_k==obj.below_right_k);
                fiveBotLeftIdx=(obj.above_k==0)&(obj.above_left_k==obj.below_left_k)&(obj.below_left_k==obj.below_right_k);
                obj.fiveBotLeft=obj.thisNode(fiveBotLeftIdx);
                
                %                 fiveBotRightIdx=(obj.above_k==0)&(obj.above_left_k==0)&(obj.above_right_k==obj.below_right_k)&(obj.below_right_k==obj.below_left_k);
                fiveBotRightIdx=(obj.above_k==0)&(obj.above_right_k==obj.below_right_k)&(obj.below_right_k==obj.below_left_k);
                obj.fiveBotRight=obj.thisNode(fiveBotRightIdx);
                
                % left edge "above" type, "above" means it is near the
                % above corner
                %                 fiveLeftAboveIdx=(obj.right_k==0)&(obj.below_right_k==0)&(obj.above_right_k==obj.above_left_k)&(obj.above_left_k==obj.below_left_k);
                fiveLeftAboveIdx=(obj.right_k==0)&(obj.above_right_k==obj.above_left_k)&(obj.above_left_k==obj.below_left_k);
                obj.fiveLeftAbove=obj.thisNode(fiveLeftAboveIdx);
                
                %                 fiveLeftBelowIdx=(obj.right_k==0)&(obj.above_right_k==0)&(obj.below_right_k==obj.below_left_k)&(obj.below_left_k==obj.above_left_k);
                fiveLeftBelowIdx=(obj.right_k==0)&(obj.below_right_k==obj.below_left_k)&(obj.below_left_k==obj.above_left_k);
                obj.fiveLeftBelow=obj.thisNode(fiveLeftBelowIdx);
                
                % right edge "above" type, "above" means it is near the above
                % corner
                %                 fiveRightAboveIdx=(obj.left_k==0)&(obj.below_left_k==0)&(obj.above_left_k==obj.above_right_k)&(obj.above_right_k==obj.below_right_k);
                fiveRightAboveIdx=(obj.left_k==0)&(obj.above_left_k==obj.above_right_k)&(obj.above_right_k==obj.below_right_k);
                obj.fiveRightAbove=obj.thisNode(fiveRightAboveIdx);
                
                %                 fiveRightBelowIdx=(obj.left_k==0)&(obj.above_left_k==0)&(obj.below_left_k==obj.below_right_k)&(obj.below_right_k==obj.above_right_k);
                fiveRightBelowIdx=(obj.left_k==0)&(obj.below_left_k==obj.below_right_k)&(obj.below_right_k==obj.above_right_k);
                obj.fiveRightBelow=obj.thisNode(fiveRightBelowIdx);
                
                
                %----------------------------------------------------------------------------------------
                % type2 or type3 or type5
                %----------------------------------------------------------------------------------------
                
                t2t3TopIdx=(obj.left_k==obj.below_k)&(obj.below_k==obj.right_k)&((obj.above_k==2*obj.below_k)|obj.above_k==0);
                % get rid of possible type5 nodes
                %                 t2t3TopIdx=t2t3TopIdx&(~(fiveBotLeftIdx|fiveBotRightIdx));
                obj.t2t3t5Top=obj.thisNode(t2t3TopIdx);
                
                t2t3BotIdx=(obj.left_k==obj.above_k)&(obj.above_k==obj.right_k)&((obj.below_k==2*obj.above_k)|obj.below_k==0);
                % get rid of possible type5 nodes
                %                 t2t3BotIdx=t2t3BotIdx&(~(fiveTopLeftIdx|fiveTopRightIdx));
                obj.t2t3t5Bot=obj.thisNode(t2t3BotIdx);
                
                t2t3LeftIdx=(obj.right_k==obj.below_k)&(obj.below_k==obj.above_k)&((obj.left_k==2*obj.right_k)|obj.left_k==0);
                % get rid of possible type5 nodes
                %                 t2t3LeftIdx=t2t3LeftIdx&(~(fiveRightAboveIdx|fiveRightBelowIdx));
                obj.t2t3t5Left=obj.thisNode(t2t3LeftIdx);
                
                
                t2t3RightIdx=(obj.left_k==obj.below_k)&(obj.below_k==obj.above_k)&((obj.right_k==2*obj.left_k)|obj.right_k==0);
                % get rid of possible type5 nodes
                %                 t2t3RightIdx=t2t3RightIdx&(~(fiveLeftAboveIdx|fiveLeftBelowIdx));
                obj.t2t3t5Right=obj.thisNode(t2t3RightIdx);
                
                
                topBoardIdx=isnan(obj.above_k);
                obj.topBoard=obj.thisNode(topBoardIdx);
                
                botBoardIdx=isnan(obj.below_k);
                obj.botBoard=obj.thisNode(botBoardIdx);
                
                leftBoardIdx=isnan(obj.left_k);
                obj.leftBoard=obj.thisNode(leftBoardIdx);
                
                rightoardIdx=isnan(obj.right_k);
                obj.rightBoard=obj.thisNode(rightoardIdx);
                
                % determing the cap size
                kMat=nodeCon(:,12:15);
                tmp=kMat;
                tmp(tmp==0)=Inf;
                obj.min_k=min(tmp,[],2);
                
            end
        end
        
        
        function obj=assignViaLocNum(obj,roundViaFLag,viaCount,k)
            
            for i=1:viaCount
                if k==4 && roundViaFLag==1
                    b = ((-i+1)*(k+1)^2-1:-1:-i*(k+1)^2+(k/2-1)*4)'; % element that needs to be replaced
                else
                    b = ((-i+1)*(k+1)^2-1:-1:-i*(k+1)^2)'; % element that needs to be replaced
                end
                [idx, ~]=find(bsxfun(@eq,obj.thisNode',b)');
                obj.thisNode (idx)=i;
            end
            
        end
        
        
        function obj=assignT2T3Nodes(obj,conMat)
            
            
            [numRow,numCol]=size(conMat);
            
            refNodes=[obj.oneTopLeft;
                obj.oneBotLeft;
                obj.oneTopRight;
                obj.oneBotRight;
                obj.fourTopLeft;
                obj.fourBotLeft;
                obj.fourTopRight;
                obj.fourBotRight];
            fiveBot=[obj.fiveBotLeft;obj.fiveBotRight];
            fiveTop=[obj.fiveTopLeft;obj.fiveTopRight];
            fiveRight=[obj.fiveRightAbove;obj.fiveRightBelow];   
            fiveLeft=[obj.fiveLeftAbove;obj.fiveLeftBelow];
            
            % iterate through all connection matrix
            for i=1:numRow
                
                %get rid of zeros;
                getRow=conMat(i,:);
                getRow=getRow(getRow~=0);
                % split the row into many vectors
                rows=splitBasedOnRef(getRow,refNodes);
                
                for j=1:size(rows,2)
                    
                    thisRow=rows{j};
                    
                    % alternating elements
                    thisRowOdd=thisRow(1:2:end);
                    thisRowEven=thisRow(2:2:end-1);
                    
                    t2t3t5OddTopinThisRowIdx=ismember(thisRowOdd,obj.t2t3t5Top);
                    t2t3t5EvenTopinThisRowIdx=ismember(thisRowEven,obj.t2t3t5Top);
                    
                    t2t3t5OddBotinThisRowIdx=ismember(thisRowOdd,obj.t2t3t5Bot);
                    t2t3t5EvenBotinThisRowIdx=ismember(thisRowEven,obj.t2t3t5Bot);
                    
                    %check if it is all zero, if not go ahead
                    if ~all(~[t2t3t5OddTopinThisRowIdx t2t3t5EvenTopinThisRowIdx])
                        % get type2 or type3 in this row
                        % odd is type3 ; even is type 2
                        threeTopCandi=thisRowOdd(t2t3t5OddTopinThisRowIdx);
                        twoTopCandi=thisRowEven(t2t3t5EvenTopinThisRowIdx);
                        
                        % get rid of type5 nodes, get the real type2 and type3
                        twoTopReal=setdiff(twoTopCandi,fiveBot);
                        threeTopReal=setdiff(threeTopCandi,fiveBot);
                        % concatenate to obj
%                         obj.twoTop=sort(vertcat(obj.twoTop,twoTopReal'));
%                         obj.threeTop=sort(vertcat(obj.threeTop,threeTopReal'));

                        obj.twoTop=sort([obj.twoTop;twoTopReal']);
                        obj.threeTop=sort([obj.threeTop;threeTopReal']);

                    end
                    
                    if ~all(~[t2t3t5OddBotinThisRowIdx t2t3t5EvenBotinThisRowIdx])
                        % get type2 or type3 in this row
                        % odd is type3 ; even is type 2
                        threeBotCandi=thisRowOdd(t2t3t5OddBotinThisRowIdx);
                        twoBotCandi=thisRowEven(t2t3t5EvenBotinThisRowIdx);
                        
                        % get rid of type5 nodes, get the real type2 and type3
                        twoBotReal=setdiff(twoBotCandi,fiveTop);
                        threeBotReal=setdiff(threeBotCandi,fiveTop);
                        
%                         obj.twoBot=sort(vertcat(obj.twoBot,twoBotReal'));
%                         obj.threeBot=sort(vertcat(obj.threeBot,threeBotReal'));
                        
                        obj.twoBot=[obj.twoBot;twoBotReal'];
                        obj.threeBot=[obj.threeBot;threeBotReal'];
                        
                    end
                    
                end
            end
            
            for i=1:numCol
                
                
                %get rid of zeros;
                getCol=conMat(:,i);
                getCol=getCol(getCol~=0);
                % split the column into many vectors
                cols=splitBasedOnRef(getCol,refNodes);
                
                
                for j=1:size(cols,2)
                    
                    thisCol=cols{j};
                    
                    thisColOdd=thisCol(1:2:end);
                    thisColEven=thisCol(2:2:end-1);
                    
                    t2t3t5OddLeftinThisColIdx=ismember(thisColOdd,obj.t2t3t5Left);
                    t2t3t5EvenLeftinThisColIdx=ismember(thisColEven,obj.t2t3t5Left);
                    
                    t2t3t5OddRightinThisColIdx=ismember(thisColOdd,obj.t2t3t5Right);
                    t2t3t5EvenRightinThisColIdx=ismember(thisColEven,obj.t2t3t5Right);
                    
                    %check if it is all zero, if not go ahead
                    if ~all(~[t2t3t5OddLeftinThisColIdx;t2t3t5EvenLeftinThisColIdx])
                        
                        threeLeftCandi=thisColOdd(t2t3t5OddLeftinThisColIdx);
                        twoLeftCandi=thisColEven(t2t3t5EvenLeftinThisColIdx);
                        
                        twoLeftReal=setdiff(twoLeftCandi,fiveRight);
                        threeLeftReal=setdiff(threeLeftCandi,fiveRight);
                        
%                         obj.twoLeft=sort(vertcat(obj.twoLeft,twoLeftReal));
%                         obj.threeLeft=sort(vertcat(obj.threeLeft,threeLeftReal));
                        obj.twoLeft=[obj.twoLeft;twoLeftReal];
                        obj.threeLeft=[obj.threeLeft;threeLeftReal];
                    end
                    
                    if ~all(~[t2t3t5OddRightinThisColIdx;t2t3t5EvenRightinThisColIdx])
                        
                        threeRightCandi=thisColOdd(t2t3t5OddRightinThisColIdx);
                        twoRightCandi=thisColEven(t2t3t5EvenRightinThisColIdx);
                        
                        twoRightReal=setdiff(twoRightCandi,fiveLeft);
                        threeRightReal=setdiff(threeRightCandi,fiveLeft);
                        
%                         obj.twoRight=sort(vertcat(obj.twoRight,twoRightReal));
%                         obj.threeRight=sort(vertcat(obj.threeRight,threeRightReal));

                        obj.twoRight=[obj.twoRight;twoRightReal];
                        obj.threeRight=[obj.threeRight;threeRightReal];
                    end
                    
                end
            end
            
            
            obj.twoTop=sort(obj.twoTop);
            obj.twoBot=sort(obj.twoBot);            
            obj.twoLeft=sort(obj.twoLeft);
            obj.twoRight=sort(obj.twoRight);
            obj.threeTop=sort(obj.threeTop);
            obj.threeBot=sort(obj.threeBot);            
            obj.threeLeft=sort(obj.threeLeft);
            obj.threeRight=sort(obj.threeRight);

        end
    end
    
end
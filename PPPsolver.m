function [Z,arg]=PPPsolver(Domain,SolverOptions,TDsettings,ESLflag,roundViaFlag,...
    viewMesh,matprty,meshSettings,f)
%  Note that x direction and y direction are defined as
%  the following parttern for the convienence of matching matrix indice
%  x-row y-column
%   .----------y
%   |
%   |
%   |
%   |
%   |
%   x

% SolverOptions=0 : L only
% SolverOptions=1 : L and C
% SolverOptions=2 : L R C G, set to this for TD
% SolverOptions=3 : R only

% define copper thicness, conductivity, relative permittivity, and loss tagent
% all length are in mm
tc=matprty.tc;%
sigma=matprty.sigma;
er=matprty.er;
tan=matprty.tan;
%
planeSizeX=matprty.planeSizeX;
planeSizeY=matprty.planeSizeY;
h=matprty.h;
viaSize=matprty.viaSize;
short=matprty.short;
source=matprty.source;
% number of nodes per side length of a via, e.g.. 2,3,4
k=meshSettings.k;
subMeshLevel=meshSettings.subMeshLevel;
% T=1, submesh region = 3x3 uniform mesh length 
% T=2, submesh region = 5x5 uniform mesh length 
meshT=meshSettings.meshT;

% set a dummy f for TD simulation and dc solution
if strcmp(Domain,'TD') || SolverOptions==0 ; f=1;end  
if SolverOptions==2; s=sqrt(-1)*2*pi*f;end       % let s=jw
if SolverOptions==0||SolverOptions==3 ; s=1; end % dc solution
if SolverOptions==1; s=sqrt(-1)*2*pi*f;tan=0;end % L and C solution

FreqSweepFlag=0;
if nargin>8
    if numel(f)>1
        FreqSweepFlag=1; % do freq sweep if length of f>1
    end
end

%%
% =========================================================================
% Generate mesh 
% =========================================================================
tic;
t1=toc;
% roundViaFlag -- a rounded via model
% viewMatFlag -- see mesh matrix
[mat,lv_mat,shortNum,sourceNum,smallestMeshSize,uniMesh]=...
    meshGen(roundViaFlag,viewMesh,planeSizeX,planeSizeY,...
    short,source,viaSize,k,subMeshLevel,meshT);
% retrieve node connection information
nodeCon=nodeNeighbors(mat,lv_mat,uniMesh);
t2=toc;
fprintf('Meshing Time: %2.2f sec \n', t2-t1);
%%
% =========================================================================
%  Construct the node object and branch objects
% =========================================================================
t1=toc;
% number of vias
numVia=shortNum+sourceNum;
nodes=Node(nodeCon,smallestMeshSize);
% assign type2 and type3 nodes
nodes=nodes.assignT2T3Nodes(rot90(mat));
branchX=BranchX(nodes,smallestMeshSize);
branchY=BranchY(nodes,smallestMeshSize);

% =========================================================================
% save data for plotting latter
savefile='./plotData/nodesAndBranches.mat';
save(savefile,'nodes','branchX', 'branchY')
% =========================================================================

% assign 11111 2222 3333....to vias
nodes=nodes.assignViaLocNum(roundViaFlag,numVia,k);
branchX=branchX.assignViaBranchLocNum(roundViaFlag,numVia,k);
branchY=branchY.assignViaBranchLocNum(roundViaFlag,numVia,k);

t2=toc;
fprintf('Generate nodes and branches : %2.2f sec \n', t2-t1);
%%
% =========================================================================
%  Construct L and C matrices. Lossy terms -- R and G are included if
%  necessary
% =========================================================================
C=C_mat(SolverOptions,nodes,h,smallestMeshSize,er,shortNum,sourceNum,k);

%find via branches
viaBranchX=branchX.startNode==branchX.endNode;
viaBranchY=branchY.startNode==branchY.endNode;

t1=toc;

[L, numXbranch, numYbranch]=L_mat(SolverOptions,branchX,branchY,viaBranchX,viaBranchY,h);
L=L+tril(L,-1)';

t2=toc;
fprintf('Generate L matrix : %2.2f sec \n', t2-t1);
% =========================================================================
% R matrix and G matrix
% =========================================================================
% R is frequency dependent. Let's do R*t here.
Rt=Rt_mat(SolverOptions,branchX,branchY,viaBranchX,viaBranchY,sigma); 
% add the loss tangent to C matrix which is G matrix, note that -1i will be
% changed to 1 after multiplying s=jw !
C_and_G=C*(1-sqrt(-1)*tan);
%%
% =========================================================================
%  Construct the  connection matrix
% =========================================================================
% Connection matrix -- A
connMatA=KVL(nodes,branchX,branchY);
% Transpose of weighted connection matrix -- A^(t)
connMatA_T=KCL (subMeshLevel ,connMatA,nodes,branchX,branchY);
allViaBranch=[viaBranchX ; viaBranchY];
connMatA(allViaBranch,:)= []; % delete the branches representing vias
connMatA_T(:,allViaBranch)= [];
%%
% =========================================================================
% Construct voltage, current reference vectors, current injecion
% vector and Zc matrix
% =========================================================================
numNodes=max(nodes.thisNode);
numBranch=numXbranch+numYbranch;
%--------------------------------------------------------------------------
% Define the element in the current and voltage matrix for shorts and
% sources
%--------------------------------------------------------------------------
I_ref = zeros (numNodes+sourceNum+numBranch, numVia) ;
for i = 1 : numVia
    if i<=shortNum
        I_ref ( i, i ) = 1 ;
    else
        j=i-shortNum;
        I_ref(j+numNodes,j+shortNum)=1;
    end
end
V_ref = I_ref';
%--------------------------------------------------------------------------
% Define source current with the value of 1, Ipi as the value of -1
%--------------------------------------------------------------------------
Is = zeros ( numVia , 1 );
for i = 1 : sourceNum
    Is ( shortNum + i ) = 1 ;
    Is (numNodes+i)=-1;
end
%--------------------------------------------------------------------------
% Zc matrix
%--------------------------------------------------------------------------
Zc=Zc_mat(ESLflag,viaSize,h,short,source,shortNum,sourceNum,numVia);
%--------------------------------------------------------------------------
% source vector--b
%--------------------------------------------------------------------------
b = 1e3*sparse([Is;zeros(numBranch+shortNum+sourceNum,1)]);%mA

fprintf('Number of unknowns to be solved: %d \n', length(b));
%--------------------------------------------------------------------------
% begin solving AX=b
%--------------------------------------------------------------------------
if strcmp(Domain,'FD')
    if ~FreqSweepFlag
        [unkX,A]=oneFreqSolver(s,C_and_G,L,Rt,connMatA,connMatA_T,I_ref,...
            V_ref,b,numBranch,sourceNum,Zc,numNodes,sigma,tc);
        [Z,arg]=dispAns(SolverOptions,ESLflag,unkX,numVia,numNodes,sourceNum);
   
        fprintf('Sparsity of MNA matrix = %3.3f%% \n',nnz(A)/numel(A)*100);
        
        I = abs(unkX(numNodes+sourceNum+1:end-numVia)); % extract current array
        % save results for plotting
        savefile='./plotData/currentPlotting.mat';
        save(savefile, 'viaBranchX', 'viaBranchY', 'I',...
            'numXbranch','numYbranch','planeSizeX','planeSizeY');
    else     
        Z=zeros(numel(f),1);
        arg=zeros(numel(f),1);
        for i=1:length(f)
            fprintf('f= %f MHz, %d out of %d frequency samples\n',1e3*f(i),i,length(f));
            [unkX,A]=oneFreqSolver(s(i),C_and_G,L,Rt,connMatA,connMatA_T,I_ref,...
                V_ref,b,numBranch,sourceNum,Zc,numNodes,sigma,tc);
            [Z(i),arg(i)]=dispAns(SolverOptions,ESLflag,unkX,numVia,numNodes,sourceNum);
            fprintf('\n');         
            if i==1
                fprintf('Sparsity of MNA matrix = %3.3f%% \n',nnz(A)/numel(A)*100);
            end       
        end
    end
elseif strcmp(Domain,'TD')
    savefile='./plotData/TDCurrent.mat';
    save(savefile, 'viaBranchX', 'viaBranchY', ...
    'numXbranch','numYbranch','planeSizeX','planeSizeY');
%-------------?o-------------------------------------------------------------
% For time domain solver
%--------------------------------------------------------------------------    
    connMatA=-connMatA;
    R=2*Rt/tc;    

    unkX = TDsolver(TDsettings,C,L,R,connMatA,connMatA_T,...
    numBranch,sourceNum,shortNum,numVia,Zc,numNodes);
    
    %  return the unknown vector, 3 dimemsion
    Z=unkX;
    arg=[];

end



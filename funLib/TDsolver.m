function unkX = TDsolver(TDsettings,capMat,indMat,resMat,connMatA,connMatA_T,...
    numBranch,sourceNum,shortNum,numVia,Zc,numNodes)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *************************Desription**************************************
% unkX = unknown vector( Voltages and Currents ) = size: (m+n)x1x1
% varargin = unkX0 = initial condition of unknown vector(default is zeros)
% capMat = C matrix
% condMat = G matrix
% indMat = L matrix
% resMat = R matrix
% connMatA = Connectivity matrix
% xAxis = time axis for TD and Frequency axis
% b = Source matrix
% ***********Time Domain Formulation of Solution in Matricx Form (2.25) ********
%   [C/hp+G  A'   ] * unkX(:,1,p) =  [C/hp    0 ] [unkX(:,1,p-1)+unkX(:,1,p-2)]..
%   [  A   -L/hp-R]                  [ 0   -L/hp]
%                                  + b(:,1,p)
%  [C/hp    0 ]
%  [ 0   -L/hp] is the LC scaling matrix for RHS
% *************************************************************************
%%
ESR=TDsettings.ESR;
ESL=TDsettings.ESL;
Cd=TDsettings.Cd;
Rs=TDsettings.Rs;
sourceType=TDsettings.sourceType;
tr=TDsettings.riseTime;
tEnd=TDsettings.tEnd;
Nsteps=TDsettings.Nsteps;
width=TDsettings.width;
%%
tic
if isempty(capMat)
    capMat = 0;
end
if isempty(indMat)
    indMat = 0;
end

if isempty(resMat)
    resMat = 0;
end

condMat=0;

% ESL and ESR matrices
ESL_mat=zeros(size(Zc));
ESR_mat=zeros(size(Zc));
for i=1:shortNum
    ESL_mat(i,i)=ESL;
    ESR_mat(i,i)=ESR;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time domain Solver
% define sourceVector 
t=linspace(0,tEnd*1e-9,Nsteps);%sec
if strcmp(sourceType,'voltage')
    % unit step function
    sourceVector=zeros(numNodes+numBranch+3*numVia,1,Nsteps);
    stepSource(1,1,:)=stepFuction(t,tr,1);
    sourceVector(end,1,:)=stepSource;
elseif strcmp(sourceType,'current')
    triSource(1,1,:)=triPulse(t,tr,width,1);
    sourceVector=zeros(numNodes+numBranch+3*numVia,1,Nsteps);
    sourceVector(shortNum+1,1,:)=triSource;
    sourceVector(end-numVia,1,:)=-triSource;
    sourceVector= sourceVector/1e-3;% [Is(mA) Vs]= [Is(A)/1e-3 V]
end
b=sourceVector;
%%
xAxis=t;
xAxis= xAxis/1e-9;%ns
[numRow,~,numTimeSamples] = size(b);
numCaps = size(capMat,1);
numInds = size(indMat,1);
timeX = xAxis;
h(2:numTimeSamples) = timeX(2:numTimeSamples)-timeX(1:numTimeSamples-1);
h(1)=h(2);
unkX = zeros(numRow,1,numTimeSamples);
% LC scaling matrix for RHS
scaLC=sparse(numRow,numRow);
V_ref=scaLC;
scaLC(1:numNodes,1:numNodes)=-(capMat);
scaLC(numNodes+1:numNodes+numBranch...
    ,numNodes+1:numNodes+numBranch)=indMat;
scaLC(numNodes+numBranch+1:numNodes+numBranch+size(Zc,1)...
    ,numNodes+numBranch+1:numNodes+numBranch+size(Zc,1))=(ESL_mat+Zc);
spaCd=sparse([numRow-numVia+1:numRow-numVia+shortNum numRow-numVia+1:numRow-numVia+shortNum]...
    ,[1:shortNum (1:shortNum)+numNodes+numBranch+size(Zc,1)],...
    [Cd*ones(shortNum,1) -Cd*ones(shortNum,1)],numRow,numRow);
scaLC=scaLC+(-spaCd);

% define BD2 method coefficients
kp=1.5/h(1);
kpMins1=-2/h(1);
kpMins2=0.5/h(1);

MNA_BD2 =  [ (kp*capMat+condMat), connMatA_T;
    connMatA        ,-(resMat+kp*indMat)];

MNA_BD2=[MNA_BD2                        ,  zeros(numCaps+numInds,size(Zc,1));
    zeros(size(Zc,1),numCaps+numInds), (-ESR_mat-(ESL_mat+Zc)*kp)];

MNA_BD2=[MNA_BD2                      ,  zeros(size(MNA_BD2,1),numVia*2);
    zeros(numVia*2,size(MNA_BD2,1)),  zeros(numVia*2,numVia*2)];

MNA_BD2=MNA_BD2+sparse(numRow-numVia+1:numRow-sourceNum,numRow-numVia+1:numRow-sourceNum,-1,numRow,numRow);
MNA_BD2(end,end)=-Rs;

for i=1:numVia
    V_ref(end-numVia+i,i)=1;
    V_ref(end-numVia+i,end-numVia*2+i)=-1;
    V_ref(end-numVia*2+i,end-numVia*2-size(Zc,1)+i)=1;
end

% change some entries in V_ref to Cd
for i=1:shortNum
    V_ref(end-numVia+i,i)=kp*Cd;
    V_ref(end-numVia+i,end-numVia*2+i)=-kp*Cd;
end

A_BD2=MNA_BD2+V_ref+V_ref';

if strcmp(sourceType,'current')
    A_BD2(end,end)=0; % delete Rs due to current excitation;
    A_BD2(end-numVia,end)=0;
end

% initial two steps
unkX0 = zeros(size(b,1),1);
unkX(:,1,1) = A_BD2\(scaLC*(kpMins1*unkX0(:)   +0) + sparse(b(:,1,1)));
unkX(:,1,2) = A_BD2\(scaLC*(kpMins1*unkX(:,1,1)+kpMins2*unkX0(:)) + sparse(b(:,1,2)));

eps=1e-6;% a small number
for itime = 3:numTimeSamples
    
    fprintf('Iteration = %d of %d \n',itime,numTimeSamples);
    
    if abs(h(itime)-h(itime-1))>eps
        
        fprintf('take a new step size!')
        
        % update coefficients
        kp=1.5/h(itime);
        kpMins1=-2/h(itime);
        kpMins2=0.5/h(itime);
        
        MNA_BD2 =  [ (kp*capMat+condMat), connMatA_T;
            connMatA        ,-(resMat+kp*indMat)];
        
        MNA_BD2=[MNA_BD2                        ,  zeros(numCaps+numInds,size(Zc,1));
            zeros(size(Zc,1),numCaps+numInds), (-ESR_mat-(ESL_mat+Zc)*kp)];
        
        MNA_BD2=[MNA_BD2                      ,  zeros(size(MNA_BD2,1),numVia*2);
            zeros(numVia*2,size(MNA_BD2,1)),  zeros(numVia*2,numVia*2)];
        
        MNA_BD2=MNA_BD2+sparse(numRow-numVia+1:numRow-sourceNum,numRow-numVia+1:numRow-sourceNum,-1,numRow,numRow);
        MNA_BD2(end,end)=-Rs;
        
        for i=1:numVia
            V_ref(end-numVia+i,i)=1;
            V_ref(end-numVia+i,end-numVia*2+i)=-1;
            V_ref(end-numVia*2+i,end-numVia*2-size(Zc,1)+i)=1;
        end
        
        % change some entries in V_ref to Cd
        for i=1:shortNum
            V_ref(end-numVia+i,i)=kp*Cd;
            V_ref(end-numVia+i,end-numVia*2+i)=-kp*Cd;
        end
        A_BD2=MNA_BD2+V_ref+V_ref';
    end
    % begin solving AX=b
    unkX(:,:,itime) = A_BD2\(scaLC*(kpMins1*unkX(:,:,itime-1)+kpMins2*unkX(:,:,itime-2))+sparse(b(:,:,itime)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I vector, 3 dimension, I(:,:,t)
It = abs(unkX(numNodes+sourceNum+1:end-numVia,1,:));
% save results for plotting
save('./plotData/TDCurrent.mat','It', '-append');
%% Denormalize the unknown Vector into [Volts Amperes]
% current was in mA
unkX(numNodes+1:end-numVia*2,1,:) = unkX(numNodes+1:end-numVia*2,1,:)/1e3;
unkX(end-numVia+1:end,1,:)=unkX(end-numVia+1:end,1,:)/1e3;
toc

function [unkX,A]=oneFreqSolver(s,C_and_G,L,Rt,connMatA,connMatA_T,I_ref,...
    V_ref,b,numBranch,sourceNum,Zc,numNodes,sigma,tc)
t1=toc;
%--------------------------------------------------------------------------
% R matrix
% skin depth
mu0=pi*4e-7;
delta=sqrt(1/(abs(s)/2*1e9*mu0*sigma))*1e3; % in mm
if delta>tc
    t=tc;
else t=delta;
end
% divided by t from the previous step to get R matrix. Because of differential
% elements we have to  multiply by 2 in the end
R=2*Rt/t; 
%--------------------------------------------------------------------------
% assemble the MNA matrix without ESL's
pre_MNA=[s*C_and_G, connMatA_T;...
         connMatA,   s*L+R];
% pad additional zeros to MNA. Please refer to thesis equation (2.20)
% pad rows
MNA=[pre_MNA(1:numNodes,:) ;zeros(sourceNum,numNodes+numBranch); ...
    pre_MNA(numNodes+1:end,:)];
% pad comlums
MNA=[MNA(:,1:numNodes) zeros(numNodes+numBranch+1,sourceNum)...
    MNA(:,numNodes+1:end)];
A=[MNA,  I_ref;
    V_ref, -s*Zc];
t2=toc;
fprintf('Assemble MNA matrix : %2.2f sec \n', t2-t1);
% =========================================================================
%  Solve MNA matrix
% =========================================================================
disp('Solving MNA matrix...')
t1=toc;
unkX = A \ b ;
t2=toc;
fprintf('Solving time %f sec\n', t2-t1);
end
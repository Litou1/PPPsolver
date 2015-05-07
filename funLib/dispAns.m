% extract and display solution from unkX array
function [result,arg]=dispAns(SolverOptions,ESLflag,unkX,numVia,numNodes,...
    sourceNum)
V=unkX(1:numNodes+sourceNum);
% get (V4-Vp), thesis page 19
result=(V(numVia)-V(numNodes+1));
arg=0;
switch SolverOptions
    case 0
        if ESLflag==0
            result=1e3*result;
            fprintf('Lplane = %16.15f pH\n', full(result));
        else
            result=1e3*result;
            fprintf('Ltot = %16.15f pH\n', full(result));
        end
    case {1,2}
        arg=angle(full(result))*180/pi;
        result=abs(result);
        fprintf('Z = %16.15f ohms; phase=%2.3f \n', full(result),arg);
    case {3}
        fprintf('R = %16.15f ohms\n', full(result));
end
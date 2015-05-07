function t=solThknssCoef(d_to_uniform,subMeshLevel,meshT,f2level,f1Lmt,f2Lmt,...
    nthLmt)
% We want to solve the equation,
% d_to_uniform+meshT*2^n=a1*2^(n-1)+a2*2^(n-2)+....+a1*2+a0
% with all q>=meshT, meshT is the relative thickness for each level
% (measured by its current level)
%
% We know 1+2+...2^(n-1)=2^n-1,
% Multiply by meshT and subtract from both side,
% d_to_uniform+meshT=(a1-meshT)*2^(n-1)+(a2-meshT)*2^(n-2)+...
%     +(a1-meshT)*2+(a0-meshT)
% denote ai-meshT=qi for the following code
if subMeshLevel==1
    t=d_to_uniform+meshT*2;
    return
end
viaNum=length(d_to_uniform);
r=double(d_to_uniform+meshT);  % remained number of number of cells to be divided
% construct an array of each submesh's thickness
t=zeros(viaNum,subMeshLevel);
for n=1:viaNum  % for each via
    t0=meshT*ones(1,subMeshLevel);
    t1=zeros(1,subMeshLevel);% at least 1 for each level
    for i=1:subMeshLevel
        j=subMeshLevel-i;
        q=floor(r(n)/2^j);
        %   eg. 5 submesh levels
        %   r=q1*2^(5-1)+q2*2^(5-2)+q3*2^(5-3)+q4*2^(5-4)+q5*2^(5-5)
        %   r=q1*2^4    +q2*2^3    +q3*2^2    +q4*2^1    +q5*2^0
        % this if statement controls family I mesh levels
        % q>f1Lmt && (ismember(i,1:subMeshLevel-f2level))
        if q>f1Lmt && (i<=subMeshLevel-f2level)
            q=f1Lmt;
        end
        % this if statement controls family II mesh levels except for the 
        % smallest level
        % if q>f2Lmt && (ismember(i,subMeshLevel-f2level+1:subMeshLevel-1))
        if q>f2Lmt && (i>subMeshLevel-f2level) && (i<subMeshLevel)
            q=f2Lmt;
        end
        r(n)=r(n)-q*2^j;% update the remainder
        t1(i)=q;
    end
    t(n,:)=t1+t0;% this is the overall thickness
    % change the second smallest mesh level to the last level
    % iff the current value is smaller than the the user input
    if t(n,end)<=nthLmt
        % if there are more (n-1)th meshes than what we need
        % change at most floor((nthLmt-t(n,end))/2) meshes to the nth level
        numOfNminus1ToChange=floor((nthLmt-t(n,end))/2); % for n-1 level
        if t(n,end-1)-1>numOfNminus1ToChange
            t(n,end-1)=t(n,end-1)-numOfNminus1ToChange;
            t(n,end)=t(n,end)+2*numOfNminus1ToChange;
        else
            if t(n,end-1)>=2 % require at least 2 meshes for (n-1)th level
             % change all (addtional-1) meshes of (n-1)th level to the nth level
                t(n,end)=t(n,end)+2*(t(n,end-1)-2);
                t(n,end-1)=2;
            end
        end
        
    end
end
end
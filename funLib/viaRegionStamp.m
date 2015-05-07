function mat=viaRegionStamp(roundViaFLag,num,k,ctrNode,mat)
% This is a helper function for via nodes stamping
%
% assign -1 -2 -3... for short 1
%        -4 -5 -6... for short 2
% We will change the numbers to 1111 2222 ... by calling the method later on
% in Node class and Branch class

% assign -4 -5 -6... for source
% we will change the numbers to 3333 4444 ... by calling the method later on
% in Node class and Branch class

for i=1:num
    if k==4 && roundViaFLag==1
        % it's a little tricky to deal with rounded via
        % first generate viaIndMat as
        %      1     0     0     0     1
        %      0     0     0     0     0
        %      0     0     0     0     0
        %      0     0     0     0     0
        %      1     0     0     0     1
        viaIndMat=zeros(k+1,k+1);
        viaIndMat(1,1)=1;
        viaIndMat(1,k+1)=1;
        viaIndMat(k+1,1)=1;
        viaIndMat(k+1,k+1)=1;
        % negate it
        viaIndMatNeg=~viaIndMat;
        % assign vaia node numbers
        %      1    -4    -9   -14     1
        %     -1    -5   -10   -15   -19
        %     -2    -6   -11   -16   -20
        %     -3    -7   -12   -17   -21
        %      1    -8   -13   -18     1   
        viaIndMat(logical(viaIndMatNeg))=...
            (-i+1)*(k+1)^2-1:-1:(-i)*(k+1)^2+(k/2-1)*4;
        mat(ctrNode(i,1)-k/2:ctrNode(i,1)+k/2,...
            ctrNode(i,2)-k/2:ctrNode(i,2)+k/2)=viaIndMat;
        
    else
        mat(ctrNode(i,1)-k/2:ctrNode(i,1)+k/2,...
            ctrNode(i,2)-k/2:ctrNode(i,2)+k/2)=...
            reshape(((-i+1)*(k+1)^2-1:-1:(-i)*(k+1)^2)',k+1,[]);
    end
end

end
function [L, numXbranch, numYbranch]=L_mat(option,branchX,branchY,viaBranchX,...
    viaBranchY,h)
% This function generates L matrix, returns in uH!
%
% viaBranchX and viaBranchY are branches contained in a via
% Examples:
% [L, numXbranch, numYbranch]=L_mat(branchX,branchY,viaBranchX,viaBranchY,h)
%
%%
%branch count w/o deleting via branch
nx=numel(branchX.thisBranch);
ny=numel(branchY.thisBranch);
if option~=3
    Lx=L_matHelper(nx,branchX,h,'X');
    Ly=L_matHelper(ny,branchY,h,'Y');
else
    Lx=zeros(nx);
    Ly=zeros(ny);
end
% delete branches at vias
Lx (viaBranchX,:) = [];
Lx (:, viaBranchX') = [];
Ly (viaBranchY,:) = [];
Ly (:, viaBranchY') = [];
numXbranch=size(Lx,1);
numYbranch=size(Ly,1);
L = [             Lx,             sparse(numXbranch,numYbranch);...
    sparse(numYbranch,numXbranch),           Ly];

end

function [X,Y,Z]=J_TD(I,branchX,branchY,viaBranchX,viaBranchY,numXbranch,numYbranch,planeSizeX,planeSizeY)

%% current abs value and contour plot

Ix = ones(numel(branchX.thisBranch),1);
Ix(viaBranchX) = 0;
% Ix(Ix==1) = I(1:numXbranch);
Ix(Ix==1) = I(1:numXbranch,1);
 
Iy = ones(numel(branchY.thisBranch),1);
Iy(viaBranchY) = 0;
Iy(Iy==1) = I(numXbranch+1:numXbranch+numYbranch);

widthX = branchX.ey-branchX.sy;
widthY = branchY.ex-branchY.sx;

Jx=Ix./widthX;% surface current density mA/mm=A/m
Jy=Iy./widthY;

% interpolated points
INp=250;

Xx=branchX.center(:,1);
Xy=branchX.center(:,2);

Yx=branchY.center(:,1);
Yy=branchY.center(:,2);

xlin=linspace(0,planeSizeX,INp);
ylin=linspace(0,planeSizeY,INp);

[X,Y]=meshgrid(xlin,ylin);

U=griddata(Xx,Xy,Jx,X,Y,'cubic');
V=griddata(Yx,Yy,Jy,X,Y,'cubic');

Z=sqrt(U.^2+V.^2);


end

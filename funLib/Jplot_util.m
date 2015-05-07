% close all
load './plotData/nodesAndbranches.mat'
load './plotData/currentPlotting.mat'
% load './plotData/TDcurrent.mat'

%% current abs value and contour plot

Ix = ones(numel(branchX.thisBranch),1);
Ix(viaBranchX) = 0;
Ix(Ix==1) = I(1:numXbranch);
% Ix(Ix==1) = It(1:numXbranch,1,300);
 
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


hFig = figure;
set(hFig, 'Position', [500 500 500 500])
% colormap(flipud(gray))
colormap('default')
surf(X,Y,Z,'faceAlpha',.85);

shading interp
set(gca,'layer','bot')

% quiver(X,Y,U,zeros(size(V)));
% hold on 
% for i=1:numel(branchX.thisBranch)
%     x=[branchX.sx(i) branchX.ex(i) branchX.ex(i) branchX.sx(i) ];
%     y=[branchX.sy(i) branchX.sy(i) branchX.ey(i) branchX.ey(i) ];
%     plotc(x,y,'g');
% end

 
% quiver(X,Y,zeros(size(U)),V);
% hold on
% for j=1:numel(branchY.thisBranch)
%     x=[branchY.sx(j) branchY.ex(j) branchY.ex(j) branchY.sx(j) ];
%     y=[branchY.sy(j) branchY.sy(j) branchY.ey(j) branchY.ey(j) ];
%     plotc(x,y,'r');
% end


xlabel('mm','FontSize',20);
ylabel('mm','FontSize',20);


grid on
colorbar
% caxis([0 0.1]);
xlim([0 planeSizeX]);
ylim([0 planeSizeY]);
% zlim([0 1]);
view(-45,45)
 ch=colorbar;
     pause(1);
     delete(ch);
print('Jplot', '-dpng', '-r200')

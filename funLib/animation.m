%
% save time domain current density data to a gif
%
close all
load './plotData/nodesAndBranches.mat'
load './plotData/TDcurrent.mat'
OFFSET=5;
filename = './animation.gif';
for n = 2:20
    n
    [X,Y,Z]=J_TD(It(:,1,n),branchX,branchY,viaBranchX,viaBranchY,numXbranch,numYbranch,planeSizeX,planeSizeY);
    surf(X,Y,Z,'faceAlpha',.85);
    
    shading interp
    set(gca,'layer','bot')

    xlabel('mm','FontSize',20);
    ylabel('mm','FontSize',20);

    grid on
    colorbar
    % caxis([0 0.1]);
    zlim([0 100]);
    xlim([OFFSET planeSizeX-OFFSET]);
    ylim([OFFSET planeSizeY-OFFSET]);
    % zlim([0 1]);
    view(-45,45)
    ch=colorbar;
%     pause(1);
    delete(ch);
    
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n == 2;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end






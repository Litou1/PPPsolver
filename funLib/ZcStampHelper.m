function Zc= ZcStampHelper(PorG,viaOnes,viaTwos,oneNum,twoNum,viaSize,h,Zc)
% PorG -- Power or Ground via, string
for i=1:oneNum
    % fix a I via
    x1=viaOnes(i,1);
    y1=viaOnes(i,2);
    zs1=0;
    ze1=h;
    
    % Lp one-one, self
    for j=setdiff(1:oneNum,i)
        x2=viaOnes(j,1);
        y2=viaOnes(j,2);
        zs2=0;
        ze2=h;
        Lvia12Pos=filamentsLp12(zs1,ze1,x1,y1,zs2,ze2,x2,y2);
        % uncomment this if want to use 2-filament approximation
        % Lvia12Pos=fourFilamentsLp12(zs1,ze1,x1,y1,zs2,ze2,x2,y2,viaSize);
        if PorG=='G'
            Zc(i,j)=Lvia12Pos;
        else
            Zc(i+twoNum,j+twoNum)=Lvia12Pos;
        end
    end
    
    % Lp one-two, mutual
    for k=1:twoNum
        x2=viaTwos(k,1);
        y2=viaTwos(k,2);
        zs2=0;
        ze2=h;
        Lvia12Neg=filamentsLp12(zs1,ze1,x1,y1,zs2,ze2,x2,y2);
        %Lvia12Neg=fourFilamentsLp12(zs1,ze1,x1,y1,zs2,ze2,x2,y2,viaSize);
        if PorG=='G'
            Zc(i,k+oneNum)=Lvia12Neg;
        else
            Zc(i+twoNum,k)=Lvia12Neg;
        end
    end
    
end

end
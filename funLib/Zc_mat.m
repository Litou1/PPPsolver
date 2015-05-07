function Zc=Zc_mat(ESLflag,viaSize,h,short,source,shortNum,sourceNum,numVia)
if ESLflag==1;
    Zc=zeros(numVia,numVia);
    if h>viaSize
        Lvia11=tubewireLp11(viaSize,h);
    else
        Lvia11=roundwireLp11(viaSize,h);
    end
    %stamp the main diagonal terms
    for i=1:numVia
        Zc(i,i)=Lvia11;
    end
    %stamp the off diagonal terms
    Zc= ZcStampHelper('G',short,source,shortNum,sourceNum,viaSize,h,Zc);
    Zc= ZcStampHelper('P',source,short,sourceNum,shortNum,viaSize,h,Zc);
else
    Zc=zeros(numVia,numVia);
end
% normalized to uH
Zc=Zc*1e6;
end
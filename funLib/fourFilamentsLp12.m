% input units are in mm, H
% computes mutual-LpKm for 2 pairs of filaments using 2-filaments approximation
% for each pair
function LpKM=fourFilamentsLp12(zs1,ze1,x1,y1,zs2,ze2,x2,y2,d)
W1=0.2113248654;
W2=0.7886751346;
W=(W2-W1)/2;
r=d/2;
a=r*W;
two_x1=[x1-a x1+a];
two_y1=[y1-a  y1+a];
two_x2=[x2-a x2+a];
two_y2=[y2-a  y2+a];

[X1, Y1]=meshgrid(two_x1,two_y1);
[X2, Y2]=meshgrid(two_x2,two_y2);
X1=X1(:);Y1=Y1(:);
X2=X2(:);Y2=Y2(:);

LpKiM=zeros(4,1);
for i=1:4
    for j=1:4
        LpKiMj=filamentsLp12(zs1,ze1,X1(i),Y1(i),zs2,ze2,X2(j),Y2(j));
        LpKiM(i)=LpKiM(i)+LpKiMj;
    end    
end
LpKM=1/sum(1./LpKiM);
end
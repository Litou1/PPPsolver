% From Dr. Ruehli's PEEC book
% appendix D7 
% note that x y z definition is different between the book 
% and the program
% input units are in mm, H
function Lp12=filamentsLp12(xs1,xe1,y1,z1,xs2,xe2,y2,z2)
xs1=xs1/1e3;
xe1=xe1/1e3;
y1=y1/1e3;
z1=z1/1e3;
xs2=xs2/1e3;
xe2=xe2/1e3;
y2=y2/1e3;
z2=z2/1e3;

a1=xe2-xs1;
a2=xs2-xs1;
a3=xs2-xe1;
a4=xe2-xe1;
a=[a1 a2 a3 a4];
p=sqrt((y2-y1)^2+(z2-z1)^2);
sum=0;
for k=1:4
    sum=sum+(-1)^(k+1)*(a(k)*log(a(k)+sqrt(a(k)^2+p^2))-sqrt(a(k)^2+p^2));
end
 Lp12=sum*1e-7;
end
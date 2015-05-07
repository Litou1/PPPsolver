% From Dr. Ruehli's PEEC book
% appendix D11 
% tube conductor
% input units are in mm, H
% d-diameter
% l-length of the wire
function Lp11=roundwireLp11(d,l)
d=d/1e3;
l=l/1e3;
a=d/2;
Lp11=2e-7*l*(log(l/a+sqrt((l/a)^2+1))-sqrt(1+(a/l)^2)+a/l);
end
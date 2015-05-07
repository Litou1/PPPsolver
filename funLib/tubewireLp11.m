% From Dr. Ruehli's PEEC book
% appendix D11 
% tube conductor
% input units are in mm, H
% d-diameter
% l-length of the wire
function Lp11=tubewireLp11(d,l)
d=d/1e3;
l=l/1e3;
a=d/2;
k=2*a/l;
Lp11=pi*1e-7*l*((k^2/480+k^4/1280+1/3600)*pi^3 ...
                +(1/18-k^2/24)*pi...
                +(-2*log(l)+6*log(2)+2+2*log(a)-4*log(k*pi))/pi+8*a/l/pi^2);
end
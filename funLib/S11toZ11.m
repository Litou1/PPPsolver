%load S11 and convert to Z11

function [f Z]=S11toZ11(filename)
data=importdata(filename,' ',5);
Z0=50;
data=data.data;
f=data(:,1);
S11=data(:,2);
phase=data(:,3);
S=S11.*exp(1i*phase/180*pi);
Z = Z0 * (1 + S)./(1 - S);

end
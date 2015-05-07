function y=stepFuction(t,tr,A)
% time ns;rise time in ns; amplitute
tr=tr/1e9;% convert to sec
y=A*t/tr.*(t<=tr) + A*(t>tr);
end
function y=triPulse(t,tr,width,A)
% time vector [ns];rise time [ns];width [ns]; amplitute
numOfperoids=t(end)/(width*1e-9);
t=linspace(0,width*1e-9,numel(t)/numOfperoids);
tr=tr/1e9;% convert to sec
y=A*(t/tr.*(t<=tr) + (1-1*(t-tr)/tr).*(t>tr & t<=2*tr));
y=repmat(y,1,numOfperoids);
end
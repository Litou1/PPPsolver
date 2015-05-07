function Lskm = Lp2fMutApprox(delta_k,delta_m,h,r13,r14,r23,r24)
% 1 2 representing conductor k
% 3 4 representing conductor m
% difference inductance LpAB=1/4(Lp13+Lp14+Lp23+Lp24)
% Lp1fMutApprox calculates Lskm=2*(Lpkm-Lpkm')

Lsk1m1=Lp1fMutApprox(delta_k,delta_m,r13,h);
Lsk1m2=Lp1fMutApprox(delta_k,delta_m,r14,h);
Lsk2m1=Lp1fMutApprox(delta_k,delta_m,r23,h);
Lsk2m2=Lp1fMutApprox(delta_k,delta_m,r24,h);
% Lskm=2*(Lpkm-Lpkm')=2*(1/4*(...sum of all above))
Lskm=1/4*(Lsk1m1+Lsk1m2+Lsk2m1+Lsk2m2);
end
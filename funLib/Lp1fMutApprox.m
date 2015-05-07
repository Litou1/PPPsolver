%%% the approximation equation
%%% delta is the cell's length, unit is mm.
%%% r is the distance between two cells.
%%% return the inductance in uH
function Lp = Lp1fMutApprox(delta_k,delta_m,r,h)
Lp = zeros(length(delta_k),1);
Lp(:,1)=1e-4*delta_k.*delta_m.*h^2./(r.^3);
end














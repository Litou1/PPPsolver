function Lp = LpRectMutZeroT(xs1,xe1,ys1,ye1,z1,xs2,xe2,ys2,ye2,z2)
% usage:
%   Partial mutual inductance for two parallel current sheets in x direction
%   Lp = LpRectMutZeroT(xs1,xe1,ys1,ye1,xs2,xe2,ys2,ye2,z2)
%   
%    
%   X front, y right, z up
%   
%       
%       -------  xs1
%      / #1   / 
%     /      /          hight = z1
%    -------  xe1
%   ys1    ye1
%
%  Note that the same scheme is used for conductor #2 at hight = z2
%  The faces are parallel to x,y
%
% See Ruehli and Brennan, IEEE Trans. MTT-21, Feb. 1973, pp. 76-82
% Note corrections on z^2 from paper and (-1)^(k+m). The integrals for 
% potential coeff.  and partial inductance are same for rectangular shapes.
%
% Partial inductances are in uH and all geometrical parameters in mm!
ak(1,:) = (xs2 - xe1)*1;
ak(2,:) = (xs2 - xs1)*1;
ak(3,:) = (xe2 - xs1)*1;
ak(4,:) = (xe2 - xe1)*1;

bm(1,:) = (ys2 - ye1)*1;
bm(2,:) = (ys2 - ys1)*1;
bm(3,:) = (ye2 - ys1)*1;
bm(4,:) = (ye2 - ye1)*1;

SMALL_NO = 1e-37;
zDist = (z2 - z1)*1 + SMALL_NO;
zSq = zDist * zDist;


Lp = 0;
for k=1:4
  akSq = ak(k,:).* ak(k,:);
  coeffK = 0.5 * (akSq - zSq);  
  for m=1:4
    bmSq = bm(m,:).* bm(m,:);
    coeffM = 0.5 * (bmSq - zSq); 
    rho = sqrt(akSq + bmSq + zSq);
    
    prod1 = coeffM .* ak(k,:) .* log(ak(k,:) + rho + SMALL_NO);
    prod2 = coeffK .* bm(m,:) .* log(bm(m,:) + rho + SMALL_NO);
    prod3 = 1/6 * (bmSq - 2*zSq + akSq).* rho; 
    prod4 = bm(m,:)*zDist.*ak(k,:) .* atan(ak(k,:).*bm(m,:)./(rho*(zDist+SMALL_NO)));

    sign = 1 - rem(k+m,2)*2;
    Lp = Lp + sign * (prod1 + prod2 - prod3 - prod4);
  end
end
% Finally, for zero thickness inductance with current in x, 
% normalize by the x cross-sectional area which is just the width
% for zero thickness 
crossSection = (ye1-ys1)*1 .* (ye2-ys2)*1;
% save Lpp Lp
% save css crossSection
% note that mu_0 = 1e-4;   
% return
  Lp = 1e-4 * Lp./ (crossSection)';
end














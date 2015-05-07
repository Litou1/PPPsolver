function L=L_matHelper(n,branch,h,XorY)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% self terms in the main diagnonal Lpkk=2(Lpkk-Lpkk'), X is X branches and
% Y is Y branches
if XorY=='X'
    LSelf=abs(2*(LpRectMutZeroT(branch.sx,branch.ex,branch.sy,branch.ey,0,...
        branch.sx,branch.ex,branch.sy,branch.ey,0)...
        -LpRectMutZeroT(branch.sx,branch.ex,branch.sy,branch.ey,0,...
        branch.sx,branch.ex,branch.sy,branch.ey,h)));
elseif XorY=='Y'
    LSelf=abs(2*(LpRectMutZeroT(branch.sy,branch.ey,branch.sx,branch.ex,0,...
        branch.sy,branch.ey,branch.sx,branch.ex,0)...
        -LpRectMutZeroT(branch.sy,branch.ey,branch.sx,branch.ex,0,...
        branch.sy,branch.ey,branch.sx,branch.ex,h)));
end
L=sparse(branch.thisBranch,branch.thisBranch,LSelf,n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutual terms
sectionLen=branch.size;
ds1=zeros(n,n);% approximation criterion limit
ds2=zeros(n,n);% zero estimate criterion limit
% combination subscript indices
% t1=tic;
% comb2Idx=load('40200choose2.mat');
% t2=toc;
% loadTime=t2-t1;
comb2Idx=VChooseK(1:n,2);
% linear index
comb2linIdx=(comb2Idx(:,1)-1)*n+comb2Idx(:,2);
clear comb2Idx
% section 2-element length combination (n choose 2)
% here I use VChooseK.c by Jan Simon
% http://www.mathworks.com/matlabcentral/fileexchange/26190-vchoosek
sLenComb2=VChooseK(sectionLen,2);
normdist=max([sLenComb2 h*ones(n*(n-1)/2,1)],[],2);

% minimum section length less than 0.5mm indices
condIdx=min(sLenComb2,[],2)<.5;
clear sLenComb2

%%
% see thesis, page 33, table 3.1

ds1(comb2linIdx(condIdx))=normdist(condIdx)*5;
ds2(comb2linIdx(condIdx))=normdist(condIdx)*6;
% minimum section length larger than 0.5mm indices
ds1(comb2linIdx(~condIdx))=normdist(~condIdx)*3;
ds2(comb2linIdx(~condIdx))=normdist(~condIdx)*4;


% ds1(comb2linIdx(condIdx))=normdist(condIdx)*10;
% ds2(comb2linIdx(condIdx))=normdist(condIdx)*11;
% % minimum section length larger than 0.5mm indices
% ds1(comb2linIdx(~condIdx))=normdist(~condIdx)*8;
% ds2(comb2linIdx(~condIdx))=normdist(~condIdx)*9;
%%
clear comb2linIdx normdist condIdx

t1=toc;
distR=pdist(branch.center);% find center to center distances between point pairs
distMatR=trilform(distR);% convert it to matrix form
% closed-form indices
[ic,jc]=find((distMatR<ds1) & (distMatR~=0));
% approximation indices and its center to center distance r
[i_apprx,j_apprx]=find((ds1<=distMatR)&(distMatR<ds2));
r13_idx=(ds1<=distMatR)&(distMatR<ds2);
r14_idx=(ds1<=distMatR)&(distMatR<ds2);
r23_idx=(ds1<=distMatR)&(distMatR<ds2);
r24_idx=(ds1<=distMatR)&(distMatR<ds2);
clear distR distMatR

distR13=pdist(branch.node1);
distMatR13=trilform(distR13);
clear distR13
r13=distMatR13(r13_idx);
clear distMatR13 r13_idx

distR24=pdist(branch.node2);
distMatR24=trilform(distR24);
clear distR24
r24=distMatR24(r24_idx);
clear distMatR24 r24_idx

% pdist2 returns distances in "matrix form"!
distMatR14=pdist2(branch.node1,branch.node2);
r14=distMatR14(r14_idx);
clear distMatR14 r14_idx

distMatR23=pdist2(branch.node2,branch.node1);
r23=distMatR23(r23_idx);
clear distMatR23 r23_idx

t2=toc;
distTime=t2-t1;
% distMatR=tril(distMatR,-1); % only look at the lower triangular part
% distMatR13=tril(distMatR13,-1);
% distMatR14=tril(distMatR14,-1);

% clear distMatR distMatR13 distMatR14 distMatR23 distMatR24
%--------------------------------------------------------------------------
%Partial-mutual inductance calculation using closed-form formula
%
%        Lpkm=2(Lpkm-Lpkm')
%--------------------------------------------------------------------------
t1=toc;
if XorY=='X'
    LpMutClosed=abs(2*(LpRectMutZeroT(branch.sx(ic),branch.ex(ic),...
        branch.sy(ic),branch.ey(ic),0,...
        branch.sx(jc),branch.ex(jc),branch.sy(jc),branch.ey(jc),0)...
   -LpRectMutZeroT(branch.sx(ic),branch.ex(ic),branch.sy(ic),branch.ey(ic),0,...
        branch.sx(jc),branch.ex(jc),branch.sy(jc),branch.ey(jc),h)));
elseif XorY=='Y'
    LpMutClosed=abs(2*(LpRectMutZeroT(branch.sy(ic),branch.ey(ic),...
        branch.sx(ic),branch.ex(ic),0,...
        branch.sy(jc),branch.ey(jc),branch.sx(jc),branch.ex(jc),0)...
   -LpRectMutZeroT(branch.sy(ic),branch.ey(ic),branch.sx(ic),branch.ex(ic),0,...
        branch.sy(jc),branch.ey(jc),branch.sx(jc),branch.ex(jc),h)));
end
t2=toc;
zeroT_time=t2-t1;
L=L+sparse(ic,jc,LpMutClosed,n,n);
%--------------------------------------------------------------------------
%Partial-mutual inductance calculation using approximation
%-------------------------------------------------------------------------
t1=toc;
LpMutApprx = Lp2fMutApprox(branch.size(i_apprx),branch.size(j_apprx),...
    h,r13,r14,r23,r24);
t2=toc;
approx_time=t2-t1;
L=L+sparse(i_apprx,j_apprx,LpMutApprx,n,n);
end
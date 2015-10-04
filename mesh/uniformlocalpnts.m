function [plocal,tlocal]=uniformlocalpnts(porder)
%UNIFORMLOCALPNTS 2-D Mesh Generator for the master element.
%   [PLOCAL,TLOCAL]=UNIFORMLOCALPNTS(PORDER)
%
%      PLOCAL:    Node positions (NPL,2)
%      TLOCAL:    Triangle indices (NT,3)
%      PORDER:    Order of the complete polynomial 
%                 NPL = (PORDER+1)*(PORDER+2)/2
%
if porder==0
  plocal=[1/3 1/3 1/3];
  tlocal=[1 2 3];
else
  [u,v]=ndgrid((0:porder)/porder,(0:porder)/porder);
  plocal=[u(:),v(:)];
  plocal=[1-sum(plocal,2),plocal];
  plocal=plocal(plocal(:,1)>=0,:);
  if nargout>=2
    tlocal=zeros(0,3);
    loc=0;
    for ii=1:porder
      jj=porder+1-ii;
      t1=[1:jj; 2:jj+1; jj+1+(1:jj)]';
      t2=[2:jj; jj+1+(2:jj); jj+1+(1:jj-1)]';
      tlocal=[tlocal; t1+loc; t2+loc];
      loc=loc+jj+1;
    end
  end
end


function nfs=shape1d(porder,plocal,pts)
%SHAPE1D Calculates the nodal shapefunctions and its derivatives for 
%        the master 1D element [0,1]
%   NFS=SHAPE1D(PORDER,PLOCAL,PTS)
%
%      PORDER:    Polynomial order
%      PLOCAL:    Node positions (NP) (NP=PORDER+1)
%      PTS:       Coordinates of the points where the shape fucntions
%                 and derivatives are to be evaluated (npoints)
%      NSF:       shape function adn derivatives (np,2,npoints)
%                 nsf(:,1,:) shape functions 
%                 nsf(:,2,:) shape fucntions derivatives w.r.t. x
%

[nf,nfx]=koornwinder(pts,porder);
A=koornwinder(plocal(:,2),porder);
nfs=[nf/A,nfx/A];
nfs=reshape(nfs,[size(nf),2]);
nfs=permute(nfs,[2,3,1]);
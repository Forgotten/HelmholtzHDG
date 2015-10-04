function nfs=shape2d(porder,plocal,pts)
%SHAPE2D Calculates the nodal shapefunctions and its derivatives for 
%        the master triangle [0,0]-[1,0]-[0,1]
%   NFS=SHAPE2D(PORDER,PLOCAL,PTS)
%
%      PORDER:    Polynomial order
%      PLOCAL:    Node positions (NP,2) (NP=(PORDER+1)*(PORDER+2)/2)
%      PTS:       Coordinates of the points where the shape fucntions
%                 and derivatives are to be evaluated (npoints,2)
%      NFS:       shape function adn derivatives (np,3,npoints)
%                 nsf(:,1,:) shape functions 
%                 nsf(:,2,:) shape fucntions derivatives w.r.t. x
%                 nsf(:,3,:) shape fucntions derivatives w.r.t. y
%

[nf,nfx,nfy]=koornwinder(pts,porder);
A=koornwinder(plocal(:,2:3),porder);
nfs=[nf/A,nfx/A,nfy/A];

nfs=reshape(nfs,[size(nf),3]);
nfs=permute(nfs,[2,3,1]);
function [x,w]=gaussquad1d(pgauss)
%GAUSSQUAD1D Calculates the Gauss integration points in 1D for [0,1]
%   [X,W]=GAUSSQUAD1D(PGAUSS)
%
%      X:         Coordinates of the integration points 
%      W:         Weights  
%      PGAUSS:         Order of the polynomila integrated exactly
%
n=ceil((pgauss+1)/2);
P=jacobi(n,0,0);
x=sort(roots(P));

A=zeros(n,n);
for i=1:n
  P = jacobi(i-1,0,0);
  A(i,:)=polyval(P,x)';
end
w=A\[2;zeros(n-1,1)];

x=(x+1)/2;
w=w/2;

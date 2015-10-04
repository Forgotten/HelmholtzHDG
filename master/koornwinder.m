function [f,fx,fy]=koornwinder(x,porder)
%KOORNWINDER Vandermonde matrix for Koornwinder polynomials (works in 1D
%            and 2D). In 1D the Legendre polynomilas normalized to [0,1]
%            are used. In 2D the koornwinder basis is normalized to the 
%            (0,0),(1,0),(0,1) triangle. The basis are orthonormal.
%   [F,FX,FY]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)/2
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%
% - Written by: J. Peraire
%

dim=size(x,2);

switch dim
 case 1
     [f,fx]=koornwinder1d(x,porder);
 case 2
     [f,fx,fy]=koornwinder2d(x,porder);
 otherwise
     error('Strange case.')
end

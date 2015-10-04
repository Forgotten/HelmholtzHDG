function [f,fx,fy] = koornwinder2d(x,p)
%KOORNWINDER2D Vandermonde matrix for Koornwinder polynomials in 
%              the master triangle [0,0]-[1,0]-[0,1]
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
x = 2*x-1;
npol=prod((p+(1:2))./(1:2));
f=zeros(size(x,1),npol);
fx=zeros(size(x,1),npol);
fy=zeros(size(x,1),npol);

pq = pascalindex(npol);

xc = x;
xc(:,2) = min( 0.99999999, xc(:,2));

e(:,1) = 2*(1+xc(:,1))./(1-xc(:,2))-1;
e(:,2) = xc(:,2);
ii = find(x(:,2) == 1);
% Correct values for function evaluation
e(ii,1) = -1;
e(ii,2) =  1;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    
% Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1));

    f(:,ii) = fc*pval.*qval;
end

% Use displaced coordinate for derivative evaluation
e(:,1) = 2*(1+xc(:,1))./(1-xc(:,2))-1;
e(:,2) = xc(:,2);
de1(:,1) = 2./(1-xc(:,2));
de1(:,2) = 2*(1+xc(:,1))./(1-xc(:,2)).^2;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end

    dpp = polyder(pp);
    dqp = polyder(qp);

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));

    dpval = polyval(dpp,e(:,1));
    dqval = polyval(dqp,e(:,2));
    
% Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1));

    fx(:,ii) = fc*dpval.*qval.*de1(:,1);
    fy(:,ii) = fc*(dpval.*qval.*de1(:,2) + pval.*dqval);
end
fx = 2*fx;
fy = 2*fy;


function pq = pascalindex(p)
l=1;
for i=0:p
    for j=0:i
        pq(l,1)=i-j;
        pq(l,2)=j;
        l = l+1;
        if l>p 
           return;
        end
    end
end
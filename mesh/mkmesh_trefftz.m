function mesh = mkmesh_trefftz(m,n,porder,tparam)
%MKMESH_TREFFTZ Creates 2D mesh data structure for Trefftz airfoil.
%   MESH=MKMESH_TREFFTZ(M,N,PORDER,TPARAM)
%
%      MESH:      Mesh structure
%      M:         Number of points in the radial direction (default=15)
%      N:         Number of points in circumferential direction
%                 (default=30)
%      PORDER:    Polynomial Order of Approximation (default=3)
%      TPARAM:    Trefftz foil parameters
%                 TPARAM(1) = left x-shift of circle center 
%                             (trailing edge at (1,0)). (default=0.1)
%                 TPARAM(2) = y-shift of circle center. (default=0.05)
%                 TPARAM(3) = K-T exponent (=< 2) (2:Jukowski). (default=1.98)                      
%
%   See also: SQUAREMESH, FIXMESH (part of DISTMESH), MKT2F, SETBNDNBRS, 
%             UNIFORMLOCALPNTS, CREATENODES, TREFFTZ
%
% - Written by: J. Peraire
%
if nargin<2, m=15; n=30; end
if nargin<3, porder=3; end
if nargin<4, tparam=[0.1,0.05,1.98]; end

n = 2*ceil(n/2);
[p0,t0] = squaremesh(m,n/2,0);
[p1,t1] = squaremesh(m,n/2,1);
np = size(p0,1);
t1 = t1+np;
p1(:,2) = p1(:,2)+1.0;
mesh.p = [p0; p1];
mesh.t = [t0; t1];
%[mesh.p,mesh.t]=fixmesh(mesh.p,mesh.t);
clear p0; clear p1; clear t0; clear t1;

mesh.porder = porder;
[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.dgnodes = createnodes(mesh);

%First map to a cricle
mesh.p(:,1) = 2*mesh.p(:,1);
mesh.p(:,2) = pi*mesh.p(:,2);
z = mesh.p(:,1) + i*mesh.p(:,2);
w = exp(z);
mesh.p = [real(w),imag(w)];
[mesh.p,mesh.t] = fixmesh(mesh.p, mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);
mesh.fcurved = repmat(true,size(mesh.f,1),1);
mesh.tcurved = repmat(true,size(mesh.t,1),1);

bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

%Now let's try a K-T transformation
x0 = tparam(1);
y0 = tparam(2);
n  = tparam(3);

%First Rotate to ensure that a point stays at the trailing edge
rot = atan2(y0,1+x0);
r = sqrt((1+x0)^2 + y0^2);

w = mesh.p(:,1) + i*mesh.p(:,2);
w = r*exp(-i*rot)*w + complex(-x0,y0);

%Now K-T
z = ((w-1)./(w+1)).^n;
w = ((1+z)./(1-z))*n;
mesh.p = [real(w),imag(w)];

%Now the same for the dgnodes
z = 2*mesh.dgnodes(:,1,:) + i*pi*mesh.dgnodes(:,2,:);
w = exp(z);

w = r*exp(-i*rot)*w + complex(-x0,y0);

z = ((w-1)./(w+1)).^n;
w = ((1+z)./(1-z))*n;

mesh.dgnodes(:,1,:) = real(w);
mesh.dgnodes(:,2,:) = imag(w);


function mesh = mkmesh_square(m,n,porder,parity)
%MKMESH_SQUARE Creates 2D mesh data structure for unit square.
%   MESH=MKMESH_SQUARE(M,N,PORDER,PARITY)
%
%      MESH:      Mesh structure
%      M:         Number of points in the horizaontal direction 
%      N:         Number of points in the vertical direction
%      PORDER:    Polynomial Order of Approximation (default=1)
%      PARITY:    Flag determining the the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%
%     Mesh struct fields
%     NP:        Number of nodes
%     NT:        Number of triangles
%     P:         Node positions (NP,2)
%     T:         Triangle indices (NT,3)
%   See also: SQUAREMESH, MKT2F, SETBNDNBRS, UNIFORMLOCALPNTS, CREATENODES
%

if nargin<2, m=2; n=2; end
if nargin<3, porder=1; end
if nargin<4, parity=0; end

if m < 2 | n < 2,
    error('At least m=2,n=2 needed.');
end

[mesh.p,mesh.t] = squaremesh(m,n,parity);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>1-1e-3)', ...
           'all(p(:,2)>1-1e-3)','all(p(:,1)<1e-3)'};     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

mesh.fcurved = zeros(size(mesh.f,1),1);
mesh.tcurved = zeros(size(mesh.t,1),1);

mesh.porder = porder;
[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.dgnodes = createnodes(mesh);

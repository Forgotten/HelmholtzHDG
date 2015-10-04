function mesh = mkmesh_master(porder)
%MKMESH_MASTER Creates 2D mesh for a mesh consiting of one master element only.
%   MESH=MKMESH_MASTER(PORDER)
%
%      MESH:      Mesh structure
%      PORDER:    Polynomial Order of Approximation (default=1)
%
%   See also: UNIFORMLOCALPNTS, CREATENODES
%
% - Written by: J. Peraire
%
if nargin<1, porder=1; end;

mesh.p = [0,0;1,0;0,1];
mesh.t = [1,2,3];

[mesh.f,mesh.t2f] = mkt2f(mesh.t);
    
ii = find(mesh.f(:,4)==0);
mesh.f(ii,4) = -1;

mesh.fcurved = zeros(size(mesh.f,1),1);
mesh.tcurved = zeros(size(mesh.t,1),1);

mesh.porder = porder;
[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.dgnodes = createnodes(mesh);

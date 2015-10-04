function mesh = mkmesh_distort(mesh,wig)
%MKMESH_DISTORT Distorts a unit square mesh keeping boundaries unchanged
%   MESH = MKMESH_DISTORT(MESH,WIG)
%
%      MESH:     Mesh data structure
%                   input: mesh for the unit square created with
%                          mkmesh_square
%                   output: mesh for duct
%      WIG:      Amount of distortion (default: 0.05)
%
% - Written by: J. Peraire
%
if nargin < 2, wig = 0.05; end

dx =  -wig*sin(2*pi*(mesh.p(:,2)-0.5)).*cos(pi*(mesh.p(:,1)-0.5));
dy =  wig*sin(2*pi*(mesh.p(:,1)-0.5)).*cos(pi*(mesh.p(:,2)-0.5));
mesh.p = mesh.p + [dx,dy];
dx =  -wig*sin(2*pi*(mesh.dgnodes(:,2,:)-0.5)).*cos(pi*(mesh.dgnodes(:,1,:)-0.5));
dy =  wig*sin(2*pi*(mesh.dgnodes(:,1,:)-0.5)).*cos(pi*(mesh.dgnodes(:,2,:)-0.5));
mesh.dgnodes(:,1,:) = mesh.dgnodes(:,1,:) + dx;
mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) + dy;
mesh.fcurved = ones(size(mesh.fcurved));
mesh.tcurved = ones(size(mesh.tcurved));
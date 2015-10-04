function scaplot(mesh,u,clim,nref,pltmesh,surf)
%SCAPLOT  Plot Scalar function
%    SCAPLOT(MESH,U,CLIM,NREF,PLTMESH,SURF)
%
%    MESH:       Mesh structure
%    U(NPL,NT):  Scalar fucntion to be plotted
%                NPL = size(mesh.plocal,1)
%                NT = size(mesh.t,1)
%    CLIM:       CLIM(2) Limits for thresholding (default: no limits)
%    NREF:       Number of refinements used for plotting (default=0)
%    PLTMESH:    0 - do not plot mesh
%                1 - plot mesh with straight edges 
%                2 - plot mesh with curved edges (slow)
%    SURF:       0 - Normal 2D view
%                1 - 3D View
%

nt=size(mesh.t,1);
npl=size(mesh.plocal,1);

u=squeeze(u);

porder=mesh.porder;

if porder==0
  [plocal,tlocal]=uniformlocalpnts(1);
  dgnodes=zeros(3,2,nt);
  dgnodes(:,1,:)=reshape(mesh.p(mesh.t',1),3,1,nt);
  dgnodes(:,2,:)=reshape(mesh.p(mesh.t',2),3,1,nt);
  u=repmat(u,[3,1]);
else
  plocal=mesh.plocal;
  tlocal=mesh.tlocal;
  dgnodes=mesh.dgnodes;
end

if nargin<4 | isempty(nref), nref=ceil(log2(max(porder,1))); end

if nref>0
  A0=koornwinder(plocal(:,1:2),porder);
  [plocal,tlocal]=uniref(plocal,tlocal,nref);
  A=koornwinder(plocal(:,1:2),porder)/A0;
  npln=size(plocal,1);
  sz=size(dgnodes); if length(sz)==2, sz = [sz,1]; end
  dgnodes=reshape(A*reshape(dgnodes,npl,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
  u=A*u;
end

npln=size(plocal,1);
nodesvis=reshape(permute(dgnodes,[1,3,2]),[npln*nt,2]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

if nargin>=6 & ~isempty(surf)
   nodesvis = [nodesvis,reshape(u,size(nodesvis(:,1)))];
end

cla, patch('vertices',nodesvis,'faces',tvis,'cdata',u, ...
           'facecol','interp','edgec','none');
      
if nargin>=3 & ~isempty(clim)
  set(gca,'clim',clim);
end

if nargin>=5 & ~isempty(pltmesh) & pltmesh
  if pltmesh==2
    % Plot curved mesh
    e=boundedges(nodesvis,tvis);
    dgnodesltx=nodesvis(:,1);
    dgnodeslty=nodesvis(:,2);
    line(dgnodesltx(e'),dgnodeslty(e'),'color',[0,0,0]);
  else
    patch('vertices',mesh.p,'faces',mesh.t, ...
          'facecolor','none','edgecolor',[0,0,0]);
  end
end

set(gcf,'rend','z');
colorbar,axis equal,drawnow
if nargin>=6 & ~isempty(surf)
   cameramenu;
end


function e=boundedges(p,t)
%BOUNDEDGES Find boundary edges from triangular mesh
%   E=BOUNDEDGES(P,T)

% Form all edges, non-duplicates are boundary edges
edges=[t(:,[1,2]);
       t(:,[1,3]);
       t(:,[2,3])];
node3=[t(:,3);t(:,2);t(:,1)];
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:);
node3=node3(ix(qx));

% Orientation
v1=p(e(:,2),:)-p(e(:,1),:);
v2=p(node3,:)-p(e(:,1),:);
ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
e(ix,[1,2])=e(ix,[2,1]);

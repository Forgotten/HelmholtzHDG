function dgnodes=createnodes(mesh,fd,varargin)
%CREATEDGNODES Computes the Coordinates of the DG nodes.
%   DGNODES=CREATENODES(MESH,FD,FPARAMS)
%
%      MESH:      Mesh Data Structure
%      FD:        Distance Function d(x,y)
%      FPARAMS:   Additional parameters passed to FD
%      DGNODES:   Triangle indices (NPL,2,NT). The nodes on 
%                 the curved boundaries are projected to the
%                 true boundary using the distance function FD
%

if nargin < 2, fd=[]; end

p = mesh.p;
t = mesh.t;
plocal = mesh.plocal;

npl=size(plocal,1);
nt=size(t,1);

% Allocate nodes
dgnodes=zeros(npl,2,nt);
for dim=1:2
  for trinode=1:3
    dp=plocal(:,trinode)*p(t(:,trinode),dim)';
    dgnodes(:,dim,:)=dgnodes(:,dim,:)+permute(dp,[1,3,2]);
  end
end

% Project nodes on the curved boundary
if ~isempty(fd) 
  tc=find(mesh.tcurved);
  for it=tc'
    p = dgnodes(:,:,it);
    deps=sqrt(eps)*max(max(p)-min(p));
    ed = find(mesh.f(abs(mesh.t2f(it,:)),4)<0);
    for id=ed'
        e = find(mesh.plocal(:,id) < 1.e-6);
        d=feval(fd,p(e,:),varargin{:});
        dgradx=(feval(fd,[p(e,1)+deps,p(e,2)],varargin{:})-d)/deps;
        dgrady=(feval(fd,[p(e,1),p(e,2)+deps],varargin{:})-d)/deps;
        dgrad2=dgradx.^2+dgrady.^2;
        dgrad2(dgrad2==0)=1;
        p(e,:)=p(e,:)-[d.*dgradx./dgrad2,d.*dgrady./dgrad2];
    end
    dgnodes(:,:,it) = p;
  end
end
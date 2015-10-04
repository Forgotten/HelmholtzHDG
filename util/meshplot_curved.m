function hh=meshplot_curved(mesh,nref,opts)
%MESHPLOT  Plot Mesh Structure (with curved edges)
%   MESHPLOT(MESH,NREF,[OPTS])
%
%    MESH:       Mesh structure
%    NREF:       Number of refinements used for plotting (default=0)
%    OPTS:       (logical)
%      OPTS(1):  Plot dgnodes (default=0)
%      OPTS(2):  Plot triangle numbers (default=0)
%      OPTS(3):  Plot node numbers (default=0)


if nargin<2, nref=0; end
if nargin<3, opts=[0]; end
if length(opts)<2, opts=[opts,0]; end
if length(opts)<3, opts=[opts,0]; end

col=[.8,1,.8];

dgnodes=mesh.dgnodes;
tlocal=mesh.tlocal;
plocal=mesh.plocal;
porder=double(mesh.porder);

if nref>0
  A0=koornwinder(plocal(:,1:2),porder);
  [plocal,tlocal]=uniref(plocal,tlocal,nref);
  A=koornwinder(plocal(:,1:2),porder)/A0;
  dgnodes=reshape(A*reshape(dgnodes,size(A0,1),[]),size(A,1),2,[]);
end

e=boundedges(plocal(:,2:3),tlocal);
e1=segcollect(e);

clf,axis equal,axis off
nt=size(dgnodes,3);
hh=zeros(nt,1);
for it=1:nt
  px=dgnodes(:,1,it);
  py=dgnodes(:,2,it);
  hh(it)=patch(px(e1{1}'),py(e1{1}'),0.0*e1{1}', ...
               'facecolor',col,'edgecolor','k','Linew',1);
end

if opts(1)
  xx=squeeze(mesh.dgnodes(:,1,:));
  yy=squeeze(mesh.dgnodes(:,2,:));
  zz=0*xx;
  line(xx(:),yy(:),zz(:),'lines','n','marker','.','markersi',16,'col','b');
end

t = mesh.t;
p = mesh.p;

if opts(2)
    for it=1:size(t,1)
      pmid=mean(p(t(it,:),:),1);
      txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
      text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
end

if opts(3)
    for i=1:size(p,1)
      txtpars={'fontname','times','fontsize',20,'fontweight','bold', ...
               'horizontala','center','col','w', 'BackgroundColor',[0.5,0.5,0.5]};
      text(p(i,1),p(i,2),num2str(i),txtpars{:});
    end
end

if nargout<1, clear hh; end


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


function e1=segcollect(e)
%SEGCOLLECT Collect polygons from edge segments.

ue=unique(e(:));
he=histc(e(:),ue);
current=ue(min(find(he==1))); % Find an endpoint
if isempty(current) % Closed curve
  current=e(1,1);
end
e1=current;
while ~isempty(e)
  ix=min(find(e(:,1)==e1(end)));
  if isempty(ix)
    ix=min(find(e(:,2)==e1(end)));
    if isempty(ix) % >1 disjoint curves, recur
      rest=segcollect(e);
      e1={e1,rest{:}};
      return;
    end
    next=e(ix,1);
  else
    next=e(ix,2);
  end
  e1=[e1,next];
  e(ix,:)=[];
end
e1={e1};





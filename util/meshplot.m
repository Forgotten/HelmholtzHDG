function hh=meshplot(mesh,opts)
%MESHPLOT  Plot Mesh Structure (with straight edges)
%    MESHPLOT(MESH,[OPTS])
%
%    MESH:       Mesh structure
%    OPTS:       (logical)
%      OPTS(1):  Plot dgnodes (default=0)
%      OPTS(2):  Plot triangle numbers (default=0)
%      OPTS(3):  Plot node numbers (default=0)


if nargin<2 | isempty(opts), opts=[0]; end
if length(opts)<2, opts=[opts,0]; end
if length(opts)<3, opts=[opts,0]; end

p=mesh.p;
t=mesh.t;

hh=[];
pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1};
clf,hh=[hh;patch('faces',t,'vertices',p,pars{:})];
view(2),axis equal

if opts(1)
  xx=squeeze(mesh.dgnodes(:,1,:));
  yy=squeeze(mesh.dgnodes(:,2,:));
  zz=0*xx;
  line(xx(:),yy(:),zz(:),'lines','n','marker','.','markersi',16,'col','b');
end

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

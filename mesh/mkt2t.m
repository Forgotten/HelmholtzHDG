function t2t=mkt2t(t)
%MKT2T Compute Element to Element Connectivity.
%   T2T=MKT2T(T)
%
%      T:         Triangle indices (NT,3)
%      T2T:       Triangle to Trangle Connectivity (NT,3)
%                 T2T(IT,IN) is the trangle that shares an edge
%                 with triangle IT and does nont contain node T(IT,IN).
%                 When an element is adjacent to the boundary the
%                 corresponding entry in T2T is set to zero
%
nt=size(t,1);
edges=[t(:,[2,3])
       t(:,[3,1])
       t(:,[1,2])];
   
ts=[repmat(1:nt,1,3); kron(1:3,ones(1,nt))]';

edges=sort(edges,2);
[foo,foo,jx]=unique(edges,'rows');

[jx,ix]=sort(jx);
ts=ts(ix,:);
ix=find(diff(jx)==0);
ts1=ts(ix,:);
ts2=ts(ix+1,:);
t2t=zeros(nt,3);
t2t(ts1(:,1)+nt*(ts1(:,2)-1))=ts2(:,1);
t2t(ts2(:,1)+nt*(ts2(:,2)-1))=ts1(:,1);


function master=mkmaster(mesh,pgauss)
%MKMASTER  Initialize Master Element structures
%      MASTER = MKMASTER(MESH)
%
%      MESH:      Mesh data structure
%      PGAUSS:    Degree of the polynomil to be integrated exactly
%                 (default: PGAUSS = 4*MESH.PORDER)
%

if nargin < 2, 
    pgauss = 4*mesh.porder; 
end

% polinomial order of the approximation
master.porder = mesh.porder;
master.plocal = mesh.plocal;

ii=find(master.plocal(:,3)<1.e-6);
[xc,jj]=sort(master.plocal(ii,2));
master.ploc1d=master.plocal(ii(jj),1:2);

master.corner = [find(master.plocal(:,1) > 1-1.e-6);
                 find(master.plocal(:,2) > 1-1.e-6);
                 find(master.plocal(:,3) > 1-1.e-6)];

master.perm=zeros(master.porder+1,3);
for i=1:3
    ii=find(master.plocal(:,i)<1.e-6);
    [xc,jj]=sort(master.plocal(ii,mod(i+1,3)+1));
    master.perm(:,i)=ii(jj);
end
master.perm = [master.perm,flipud(master.perm)];
master.perm = reshape(master.perm,[master.porder+1,3,2]);

[master.gpts,master.gwgh] = gaussquad2d(pgauss);
[master.gp1d,master.gw1d] = gaussquad1d(pgauss);

master.shap=shape2d(master.porder,master.plocal,master.gpts);
master.sh1d=shape1d(master.porder,master.ploc1d,master.gp1d);

master.mass = squeeze(master.shap(:,1,:))*diag(master.gwgh)*squeeze(master.shap(:,1,:))';
master.conv(:,:,1) = squeeze(master.shap(:,1,:))*diag(master.gwgh)*squeeze(master.shap(:,2,:))';
master.conv(:,:,2) = squeeze(master.shap(:,1,:))*diag(master.gwgh)*squeeze(master.shap(:,3,:))';

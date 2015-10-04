function ustar = postprocess(master,mesh,master1,mesh1,u,q)
% POSTPROCESS postprocesses the HDG solution to obtain a better solution
%
%      ustar = postprocess(mesh,master,uh,qh)
%
%      MASTER:       Master structure of porder
%      MESH:         Mesh structure of porder
%      MASTER1:      Master structure of porder+1
%      MESH1:        Mesh structure of porder+1
%      U:            Approximate scalar variable
%      Q:            Approximate flux
%      USTAR:        Postprocessed scalar variable


ne    = size(mesh.dgnodes,3);
npl1  = size(mesh1.dgnodes,1);
ustar = zeros(npl1,1,ne);

shap    = squeeze(master.shap(:,1,:));
shap1   = squeeze(master1.shap(:,1,:));
shapxi1 = squeeze(master1.shap(:,2,:))*diag(master1.gwgh);
shapet1 = squeeze(master1.shap(:,3,:))*diag(master1.gwgh); 
shapxi2 = squeeze(master1.shap(:,2,:));
shapet2 = squeeze(master1.shap(:,3,:)); 

shapt = shape2d(master.porder,master.plocal,master1.gpts);
shap0 = squeeze(shapt(:,1,:));

for i=1:ne
    dg  = mesh.dgnodes(:,:,i);    
    dg1 = mesh1.dgnodes(:,:,i);    

    xxi = squeeze(master.shap(:,2,:))'*squeeze(dg(:,1));
    xet = squeeze(master.shap(:,3,:))'*squeeze(dg(:,1));
    yxi = squeeze(master.shap(:,2,:))'*squeeze(dg(:,2));
    yet = squeeze(master.shap(:,3,:))'*squeeze(dg(:,2));
    jac = xxi.*yet - xet.*yxi;
    
    xxi1 = squeeze(master1.shap(:,2,:))'*squeeze(dg1(:,1));
    xet1 = squeeze(master1.shap(:,3,:))'*squeeze(dg1(:,1));
    yxi1 = squeeze(master1.shap(:,2,:))'*squeeze(dg1(:,2));
    yet1 = squeeze(master1.shap(:,3,:))'*squeeze(dg1(:,2));
    jac1 = xxi1.*yet1 - xet1.*yxi1;
    shapx1 =   shapxi1*diag(yet1) - shapet1*diag(yxi1);
    shapy1 = - shapxi1*diag(xet1) + shapet1*diag(xxi1);
    shapx2 =   shapxi2*diag(yet1) - shapet2*diag(yxi1);
    shapy2 = - shapxi2*diag(xet1) + shapet2*diag(xxi1);

    F  = shap*(master.gwgh.*jac);
    F1 = shap1*(master1.gwgh.*jac1);
    K1 = (shapx2*diag(1./jac1)*shapx1'+shapy2*diag(1./jac1)*shapy1');
    Cx = shapx1*shap0';
    Cy = shapy1*shap0';

    K1(end,:) = F1';
    L = (Cx*q(:,1,i)+Cy*q(:,2,i));
    L(end,:) = F'*u(:,1,i);
    ustar(:,1,i) = K1\L;    
end



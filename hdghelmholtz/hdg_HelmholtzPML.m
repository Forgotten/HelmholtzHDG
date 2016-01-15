function [u,q,uhat,H,R] = hdg_HelmholtzPML(master, mesh,m, omega, tau, lambda, ...
                                        sigmaMax)
%HDG_Helmholtz Solves the Helmholtz problem using HDG
%
%   [u,q,uhat,H,R] = HDG_POI(MASTER,MESH,kappa, tau) 
%
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      m:            Squared slowness 1/c^2(x)
%      omega:        Frequency
%      tau:          Stabilization parameter
%      lambda:       Width of the pml (we suppose a squared computational
%      domain 
%      sigmaMax:     Maximum absorbtion 
%      we output H the HDG matrix to study it's sparsity pattern and the
%      right hand side R




nt  = size(mesh.t,1);
nf  = size(mesh.f,1);
% number of points in the volume
npv = size(mesh.dgnodes,1);
% number of points in the interface
nps = size(master.ploc1d,1);

dgnodes = mesh.dgnodes;
xmax = max(max(squeeze(dgnodes(:,1,:))));
xmin = min(min(squeeze(dgnodes(:,1,:))));
ymax = max(max(squeeze(dgnodes(:,2,:))));
ymin = min(min(squeeze(dgnodes(:,2,:))));

sigmaPML_x = @(x)sigmaMax*( (x-xmin-lambda).^2.*(x < xmin + lambda) + ...
                (x-(xmax-lambda)).^2.*(x > xmax - lambda))/lambda^2; 
sigmaPML_y = @(y) sigmaMax*( (y-ymin-lambda).^2.*(y < ymin + lambda) ...
                + (y-(ymax-lambda)).^2.*(y > ymax - lambda))/lambda^2;
s_x = @(x,y) (1+1i*sigmaPML_y(y)/omega)./(1+1i*sigmaPML_x(x)/omega);
s_y = @(x,y) (1+1i*sigmaPML_x(x)/omega)./(1+1i*sigmaPML_y(y)/omega);

s_xy = @(x,y) ((1+1i*sigmaPML_x(x)/omega).*(1+1i*sigmaPML_y(y)/omega));

% extracting information for the master node
% values at the dg nodes of the shape functions
shap = squeeze(master.shap(:,1,:));
% values of the derivatives with respect to xi
shapxi = squeeze(master.shap(:,2,:))*diag(master.gwgh);
% values of the derivatives with respect to eta
shapet = squeeze(master.shap(:,3,:))*diag(master.gwgh); 
sh1d = squeeze(master.sh1d(:,1,:));
perm = master.perm;

% compute the index mapping
elcon = elconnectivities(nps,mesh.t2f);

% allocate memory for the elemental matrices
A_K = zeros(2*npv,2*npv);
B_K = zeros(2*npv,npv);

% allocate memory for the global matrix and vector
H = sparse(nf*nps,nf*nps);
R = sparse(nf*nps,1);

% allocate memory 
P = zeros(3*npv,3*nps,nt);
L = zeros(3*npv,nt);


% TODO modify this thing in order to use sparse (cols, rows, values)
for i = 1:nt % loop over each element
    % computing the derivatives of (x,y) with respect to (xi, eta)
    % computing the Jacobian matrix for the change of variables
    xxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,1,i));
    xet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,1,i));
    yxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,2,i));
    yet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,2,i));
    % computing the determinant of Jacobian of the transformation
    jac = xxi.*yet - xet.*yxi;
    % derivative in x of the shape functions
    % we use tha chain rule
    % \partial_x \phi_i = \partial_xi  \phi \partial_x xi + 
    %                     \partial_eta \phi \partial_x eta 
    shapx =   shapxi*diag(yet) - shapet*diag(yxi);
    % derivatives in y of the shape functions using the same for y
    shapy = - shapxi*diag(xet) + shapet*diag(xxi);
    % the indices of this matrices are given by (I think)
    % shapx_{i,j} = \partial_x \phi_i(x_j)
    
      
    
    %computing the local Mass matrix
    Mass  = shap*diag(master.gwgh.*jac)*shap';
    
    % we extract the coordinated of the dg nodes
    xg = squeeze(master.shap(:,1,:))'*mesh.dgnodes(:,:,i);
    
    % evaluating S_x and S_y at the dgnodes
    s_x_dgnodes = s_x(xg(:,1),xg(:,2));
    s_y_dgnodes = s_y(xg(:,1),xg(:,2));
    
    Mass_x = shap*diag(master.gwgh.*jac.*(1./s_x_dgnodes))*shap';
    Mass_y = shap*diag(master.gwgh.*jac.*(1./s_y_dgnodes))*shap';
    
    % We start building the differential operators in Matrix form
    % form A_K 
    % \int_{E} Q \phi_i * \phi_j
    A_K(1:npv,1:npv) = Mass_x;
    A_K((npv+1):2*npv,(npv+1):2*npv) = Mass_y;
    
    % form B_K
    %\int_{E} \partial_x \phi_i * \phi_j
    B_K(1:npv,1:npv) = shapx*shap';
    %\int_{E} \partial_y \phi_i * \phi_j
    B_K((npv+1):2*npv,1:npv) = shapy*shap';    
    
    % this is the source term 
    % form F_K
    % evaluate the source at the Dg nodes
    fg = source(xg(:,1),xg(:,2));
    % projection into the DG space by integration
    % \int_{E} s(x,y) * \phi_j(x,y) dxdy
    F_K = shap*(master.gwgh.*jac.*fg);
    
    % we evaluate the slowness squared
    slowness_sqd = m(xg(:,1),xg(:,2));
    
    
    
    % allocate memory for the elemental matrices
    C_K = zeros(2*npv,3*nps);
    D_K = zeros(npv,npv);
    E_K = zeros(npv,3*nps);
    M_K = zeros(3*nps,3*nps);    
    
    for j = 1:3 % loop over each face of the element i
        % obtaining the indeces
        I = perm(:,j,1);
        J = ((j-1)*nps+1):j*nps;
        
        % sampling the derivatives of the shape function at the 
        % dg nodes within the boundary 
        xxi = squeeze(master.sh1d(:,2,:))'*squeeze(mesh.dgnodes(I,1,i)); 
        yxi = squeeze(master.sh1d(:,2,:))'*squeeze(mesh.dgnodes(I,2,i));  
        dsdxi = sqrt(xxi.^2+yxi.^2);
        % normal vector 
        nl = [yxi./dsdxi,-xxi./dsdxi];
        dws = master.gw1d.*dsdxi;                
                
        % form C_K
        % \int_{\partial T} \phi_i \nu_x \phi_j dS
        C_K(I,J) = C_K(I,J) + sh1d*diag(dws.*nl(:,1))*sh1d';
        % \int_{\partial T} \phi_i \nu_y \phi_j dS
        C_K(npv+I,J) = C_K(npv+I,J) + sh1d*diag(dws.*nl(:,2))*sh1d';

        % stabilization parameter 
        % \int_{\partial T} \phi_i \tau \phi_j dS
        tmp = sh1d*diag(tau*dws)*sh1d';
        
        % form D_K
        D_K(I,I) = D_K(I,I) + tmp;

        % form E_K
        E_K(I,J) = E_K(I,J) + tmp;

        % form M_K
        M_K(J,J) = M_K(J,J) + tmp;        
    end
    
    s_xy_dgnodes = s_xy(xg(:,1),xg(:,2));
    slowness_Mass = shap*diag(master.gwgh.*jac.*slowness_sqd.*s_xy_dgnodes)*shap';
    
    % \lbracket \tau u_h, w \rbracket - \omega^2 m (u_h, w)
    D_K = D_K - omega.^2*slowness_Mass; 
    
    % these are needed to compute u and q
    P(:,:,i) = [ A_K  B_K;...
                -B_K' D_K]\[C_K; E_K];
    L(:,i) = [ A_K  B_K; ...
              -B_K' D_K]\[zeros(2*npv,1); F_K];
   
    % form elemental stiffness matrix 
    H_K = M_K + [C_K' -E_K']*P(:,:,i);
            
    % form elemental residual vector 
    R_K = -[C_K' -E_K']*L(:,i); 
   
    % assemble the elemental matrices and vectors 
    % into the global matrix and vector
    imap = elcon(:,i);
    % extremely inefficient (we need to to do this with a proper list)
    H(imap,imap) = H(imap,imap) + H_K;
    R(imap) = R(imap) + R_K;
end

% Imposing Homogenous Dirichelt Boundary conditions
% remove rows and columns corresponding to the boundary nodes
m = nf - length(find(mesh.f(:,4)<0));
ind = m*nps+1:nf*nps;  
H(ind,:) = [];
H(:,ind) = [];
R(ind) = [];
        
% solve for uhat
uhat = H\R;
% we add the homegeneous Dirichlet boundary conditions
uhat = [uhat; zeros(length(ind),1)];

% compute u and q
u = zeros(npv,1,nt);
q = zeros(npv,2,nt);
for i = 1:nt
    imap = elcon(:,i);
    QU = L(:,i) + P(:,:,i)*uhat(imap);
    q(:,:,i) = reshape(QU(1:2*npv),[npv 2]);
    u(:,1,i) = QU(2*npv+1:end);    
end



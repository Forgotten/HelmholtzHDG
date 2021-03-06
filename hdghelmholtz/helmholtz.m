%DRIVER FOR THE POISSON PROBLEM 

porder = 4;
n      = 30;
omega  = 2*pi*10.2;
tau    = 1;

mesh = mkmesh_square(n,n,porder);
master = mkmaster(mesh,2*porder);

% define a constant wave speed
c = @(x,y) 1*x + 0*y +1;
m = @(x,y) 1./(c(x,y).^2);

[u,q,uhat] = hdg_Helmholtz(master, mesh, m, omega, 1/(porder*omega));

%% post processing (not valid for variable density)

% mesh1 = mkmesh_square(n,n,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% 
% ustar = postprocess(master,mesh,master1,mesh1,u,q);

figure(1); clf; scaplot(mesh,u); axis off;
% figure(2); clf; scaplot(mesh1,ustar); axis off;




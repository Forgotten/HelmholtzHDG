%DRIVER FOR THE POISSON PROBLEM 

porder = 4;
n      = 60;
omega  = 2*pi*30.2;
tau    = 1/(porder*omega);

mesh = mkmesh_square(n,n,porder);
master = mkmaster(mesh,2*porder);

% define a constant wave speed
c = @(x,y) .125*sin(8*pi*x) + .125*sin(12*pi*y) +1;
m = @(x,y) 1./(c(x,y).^2);

[u,q,uhat,H] = hdg_HelmholtzPML(master, mesh, m, omega, tau, 0.1, 500);


%% post processing (not valid for variable density)

% mesh1 = mkmesh_square(n,n,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% 
% ustar = postprocess(master,mesh,master1,mesh1,u,q);

figure(1); clf; scaplot(mesh,real(u)); axis off;
figure(2); clf; scaplot(mesh,imag(u)); axis off;
% figure(2); clf; scaplot(mesh1,ustar); axis off;




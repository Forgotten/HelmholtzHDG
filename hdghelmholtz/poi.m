%DRIVER FOR THE POISSON PROBLEM 

porder = 2;
n      = 30;
kappa  = 1;
tau    = 1;

mesh = mkmesh_square(n,n,porder);
master = mkmaster(mesh,2*porder);

[u,q,uhat] = hdg_helmholtz(master, mesh, kappa, tau);

mesh1 = mkmesh_square(n,n,porder+1);
master1 = mkmaster(mesh1,2*(porder+1));

ustar = postprocess(master,mesh,master1,mesh1,u,q);

figure(1); clf; scaplot(mesh,u); axis off;
figure(2); clf; scaplot(mesh1,ustar); axis off;




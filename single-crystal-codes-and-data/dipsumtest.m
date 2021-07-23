%dipole sum test, should calculate pseudocontact term in crystal frame

lat = [10.4974 0 0;0 6.1086 0;0 0 4.7358];

Lipos = lat * [0.5;0.5;0.5];
FePositions = [0.7180    0.2180    0.7820    0.2820; 0.7500    0.7500    0.2500    0.2500;...
    0.0248    0.4753    0.5248    0.9752];

Paramagnets = lat * FePositions;

%Dipole sum in a single unit cell

D = zeros(3,3,4);

muB = 9.274009994e-24;
T = 294.261; 
kB = 1.38064852e-23;
gam = 16.546;
S = 2;
gvalue = 2.12;
omega = -161;
mu0 = 4*pi*1e-7;
hbar = 1.054e-34;


for k = 1:length(Paramagnets)
    rdist = Paramagnets(:,k) - Lipos;
    D(:,:,k) = (3*(rdist*rdist') - ((norm(rdist)^2)*eye(3)))/(4*pi*norm(rdist)^5);
end

%7Li Pseudocontact tensor - sum over 1 unit cell

D
Dipoletensor = sum(D,3) %in angstroms^-3
pseudo = 1e6*mu0*muB^2*S*(S+1)/(3*kB*(T-omega))*Dipoletensor*(gvalue^2)*1e30


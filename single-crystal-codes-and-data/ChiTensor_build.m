function [pc1PAS, pc2PAS, D] = ChiTensor_build(lattice,positions)
%function [chi,molchi,vec,eigen,axial,rhombic] = ChiTensor(del,lattice,positions,density,molweight)

%computes the Magnetic susceptibility Tensor for the given material, and
%determines the maximum size of supercell which gives significant
%contribution to tensor
%To be used with 7Li shift tensor!

%DeL is the shift tensor LESS the contribution from the fermi contact

%lattice is (latfix) lattice vectors in angstrom units, lat = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0

%positions are atomic positions of paramagnetic sites in the unit cell (FePositionsfix)
%(rel), positions = 3 x n matrix, positions(:,n) = distance vector of nth site

%U(density) = g/cm^3 and U(molweight) = g/mol

tic;

alat = lattice(:,1);
blat = lattice(:,2);
clat = lattice(:,3);

n = 0;
D = zeros(3);
B = zeros(3);

%positions2 = [alat(1)*positions(1,:); blat(2)*positions(2,:); clat(3)*positions(3,:)];
%Lipos = [alat(1)*5.00012e-01; blat(2)*4.99969e-01; clat(3)*5.00003e-01]; %for FePO4 crystal
positions2 = lattice * positions;

Lipos = lattice * [0.5; 0.5; 0.5];

continuenow = true;
while continuenow 
    if mod(n,2) == 0
        fprintf('Checking sum in a %d x %d x %d shell...\n', n,n,n) 
    end
    for a = -n:n %dipole sum
        for b = -n:n
            for c = -n:n
                for z = 1:length(positions2)
                    ra = (positions2(:,z) + (a*alat) + (b*blat) + (c*clat)) - Lipos;
                    rb = ra/norm(ra);
                    D = D +(((3*(rb*rb')) - eye(3))/norm(ra)^3); 
                end
            end
        end
    end
    if(all(abs((D-B)) < abs(0.001*B))) %convergence condition
        continuenow = false; 
    else %starts new iteration
        B = D;
        D = zeros(3);
        n = n + 1;
    end
end
t = toc;

%chi = (4*pi)*(D\DeL); %(A^3) per particle susceptibility 
%molchi = (1e-6)*6.022e23 * chi; %(cm^3) per mol susceptibility
%molchi = molchi*(1e-24)*(density)/molweight; %unit-less PAS volume susceptibility
%[vec,eigeN] = eig(molchi); %~ component is diagonalizing matrix (use for g-tensor calc.)
%z = eigeN;


ge = 2.0023;
giso = 0;

muB = 9.274009994e-24;
S = input('Please input spin of TM site (2 for Fe, 2.5 for Mn): ');
T = 293;
k = 1.38064852e-23;
gam = input('Please input gyromagnetic ratio (17.235 for P31, 16.546 for Li7): ');
hbar = 1.054e-34;
mu0 = 4*pi*1e-7;
pc1 = mu0*muB^2*S*(S+1)/(3*k*(T))*(ge*D*1e36/4/pi);
pc2 = mu0*muB^2*S*(S+1)/(3*k*(T))*(giso*D*1e36/4/pi);

[~,pc1PAS] = eig(pc1);
[~,pc2PAS] = eig(pc2);






fprintf('The calculation converged after n = %d,  in %d seconds', n,t)
 

end





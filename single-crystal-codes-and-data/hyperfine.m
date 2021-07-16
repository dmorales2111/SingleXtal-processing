function [pc,hyp,D] = hyperfine(del,lattice,chi,phos,positions)
%takes in susceptibility tensor from ChiTensor, and lattice positions to 
%calculate hyperfine interaction for 31P shifts

%del is the shift tensor LESS the contribution from the fermi contact

%chi is the susceptibility tensor (in XTAL frame) that was calculated using
%ChiTensor

%chiPAS is the diagonalized susceptibility tensor
%vectors are the eigenvectors of chiPAS

%lat is the lattice vectors, lat = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0

%positions are atomic positions of paramagnetic sites in the unit cell
%(rel), positions = 3 x n matrix, positions(:,n) = distance vector of nth
%site

alat = lattice(:,1);
blat = lattice(:,2);
clat = lattice(:,3);

n = 0;
D = zeros(3);
B = zeros(3);

positions2 = lattice * positions;

Psite1 = lattice * phos;

continuenow = true;
while continuenow 
    if mod(n,2) == 0
        fprintf('Checking sum in a %d x %d x %d shell...\n', n,n,n) 
    end
    for a = -n:1:n %dipole sum
        for b = -n:1:n
            for c = -n:1:n
                for z = 1:length(positions2)
                    ra = (positions2(:,z) + (a*alat) + (b*blat) + (c*clat)) - Psite1;
                    rb = ra/norm(ra);
                    D = D +(((3*(rb*rb')) - eye(3))/(norm(ra)^3)); 
                end
            end
        end
    end
    if(all(abs((D-B)) < abs(0.1*B))) %convergence condition
        continuenow = false; 
    else %starts new iteration
        B = D;
        D = zeros(3);
        n = n + 1;
    end
end
t = toc;
fprintf('The calculation converged after n = %d,  in %d seconds\n', n,t)


pc = 1/(4*pi)*D*chi;
hyp = del - pc; %calculated fermi contact term

%Na = 6.02e23;
%muB = 9.274009994e-24;
%S = input('Please input spin of TM site (2 for Fe, 2.5 for Mn): ');
%T = 298;
%k = 1.38064852e-23;
%mu0 = 4*pi*10^-7;
%gam = 17.235;
%hbar = 1.054e-34;

%gtensorPAS = sqrtm(3*k*T*chiPAS/(S*(S+1)*Na*mu0*muB^2))
%gtensor = vectors*gtensorPAS*inv(vectors);

%A = inv(gtensor)*hyp*(3*hbar*k*T*gam/(S*(S+1)*muB)); %estimation of hyperfine coupling constant 


end
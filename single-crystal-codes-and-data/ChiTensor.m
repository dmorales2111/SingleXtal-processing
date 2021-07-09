function[chi,molchi,vec,eigen,axial,rhombic] = ChiTensor(del,lattice,positions,density,molweight)

%computes the Magnetic susceptibility Tensor for the given material, and
%determines the maximum size of supercell which gives significant
%contribution to tensor
%To be used with 7Li shift tensor!

%del is the shift tensor LESS the contribution from the fermi contact

%lat is the lattice vectors, lat = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0

%positions are atomic positions of paramagnetic sites in the unit cell
%(rel), positions = 3 x n matrix, positions(:,n) = distance vector of nth
%site

tic;


alat = lattice(:,1);
blat = lattice(:,2);
clat = lattice(:,3);

n = 0;
D = zeros(3);
B = zeros(3);
positions2 = lattice * positions;

Lipos = lattice * [0.5; 0.5; 0.5]; %same for all crystals

continuenow = true;
while continuenow 
    if mod(n,2) == 0
        fprintf('Checking sum in a %d x %d x %d shell...\n', n,n,n) 
    end
    for a = -n:1:n %dipole sum
        for b = -n:1:n
            for c = -n:1:n
                for z = 1:length(positions2)
                    ra = (positions2(:,z) + (a*alat) + (b*blat) + (c*clat)) - Lipos;
                    rb = ra/norm(ra);
                    D = D +(((3*(rb*rb')) - eye(3))/(norm(ra)^3)); 
                end
            end
        end
    end
    if(all(abs((D-B)) < abs(0.01*B))) %convergence condition
        continuenow = false; 
    else %starts new iteration
        B = D;
        D = zeros(3);
        n = n + 1;
    end
end
t = toc;


chi = 1e-6*4*pi*(D*del; %susceptibility per particle
molchi = 6.02e23 * chi; %molar susceptibility

[vec,eigen] = eig(molchi); %diagonalize
eigen = abs(eigen);

eigen = eigen * (density*1e6) / molweight;

rhombic = eigen(1,1) - eigen(2,2); %rhombic anisotropy, from Bertini et al.
axial = 1.5 * (eigen(3,3) - (trace(eigen)/3)); %axial anisotropy, from Bertini et al.
fprintf('The calculation converged after n = %d,  in %d seconds', n,t)
end





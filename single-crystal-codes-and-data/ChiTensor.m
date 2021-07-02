function[chi,n] = ChiTensor(del,lat,positions)

%computes the Magnetic susceptibility Tensor for the given material, and
%determines the maximum size of supercell which gives significant
%contribution to tensor
%To be used with 7Li shift tensor!
%sig is the shift tensor LESS the contribution from the fermi contact
%lat is the lattice vectors, lat = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0
%positions are atomic positions of paramagnetic sites in the unit cell
%(rel), positions = 3 x n matrix, positions(:,n) = distance vector of nth
%site

tic;


alat = lat(:,1);
blat = lat(:,2);
clat = lat(:,3);

n = 0;
D = zeros(3);
B = zeros(3);
positions2 = zeros(3,length(positions));
for a = 1:length(positions)
    positions2(:,a) = lat*positions(:,a);
end

continuenow = true;
while continuenow 
    if mod(n,2) == 0
        fprintf('Checking sum up to %d unit cells away...\n', n) 
    end
    for a = -n:1:n %dipole sum
        for b = -n:1:n
            for c = -n:1:n
                for z = 1:length(positions2)
                    ra = positions2(:,z) + (a*alat) + (b*blat) + (c*clat);
                    rb = ra/norm(ra);
                    D = D +((3*(rb*rb') - eye(3))/(norm(ra)^3)); 
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
chi = -4*pi*del*mldivide(D,eye(3));

fprintf('The calculation converged after n = %d,  in %d seconds', n,t)
end





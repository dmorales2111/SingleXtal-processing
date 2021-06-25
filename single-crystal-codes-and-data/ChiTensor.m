function[chi,n] = ChiTensor(sig,lat,positions)

%computes the Magnetic susceptibility Tensor for the given material, and
%determines the maximum size of supercell which gives significant
%contribution to tensor
%sig is the shift tensor LESS the contribution from the fermi contact, lat is the lattice vectors of the crystal, and
%positions are the atomic positions of paramagnetic sites in the unit cell

tic;
%muB = 9.274009994e-24;
%S = 0.5;
%T = 298;
%k = 1.38064852e-23;
%mu0 = 4*pi*10^-7;
%Na = 6.02e23;

%if nuc == '7Li'
%    gam = 16.546;
%elseif nuc =='31P'
%    gam = 17.235;
%end



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
        fprintf('Now testing %d by %d by %d supercell...\n', n,n,n) 
    end
    for a = -n:1:n %dipole sum
        for b = -n:1:n
            for c = -n:1:n
                for z = 1:length(positions2)
                    ra = positions2(:,z) + (a*alat) + (b*blat) + (c*clat);
                    rb = ra/norm(ra);
                    D = B +(eye(3) - 3*(rb*rb'))/(norm(ra)^3); 
                end
            end
        end
    end
    if(all(abs((D-B)) < abs(0.0001*B))) %convergence condition
        continuenow = false; 
    else %starts new iteration
        B = D;
        n = n + 1;
    end
end
t = toc;
chi = 4*pi*sig*mldivide(D,eye(3));
fprintf('The calculation converged after n = %d,  in %d seconds', n,t)
end





function [pc, hyp, g, A] = hyperfine(sig,lat,chi,positions)
%takes in susceptibility tensor from ChiTensor, and lattice positions to 
%calculate hyperfine interaction for 31P shifts

%sig is the shift tensor 

%lat is the lattice vectors, lat = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0

%positions are atomic positions of paramagnetic sites in the unit cell
%(rel), positions = 3 x n matrix, positions(:,n) = distance vector of nth
%site


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

Psite1 = [.41782;.25;.09477]*1e-10;
Psite2 = [.91789;.75;.40522]*1e-10;

continuenow = true;
while continuenow 
    if mod(n,2) == 0
        fprintf('Now testing %d by %d by %d supercell...\n', n,n,n) 
    end
    for a = -n:1:n %dipole sum
        for b = -n:1:n
            for c = -n:1:n
                for z = 1:length(positions2)
                    ra = (positions2(:,z) + (a*alat) + (b*blat) + (c*clat)) - Psite1;
                    rb = ra/norm(ra);
                    D = D +((eye(3) - 3*(rb*rb'))/(norm(ra)^3));  
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
fprintf('The calculation converged after n = %d,  in %d seconds', n,t)

pc = 1/(4*pi)*D*chi; %calcualted pseudocontact term



muB = 9.274009994e-24;
S = 0.5;
T = 298;
k = 1.38064852e-23;
mu0 = 4*pi*10^-7;
Na = 6.02e23;
gam = 17.235;

g = sqrtm(3*k*T*chi/(S*(S+1)*mu0*muB^2));

hyp = sig - pc; %calculated fermi contact term

A = hyp*3*k*T*gam/(S*(S+1)*muB)*g; %estimation of hyperfine coupling constant 


end
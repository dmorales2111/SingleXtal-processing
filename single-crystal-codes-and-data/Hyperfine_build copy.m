function [] = Hyperfine_build(LiDeL,PDeL,atom,lattice,positions,phos)

%Computes the individual paramagnetic contact and dipolar coupling terms,
%both in the crystal frame and PAS - valid for LiMnPO4 and LiFePO4 (mixed
%system will come later)

%LiDeL is the 7Li experimental shift tensor,PDeL is the 31P experimental shift
%tensor

%lattice is lattice vectors in angstrom units, lat = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0

%positions are fractional atomic positions of paramagnetic sites in the unit cell
%positions = 3 x n matrix, positions(:,n) = distance vector of nth site


%atom = "Fe" or "Mn"

Lipos = lattice * [0.5; 0.5; 0.5];
Ppos = lattice * phos;

LiD = dipsum(Lipos,lattice,positions);
PD = dipsum(Ppos,lattice,positions);


ge = 2.0023; %isotropic g value

%experimental g values for Fe and Mn respectively, found in Pigliapochi et al.
if atom == "Fe"
    giso = 2.12 - ge;
elseif atom == "Mn"
    giso = 2.00 - ge;
end

%constants
muB = 9.274009994e-24;
S = input('Please input spin of TM site (2 for Fe, 2.5 for Mn): ');
T = 293;
k = 1.38064852e-23;
gam = input('Please input gyromagnetic ratio (17.235 for P31, 16.546 for Li7): ');
hbar = 1.054e-34;
mu0 = 4*pi*1e-7;

%Dipolar terms

Lipc1 = 1e6*mu0*muB^2*S*(S+1)/(3*k*T)*(ge*LiD*1e30/4/pi); %Pigliapochi term (e) 
Lipc2 = 1e6*mu0*muB^2*S*(S+1)/(3*k*T)*(giso*LiD*1e30/4/pi); %Pigliapochi term (g)

Ppc1 = mu0*muB^2*S*(S+1)/(3*k*T)*(ge^2*PD*1e36/4/pi); %Pigliapochi term (e) 
Ppc2 = mu0*muB^2*S*(S+1)/(3*k*T)*(giso^2*PD*1e36/4/pi); %Pigliapochi term (g)

[~,~,Lidipiso1,Lidel1,Lieta1] = stats(Lipc1);
[~,~,Lidipiso2,Lidel2,Lieta2] = stats(Lipc2);

[~,~,Pdipiso1,Pdel1,Peta1] = stats(Lipc1);
[~,~,Pdipiso2,Pdel2,Peta2] = stats(Lipc2);
 

%Contact terms and hyperfine coupling constant calculations

%whole contact term, Pigliapochi terms (a), (c), (d)
Lihyp = LiDeL - Lipc1 - Lipc2;
Phyp =  PDeL - Ppc1 - Ppc2;


save calculated_params.mat LiDeL PDeL Lipc1 Lipc2 Ppc1 Ppc2 Lidipiso1 ...
    Lidel1 Lieta1 Lidipiso2 Lidel2 Lieta2 Pdipiso1 Pdel1 Peta1 Pdipiso2 ...
    Pdel2 Peta2 Lihyp Phyp

end

function D = dipsum(nucsite,lattice,positions)
tic;

alat = lattice(:,1);
blat = lattice(:,2);
clat = lattice(:,3);

metalpositions = lattice * positions;

n = 0;
D = zeros(3); %init D tensor
B = zeros(3); %init convergence check tensor
continuenow = true;
while continuenow 
    if mod(n,2) == 0
        fprintf('Checking sum in a %d x %d x %d shell...\n', n,n,n) 
    end
    for a = -n:n %dipole sum
        for b = -n:n
            for c = -n:n
                for z = 1:length(metalpositions)
                    ra = (metalpositions(:,z) + (a*alat) + (b*blat) + (c*clat)) - nucsite;
                    rb = ra/norm(ra);
                    D = D +(((3*(rb*rb')) - eye(3))/norm(ra)^3); 
                end
            end
        end
    end
    if(all(abs((D-B)) < abs(0.05*B))) %convergence condition
        continuenow = false; 
    else %starts new iteration
        B = D;
        D = zeros(3);
        n = n + 1;
    end
end
t = toc;
fprintf('The calculation converged after n = %d,  in %d seconds\n', n,t)
end

function [v,d,iso,del,eta] = stats(matrix)
%takes a rank 2 tensor, diagonalizes it, then calculates by herlerben
%convention:

%|d33 - iso| > |d11 - iso| > |d22 - iso|
%isotropic value iso = (d11 + d22 + d33)/3
%symmetric anisotropy del = (iso - d33)
%asymmetry parameter eta = (d11 - d22)/(iso - d33)

[v,d] = eig(matrix);
iso = trace(d)/3;

dif = zeros(3,2);

for k = 1:length(matrix)
    dif(k,1) = d(k,k);
    dif(k,2) = abs(d(k,k) - iso);
end

dif = sortrows(dif,2,'descend'); prin = dif(:,1);

del = prin(1) - iso;
eta = (prin(3) - prin(2))/del;
end

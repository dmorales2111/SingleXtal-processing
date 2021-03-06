function [] = Hyperfine_build(LiDeL,PDeL,atom,lattice,positions,phos)

%Chemical Shift Tensor separation program, written by D.Morales 2021

%Computes the individual paramagnetic contact (Fermi Contact) and 
%dipolar coupling (pseudocontact) terms,
%both in the crystal frame and PAS - for LiMnPO4 and LiFePO4

%Dipolar term is first computed in crystal frame, then subtracted from
%experimental shift tensors for Contact component

%LiDeL is the 7Li experimental shift tensor, PDeL is the 31P experimental shift
%tensor (crystal frame)

%atom = "Fe" or "Mn", depending on which crystal being computed

%lattice is lattice vectors in angstrom units, lattice = 3 x 3 matrix, diagonals are a,b,c
%lattice constants, off diagonals are 0

%positions are fractional atomic positions of paramagnetic sites in the unit cell
%positions = 3 x n matrix, positions(:,n) = distance vector of nth site

%phos are the fractional coordinates of the phosphorus site for 31P
%calculations



Lipos = lattice * [0.5; 0.5; 0.5]; %center Li site
Ppos = lattice * phos;

ge = 2.00231930436256; %free electron g value

%experimental g values and weiss constants for Fe and Mn respectively, found in Pigliapochi et al.
if atom == "Fe"
    giso = 2.12 - ge;
    theta = -161;
    S = 2; %4 unpaired d electrons
elseif atom == "Mn"
    giso = 2.00 - ge;
    theta = -69; %averaged values referenced from doi.org/10.1039/c1dt10515a, doi.org/10.1021/ja017838m, and doi.org/10.1016/j.jcrysgro.2009.01.090
    S =  2.5; %5 unpaired d electrons
end

%constants, in SI units
muB = 9.274009994e-24;
T = 294.261; %Lab temperature, 70 deg F
k = 1.38064852e-23;
gam = [17.235; 16.546]; %gyromagnetic ratios of 31P and 7Li in MHz/T
mu0 = 4*pi*1e-7;


%chemical shift tensor PAS
[Lvec,Lprin] = eig(LiDeL);
[Pvec,Pprin] = eig(PDeL);

%%
%Dipolar terms

LiD = dipsum(Lipos,lattice,positions);
PD = dipsum(Ppos,lattice,positions);


LiPC = 1e6*mu0*muB^2*S*(S+1)/(3*k*(T - theta))*((ge+giso)^2)*(LiD*1e30);
PPC = 1e6*mu0*muB^2*S*(S+1)/(3*k*(T - theta))*((ge+giso)^2)*(PD*1e30);

[~,~,~,Lidel1,Lieta1] = stats(LiPC);
[~,~,~,Pdel1,Peta1] = stats(PPC);

 
%%
%Contact terms and hyperfine coupling constant calculations

%Isotropic Fermi contact term
Lihyp = LiDeL - LiPC;
Phyp = PDeL - PPC;


[~,~,Liso,Liisodel,Liisoeta,Liisoprin] = stats(Lihyp);
[~,~,Piso,Pisodel,Pisoeta,Pisoprin] = stats(Phyp);

%Anisotropic Fermi contact term
Lianiso = Lihyp - Liso*eye(3);
Paniso = Phyp - Piso*eye(3);

[~,~,~,Lidel2,Lieta2] = stats(Lianiso);
[~,~,~,Pdel2,Peta2] = stats(Paniso);

%Hyperfine coupling constant estimation
LiA = Liso*1e-6*3*k*(T - theta)*gam(2)/(muB*S*(S+1))/(ge+giso)*2*pi;
PA = Piso*1e-6*3*k*(T - theta)*gam(1)/(muB*S*(S+1))/(ge+giso)*2*pi;

save calculated_params.mat LiDeL PDeL LiD PD  ...
    Lidel1 Lieta1 Lidel2 Lieta2 Pdel1 Peta1 ...
    Pdel2 Peta2 Lihyp Phyp Lprin Pprin Liso Liisodel Liisoeta Liisoprin LiA ...
    Piso Pisodel Pisoeta Pisoprin PA LiPC PPC Lvec Pvec

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
                    D = D +(1/(4*pi)*(((3*(rb*rb')) - eye(3))/norm(ra)^3)); 
                end
            end
        end
    end
    if(all(abs((D-B)) < .0005)) %convergence condition, makes sense if n < 3 due to 1/r^3 dependence
        continuenow = false; 
    else %starts new iteration
        B = D;
        D = zeros(3);
        n = n + 1;
    end
end
t = toc;
fprintf('The calculation converged after n = %d,  in %d seconds\n\n', n,t)
end

function [v,d,iso,del,eta,prin] = stats(matrix)
%takes a rank 2 tensor, diagonalizes it, then calculates by herlerben
%convention:

%|d33 - iso| > |d11 - iso| > |d22 - iso|
%isotropic value iso = (d11 + d22 + d33)/3
%anisotropy del = d33 - iso
%asymmetry parameter eta = (d22 - d11)/(iso - d33)

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

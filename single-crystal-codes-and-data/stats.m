function [d,iso,del,eta] = stats(matrix)
%takes a rank 2 tensor, diagonalizes it, then calculates
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


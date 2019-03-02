
function res = QAOA(p, n)
%We first load the Hamiltonian
H = {};
for i = 1:n-1
    H = {H{:},{{i, i+1}, h}}
end



end
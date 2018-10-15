%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% Compute lagrange for T-MMSE scheme. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v_real = cmp_lagrange(A, P, Es)

v = [];

[U, Lambda] = eig(A);
M = rank(Lambda);
lambda = real(diag(Lambda));
[tmp, ind] = sort(lambda, 'descend');
lambda = lambda(ind(1 : M));


exp1 = [1, 2 * lambda(1), lambda(1)^2];
for j = 2 : M
    exp1 = [exp1, 0, 0] + [0, exp1 * 2 * lambda(j), 0] + [0,0, lambda(j)^2 * exp1];
end
exp1 = P / Es * exp1;

exp3 = 0;
for j = 1 : M
    exp2 = 1;
    for k = 1 : M
        if(k ~= j)
            exp2 = [exp2, 0, 0] + [0, exp2 * 2 * lambda(k), 0] + [0,0, lambda(k)^2 * exp2];
%             exp2 = expand((exp2) * (expand(lambda(k) + x)^2));%
        end
    end
    exp2 = lambda(j) * exp2;%
    exp3 = exp3 + exp2;
end

exp = [exp1(1:2), exp1(3:end) - exp3];

v = roots(exp);

v_real = [];
for i = 1 : length(v)
    if(imag(v(i)) == 0)
        v_real = [v_real, v(i)];
    end
end
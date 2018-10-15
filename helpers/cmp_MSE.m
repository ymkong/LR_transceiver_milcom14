%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% Compute MSE for T-MMSE scheme. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mse_c = cmp_MSE(G, F, H, N_r, N_t, B, users, sigma, Es)

mse_c = 0;

for user_id = 1 : users
    G_i = G(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id))); % receive matrix of 'user_id'
    H_i = H(sum(N_r(1 : (user_id - 1))) + 1 : sum(N_r(1 : user_id)), :);
    F_i = F(:, sum(B(1 : (user_id - 1))) + 1 : sum(B(1 : user_id)));
    mse_c = mse_c + trace(G_i * H_i * F * F' * H_i' * G_i' + sigma^2 / Es * G_i * G_i' - F_i' * H_i' * G_i' - G_i * H_i * F_i + eye(B(user_id)));
end
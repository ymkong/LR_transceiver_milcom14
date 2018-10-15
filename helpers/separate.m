%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% Get separate channels for individual users for S-MMSE-max-SNR scheme. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_separate = separate(H, N_r)

users = size(N_r, 2);
H_separate = cell(1, users);

for user_id = 1 : users

    H_separate{user_id} = H(sum(N_r(1 : (user_id-1))) + 1 : sum(N_r(1 : user_id)), :);

end
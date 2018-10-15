%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
%
%
% calculate total precoding matrix F = [F1, F2, ..., FK]; 
% F: M_T * r
% Fi: M_T * ri
%
% initialize receive matrices
% loop:-------
% calcualte langrange multiplier (need some work)
% compute transmit matrices
% update receive matrices
% ------------
% 
%
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, G] = t_mmse(num_iteration, H, P, N_r, N_t, B, users, sigma, Es)

M_T = sum(N_t);
M_R = sum(N_r);

% initialize F and G
F = zeros(M_T, sum(B));
G = zeros(sum(B), M_R);

% store individual channel matrix
H_i = cell(1, users);
for user_id = 1 : users
    H_i{user_id} = H(sum(N_r(1 : (user_id - 1))) + 1 : sum(N_r(1 : user_id)), :); % the channel matrix for 'user_id'
end

% initialize receive matrices
for user_id = 1 : users
    G(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id))) ...
        = eye_rect(B(user_id), N_r(user_id));
end

G_last = G;

for iter = 1 : num_iteration
    % compute A
    A = zeros(M_T, M_T);
    for user_id = 1 : users
%         H_i = H(sum(N_r(1 : (user_id - 1))) + 1 : sum(N_r(1 : user_id)), :); % the channel matrix for 'user_id'
        G_i = G(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id))); % receive matrix of 'user_id'
        A = A + H_i{user_id}' * G_i' * G_i * H_i{user_id};
    end
    
    % compute lagrange multiplier
    v = cmp_lagrange(A, P, Es);
    
    % find the v that gives min Total MSE
    if(isempty(v))
        disp('error in computing lagrange multiplier!!!');
    else
        MinMSE = 1000000;
        minF = zeros(M_T, sum(B));
        minG = zeros(sum(B), M_R);
        
        for i = 1 : length(v)
            
            G_tmp = zeros(sum(B), M_R);
            
            % compute transmit matrices
            for user_id = 1 : users
                G_i = G(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id))); % receive matrix of 'user_id'
                F(:, sum(B(1 : (user_id - 1))) + 1 : sum(B(1 : user_id))) = inv(A + v(i) * eye(size(A, 2))) * H_i{user_id}' * G_i';
            end
    
            % update receive matrices
            for user_id = 1 : users
                F_i = F(:, sum(B(1 : (user_id - 1))) + 1 : sum(B(1 : user_id)));
                G_tmp(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id))) = F_i' * H_i{user_id}' * inv(H_i{user_id} * F * F' * H_i{user_id}' + sigma ^2 / Es * eye(N_r(user_id)));
            end
            
            currentMSE = cmp_MSE(G_tmp, F, H, N_r, N_t, B, users, sigma, Es);
            if(abs(currentMSE) < MinMSE)
                minF = F;
                minG = G_tmp;
                MinMSE = abs(currentMSE);
            end
        end
        
        F = minF;
        G = minG;
        
    end
    if(sum(sum(abs(G - G_last) .^ 2)) < 0.0001)
        break;
    else
        G_last = G;
    end
    
end

end
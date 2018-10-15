%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
%
% successively calculating the columns of the precoding matrix F for each of the receive antennas separately.
% columns in the precoding matrix Fi, each corresponding to one receive antenna, are calculated successively.
% max SNR design.
% See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mdl, m ] = det_s_mmse_maxsnr(act, mdl, m, mu_name, lr_func, lr_name)

if (strcmp(act, 'updateH'))
    if(~isfield(mdl.chn_info, mu_name))    
        N_r = mdl.N_r;
        N_t = mdl.N_t;
        B = mdl.B;
        NS = length(mdl.SNRdb);
                
        P = mdl.P;
        H = m.H;
        Fs = cell(1, NS);
        Gs = zeros(sum(B), sum(N_r), NS);
        users = mdl.users;
        
        Alphas = zeros(1, NS);

        for SNR_ind = 1 : NS    
            sigma = mdl.sigmas(SNR_ind);
            for user_id = 1 : users
                F_i = []; % precoding matrix for each user
                H_i_hat = [H(1 : sum(N_r(1: (user_id-1))), :); H(sum(N_r(1:user_id))+1 : end, :)];
                H_i = H(sum(N_r(1:(user_id-1)))+1 : sum(N_r(1:user_id)), :);
                
                for j = 1 : N_r(user_id)    
                    H_ij_hat = [H_i(j, :); H_i_hat];
                    
                    alpha = sigma .^2 * sum(N_r) / P;
                    
                    F_ij = inv(H_ij_hat' * H_ij_hat + alpha * eye(sum(N_t))) * H_ij_hat';
                                        
                    col = F_ij(:, 1);
                    
                    F_i = [F_i, col];
                end

                H_e = H_i * F_i;
                
                [UU,SS,VV] = svd(H_e);
                
                F = F_i * VV(:, 1);
                
                Fs{SNR_ind} = [Fs{SNR_ind}, F];

                Gs(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id)) , SNR_ind) = UU(:, 1)' ./ SS(1 , 1);
                
            end
            beta = sqrt(P / trace(Fs{SNR_ind} * Fs{SNR_ind}') / mdl.Es);
            
            Fs{SNR_ind} = beta * Fs{SNR_ind};
            
            Alphas(1, SNR_ind) = beta;
            
        end 
        
        mdl.chn_info.(mu_name).Gs = Gs;
%         mdl.chn_info.(mu_name).Ts = Ts;
%         mdl.chn_info.(mu_name).ITs = ITs;
        
        m.Fs = Fs;
        m.alphas = Alphas;
        
    end
elseif (strcmp(act, 'det'))
    
    B = mdl.B;
    G = mdl.chn_info.(mu_name).Gs(:, :, m.SNR_ind);
    m.s_hat = G * m.y / m.alphas(1, m.SNR_ind);
    
end

end
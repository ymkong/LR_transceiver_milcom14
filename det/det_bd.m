%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% precoding with Block diagonalization (BD). See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mdl, m ] = det_bd(act, mdl, m, mu_name)

if (strcmp(act, 'updateH'))
    if(~isfield(mdl.chn_info, mu_name))    
        N_r = mdl.N_r;
        N_t = mdl.N_t;
        B = mdl.B;
        NS = length(mdl.SNRdb);
        Gs = zeros(sum(B), sum(N_r), NS);
        P = mdl.P;
        H = m.H;
        Fs = cell(1, NS);
        users = mdl.users;

        for user_id = 1 : users

            H_i_hat = [H(1 : sum(N_r(1: (user_id-1))), :);H(sum(N_r(1:user_id))+1 : end, :)];
            H_i = H(sum(N_r(1:(user_id-1)))+1 : sum(N_r(1:user_id)), :);

            [UU, SS, VV] = svd(H_i_hat);

            % take the null space of H_i_hat
            VV = VV(:, rank(SS) + 1 : end);

            % Get A, and make sure VV * A meets the power constraints  
            % Es * tr(FF') <= P / users
            A = [eye(B(user_id)); zeros(size(VV, 2) - B(user_id), B(user_id))];
            beta = P ./ sqrt(mdl.users) ./ sqrt(B(user_id)) ./ sqrt(mdl.Es);
            F = beta .* VV * A;

            for SNR_ind = 1 : NS
                Fs{SNR_ind} = [Fs{SNR_ind}, F];
                sigma = mdl.sigmas(SNR_ind);
                C = inv((H_i * F)' * (H_i * F) + sigma .^ 2 * eye(B(user_id)) ./ mdl.Es);
                Gs(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id)) , SNR_ind) = C  * (H_i * F)';
            end
        end 

        mdl.chn_info.(mu_name).Gs = Gs;
        m.Fs = Fs;
    end
elseif (strcmp(act, 'det'))
    B = mdl.B;
    G = mdl.chn_info.(mu_name).Gs(:, :, m.SNR_ind);
    
    m.s_hat = G * m.y;
end

end
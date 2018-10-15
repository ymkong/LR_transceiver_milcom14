%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% precoding using Block diagonalization (BD) with joint transceiver design. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mdl, m ] = det_bd_j(act, mdl, m, mu_name)

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

        for user_id = 1 : mdl.users
            for SNR_ind = 1 : NS
                sigma = mdl.sigmas(SNR_ind);

                H_i_hat = [H(1 : sum(N_r(1: (user_id-1))), :);H(sum(N_r(1:user_id))+1 : end, :)];
                H_i = H(sum(N_r(1:(user_id-1)))+1 : sum(N_r(1:user_id)), :);

                [UU, SS, VV] = svd(H_i_hat);

                % take the null space of H_i_hat
                VV = VV(:, rank(SS) + 1 : end);
                H_e = H_i * VV;

                % design A, and make sure VV * A meets the power constraints  
                % Es * tr(FF') <= P / users
                delta = ones(B(user_id), 1) * mdl.Es;

                [V, Lambda] = eig((H_e' * H_e) ./ sigma .^ 2);
                lambda = diag(Lambda);
                % Order lambda
                [tmp, ind] = sort(lambda, 'descend');
                lambda = lambda(ind(1 : B(user_id)));
                V = V(:, ind(1 : B(user_id)));

                rho = lambda .* delta;

                for k = B(user_id) : -1 : 1
                    if (rho(k) > eps)
                        mu = (sum(sqrt(delta(1 : k)) ./ sqrt(lambda(1 : k))) / (P / users + sum(1 ./ lambda(1 : k)))) .^ 2;
                        if (mu < rho(k))
                            break;
                        end
                    end
                end

                Theta_A = diag([sqrt(pos_threshold(1 ./ sqrt(mu) ./ sqrt(lambda(1 : k)) ./ sqrt(delta(1 : k))  - 1 ./ lambda(1 : k) ./ delta(1 : k))); zeros(B(user_id) - k, 1)]);

                A = V * Theta_A;
                F = VV * A;
                Fs{SNR_ind} = [Fs{SNR_ind}, F];

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% precoding using GZI with ELR joint ELR-aided transceiver design. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mdl, m ] = det_gzi_lr_j(act, mdl, m, mu_name, lr_func, lr_name)

if (strcmp(act, 'updateH'))
    if(~isfield(mdl.chn_info, mu_name))    
        N_r = mdl.N_r;
        N_t = mdl.N_t;
        B = mdl.B;
        NS = length(mdl.SNRdb);
        
        Gs = zeros(sum(B), sum(N_r), NS);
        Ts = zeros(sum(B), sum(B), NS);
        ITs = zeros(sum(B), sum(B), NS);
        
        P = mdl.P;
        H = m.H;
        Fs = cell(1, NS);
        users = mdl.users;
        
        H_inv = H' * inv(H * H');
        
        for user_id = 1 : mdl.users
            for SNR_ind = 1 : NS
                sigma = mdl.sigmas(SNR_ind);

                %H_i_hat = [H(1 : sum(N_r(1: (user_id-1))), :);H(sum(N_r(1:user_id))+1 : end, :)];
                H_i = H(sum(N_r(1:(user_id-1)))+1 : sum(N_r(1:user_id)), :);
                
                H_i_inv = H_inv(:, sum(N_r(1:(user_id-1)))+1 : sum(N_r(1:user_id)));
                
                [tQ, tR] = qr(H_i_inv, 0);
                %[UU, SS, VV] = svd(H_i_hat);

                % take the null space of H_i_hat
                %VV = VV(:, rank(SS) + 1 : end);
                VV = tQ;
                
                H_e = H_i * VV;

                % design A (equivalent to F in SU-MIMO case)
                % make sure VV * A meets the power constraints, Es * tr(FF') <= P / users
                A = [eye(B(user_id)) * sqrt(P ./ users ./ B(user_id)); zeros(size(H_e, 2) - B(user_id), B(user_id))] ./ sqrt(mdl.Es);

                HF = H_e * A;
                HR = [HF ; eye(B(user_id)) .* sigma ./ sqrt(mdl.Es)];
                [Ht, T] = lr_func(HR);
                
                I_T = inv(T);
                R_zz = (I_T * I_T') * mdl.Es;
                [U, Delta] = eig(R_zz);
                delta = diag(Delta);
                [tmp, ind] = sort(delta, 'descend');
                delta = delta(ind(1 : B(user_id)));
                U = U(:, ind(1 : B(user_id)));
                
                [V, Lambda] = eig((H_e' * H_e) ./ sigma .^ 2);
                lambda = diag(Lambda);
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
                
                A_t = V * Theta_A * U';
                
                A = A_t * I_T;
                
                F = VV * A;
                Fs{SNR_ind} = [Fs{SNR_ind}, F];

                C = inv((H_i * F)' * (H_i * F) + sigma .^ 2 * eye(B(user_id)) ./ mdl.Es);
                Gs(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(N_r(1: (user_id-1)))+1 : sum(N_r(1:user_id)) , SNR_ind) = I_T * C  * (H_i * F)';
                Ts(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)) , SNR_ind) = T;
                ITs(sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)), sum(B(1: (user_id-1)))+1 : sum(B(1:user_id)) , SNR_ind) = I_T;
                
            end
        end 

        mdl.chn_info.(mu_name).Gs = Gs;
        mdl.chn_info.(mu_name).Ts = Ts;
        mdl.chn_info.(mu_name).ITs = ITs;
        
        m.Fs = Fs;
    end
elseif (strcmp(act, 'det'))
    B = mdl.B;
    G = mdl.chn_info.(mu_name).Gs(:, :, m.SNR_ind);
    xc = G * m.y;
    
    T = mdl.chn_info.(mu_name).Ts(: , : , m.SNR_ind);
    I_T = mdl.chn_info.(mu_name).ITs(: , : , m.SNR_ind);
    
    zc_hat = T * round(xc / 2 - I_T * (0.5 * ones(sum(B),1)) - 1j * (I_T * (0.5 * ones(sum(B),1))));
    m.s_hat = 2 * (zc_hat + (0.5 + 1j * 0.5));

    
end

end




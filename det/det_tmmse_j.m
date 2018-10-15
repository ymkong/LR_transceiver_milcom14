%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 

% precoding T-MMSE. 
% minimize the system wide MSE under the total transmit power constraint.
% Total min MSE criterion.
% iterative MMSE method.
% iterativly compute transmit and receive matrix.
% based on total min MSE criterion.
% See details in paper [kong14lr].
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mdl, m ] = det_tmmse_j(act, mdl, m, mu_name)

num_iteration = 150;

if (strcmp(act, 'updateH'))
    if(~isfield(mdl.chn_info, mu_name))    
        N_r = mdl.N_r;
        N_t = mdl.N_t;
        B = mdl.B;
        NS = length(mdl.SNRdb);        
        Gs = zeros(sum(B), sum(N_r), NS); % to store the receive matrices
        P = mdl.P; % total transmit power
        H = m.H; % total channel matrix
        Fs = cell(1, NS); % transmit matrices
        users = mdl.users; % number of users
        
        for SNR_ind = 1 : NS
            sigma = mdl.sigmas(SNR_ind);

            [F, G] = t_mmse(num_iteration, H, P, N_r, N_t, B, users, sigma, mdl.Es);
            
            Fs{SNR_ind} = F;
            Gs(:, :, SNR_ind) = G; % G is a block diagonal matrix
            
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
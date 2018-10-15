%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% This program runs the BER of various precoding schemes in a multi-user MIMO downlinks. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
% Acknowledgement: the code style here is heavily influenced by that of Qi Zhou.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R] = jtxrx_mu_sim(model)

% set random seed
rand('state', model.seed);
randn('state', model.seed);

% copy down simulation parameters
N_r = model.N_r;
N_t = model.N_t;
B = model.B;
sim_n = model.sim_n;
chn_n = model.chn_n;
SNRdb = model.SNRdb;
P = model.P;
algs = model.algs;
h_mod = model.mod;
h_dmod = modem.qamdemod(h_mod, 'OutputType', 'bit');

% number of users
if(isfield(model, 'users'))
    users = model.users;
else
    users = size(B, 2);
end

% compute more parameters
M_c = log2(h_mod.M); % bits per symbol
Es = mean(abs(h_mod.Constellation) .^ 2);
Eb = Es / M_c;
sigmas = sqrt(P ./ (10 .^ (SNRdb / 10))); % noise variance based on transmit SNR
model.M_c = M_c;
model.Es = Es;
model.Eb = Eb;
model.sigmas = sigmas;

% set up detection model
det_models = cell(length(algs), 1);
for alg_ind = 1 : length(algs)
    det_models{alg_ind} = struct('H', [], 'y', [], ...
        'sigma', [], 's_hat', [], 'SNR_ind', 0);
end

% set up error counter and helper variables
errs = zeros(length(SNRdb), length(algs));
start_SNR_ind = ones(length(algs), 1);
blocks_per_SNR = zeros(length(SNRdb), length(algs));

tic
% start simulation
for runs = 1 : sim_n
    H = (randn(sum(N_r), sum(N_t)) + randn(sum(N_r) , sum(N_t)) * 1i) / sqrt(2);
    model.chn_info = struct();
    
    % Generate the transmitted signal (s with Identify correlation matrix)
    bs = round(rand(sum(B) * M_c, 1));
    s = modulate(h_mod, bs);
    
    for alg_ind = 1 : length(algs)
        det_models{alg_ind}.H = H;
        det_models{alg_ind}.s = s;

        % Update channel matrix and obtain precoding matrix
        % Note that power constraint should be satisified, i.e., tr{F * R_ss * F} <= P        
        [model, det_models{alg_ind}] = algs{alg_ind}.func('updateH', model, det_models{alg_ind});
        
        % Check if the power constraint is satisfied
        for SNR_ind = 1 : length(SNRdb)
            P_ind = real(sum(sum(abs(det_models{alg_ind}.Fs{SNR_ind}) .^ 2))) * Es;
            if (P_ind - 1e-3 > P)
                warning('Power constraint is not satisfied!');
            elseif (P - 1e-3 > P_ind)
                warning('Power is under-utilized!');
            end
        end
    end
    
    for chn_runs = 1 : chn_n
        % Generate noise (normalized)
        n = (randn(sum(N_r), 1) + randn(sum(N_r), 1) * 1i) ./ sqrt(2);
        
        for SNR_ind = 1 : length(SNRdb)
            model.algshared = struct();
            for alg_ind = 1 : length(algs)
                if (start_SNR_ind(alg_ind) > SNR_ind)
                    continue
                else
                    blocks_per_SNR(SNR_ind, alg_ind) = blocks_per_SNR(SNR_ind, alg_ind) + 1;
                end
                
                % get noise variance
                sigma = sigmas(SNR_ind);

                % get received signal
                y = H * det_models{alg_ind}.Fs{SNR_ind} * s + n * sigma;                

                % setup detection variables                
                func = algs{alg_ind}.func;
                det_models{alg_ind}.y = y;
                det_models{alg_ind}.sigma = sigma;
                det_models{alg_ind}.SNR_ind = SNR_ind;
                det_models{alg_ind}.s = s;
                
                % Call the detection algorithms
                [model, det_models{alg_ind}] = func('det', model, det_models{alg_ind});
                
                % demodulate detected symbols
                b_hat = demodulate(h_dmod, det_models{alg_ind}.s_hat);
                
                % compute errors
                errs(SNR_ind, alg_ind) = errs(SNR_ind, alg_ind) + sum(b_hat ~= real(bs));
            end
        end
        
        % if collected errors for a certain SNR is enough, we do not have to run that SNR next time
        for alg_ind = 1 : length(algs)
            while (start_SNR_ind(alg_ind) <= length(SNRdb)) && (errs(start_SNR_ind(alg_ind), alg_ind) > model.max_no_error)
                start_SNR_ind(alg_ind) = start_SNR_ind(alg_ind) + 1;
            end
        end
    end

    % print out timing info
    if (mod(runs * chn_n, 500) == 0)
        fprintf('Iter: %d, Time: %f, TPI: %f, Remaining: %0.2f\n',...
            runs * chn_n, toc, toc / runs / chn_n, (sim_n - runs) * (toc / runs));
    end
    
end
toc

% return BER
Pe = errs ./ sum(B) ./ blocks_per_SNR./ M_c;
R = struct('Pe', Pe, 'errs', errs);

end

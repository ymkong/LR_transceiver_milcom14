%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% This program compares BER of various precoding schemes in a multi-user MIMO downlinks. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
addpath('./det', './lr', './helpers');

% setup model parameters
model = struct('seed', sum(clock .* 100), ...
    'N_r', [4, 4, 4], ... % receive antenna configuration
    'N_t', 12, ... % % number of transmit antennas
    'B', [4, 4, 4], ... % user streams configuration. If '-maxsnr' schemes and t-mmse are used, B = [1, 1, 1]
    'sim_n', 5e4, ... % number of simulations
    'chn_n', 1e1, ...
    'SNRdb', 15 : 3 : 33, ...
    'P', 1, ... % total transmit power
    'algs', [], ...
    'max_no_error', 1e3, ...
    'mod', modem.qammod('M', 4, 'SymbolOrder', 'gray', 'InputType', 'bit')); % If '-maxsnr' schemes and t-mmse are used, M needs to set accordingly (e.g., M = 256)

% number of users
users = size(model.B, 2);
model.users = users;

% setup transceiver algorithms
algs = {};
% BD
algs{length(algs) + 1} = struct('sn', 'BD', 'title', 'BD', 'func', @(act, mdl, m) det_bd(act, mdl, m, 'bd'), 'marker', 'bd-');
algs{length(algs) + 1} = struct('sn', 'BD-J', 'title', 'BD-J', 'func', @(act, mdl, m) det_bd_j(act, mdl, m, 'bd_j'), 'marker', 'gs-');
algs{length(algs) + 1} = struct('sn', 'BD-LR-J', 'title', 'BD-LR-J', 'func', @(act, mdl, m) det_bd_lr_j(act, mdl, m, 'bd_lr_j', @(H) elr_dual_c(H), 'elr'), 'marker', 'm<-');
% GZI
algs{length(algs) + 1} = struct('sn', 'GZI-LR-J', 'title', 'GZI-LR-J', 'func', @(act, mdl, m) det_gzi_lr_j(act, mdl, m, 'gzi_lr_j', @(H) elr_dual_c(H), 'elr'), 'marker', 'r*-');

% % S-MMSE
% algs{length(algs) + 1} = struct('sn', 'S-MMSE-max-SNR', 'title', 'S-MMSE-max-SNR', 'func', @(act, mdl, m) det_s_mmse_maxsnr(act, mdl, m, 's_mmse_maxsnr'), 'marker', 'k*-');
% % total MMSE (T-MMSE)
% algs{length(algs) + 1} = struct('sn', 'T-MMSE', 'title', 'T-MMSE', 'func', @(act, mdl, m) det_tmmse_j(act, mdl, m, 'tmmse_j'), 'marker', 'bd-');
% % RBD
% algs{length(algs) + 1} = struct('sn', 'RBD-max-SNR', 'title', 'RBD-max-SNR', 'func', @(act, mdl, m) det_rbd_maxsnr(act, mdl, m, 'rbd_maxsnr'), 'marker', 'r^-');

model.algs = algs;

% print a string that summarizes the current simulation paramters
fn = model_fn(model)

% Execute the simulator
[R] = jtxrx_mu_sim(model);

% plot the BER of various transceiver designs
algs = model.algs;
Pe = R.Pe;
figure
legends = cell(length(algs), 1);
for i = 1 : length(algs)
    semilogy(model.SNRdb, Pe(: , i), algs{i}.marker, 'LineWidth', 1.5)
    hold on
    legends{i} = algs{i}.title;
end
legend(legends)
axis([min(model.SNRdb), max(model.SNRdb), 1e-5, 1]);
xlabel('SNR (dB)', 'FontSize' , 12);
ylabel('BER', 'FontSize' , 12);

% name the figure
N_t = model.N_t;
N_r = model.N_r;
B = model.B;
if(users == 3)
    simName = sprintf('Nt%d-Nr(%d,%d,%d)-r(%d,%d,%d)-maxerror %0.1e %dQAM iter=%0.1e',N_t, ...
        N_r(1),N_r(2),N_r(3),B(1),B(2),B(3), model.max_no_error, model.mod.M, model.sim_n);
elseif(users == 2)
    simName = sprintf('Nt%d-Nr(%d,%d)-r(%d,%d)-maxerror %0.1e %dQAM iter=%0.1e',N_t, ...
        N_r(1),N_r(2),B(1),B(2), model.max_no_error, model.mod.M, model.sim_n);
elseif(users == 1)
    simName = sprintf('Nt%d-Nr(%d)-r(%d)-maxerror %0.1e %dQAM iter=%0.1e',N_t, ...
        N_r(1), B(1), model.max_no_error, model.mod.M, model.sim_n);
else
    disp('naming error');
end

% save the figure
title(simName)
saveas(gcf, getFigName(simName));
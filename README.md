# LR_transceiver_milcom14

This repository contains necessary MATLAB codes to reproduce BER results of LR-aided transceiver designs from [1].

To run results with BD, BD-J, BD-LR-J, GZI-LR-J,
run main_jtxrx_mu.m.

To run and compare results with S-MMSE-max-SNR, RBD-max-SNR, T-MMSE, set
model.B = [1, 1, 1];
model.mod.M = 256;
then run main_jtxrx_mu.m.

[1] Y. Kong, Q. Zhou, and X. Ma, “Lattice reduction aided transceiver design for MU MIMO downlink transmissions,” in Proc. Military Commu. Conf., Baltimore, MD, Oct 2014, pp. 556-562.

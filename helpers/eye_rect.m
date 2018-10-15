%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong14lr] Lattice reduction aided transceiver design for MU MIMO downlink transmissions 
% Get corresponding identity matrix for T-MMSE scheme. See details in paper [kong14lr]
% 
% Written by: Yiming Kong
% Date: 3/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = eye_rect(a, b)

I = zeros(a, b);

if(a <= b)
    I(1 : a, 1 : a) = eye(a);
else
    I(1 : b, 1 : b) = eye(b);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELR wrapper function
% Written by: Qi Zhou
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ H_t T info] = elr_dual_c( H, in_info )

N = size(H , 2);

if (nargin >= 2 && isfield(in_info, 'H'))
    G = in_info.G;
else
    G = (H' * H);
end


if (nargin >= 2 && isfield(in_info, 'C'))
    C = in_info.C;
else
    C = inv(G);
end

if (nargin >= 2 && isfield(in_info, 'strategy'))
    s = in_info.strategy;
else
    s = 0;
end

if (nargout >= 3)
    [T info] = elr_dual_core_c(C, eye(N), s);
    info.pre_C = C;
    info.pre_G = G;
else
    [T] = elr_dual_core_c(C, eye(N), s);
end

H_t = H * T;

end


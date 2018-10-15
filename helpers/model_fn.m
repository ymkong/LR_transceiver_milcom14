function [ fn ] = model_fn( model )
%MODEL_FN Summary of this function goes here
%   Detailed explanation goes here
pcn = {'npc', 'pc'};

algn = '';
for i = 1 : length(model.algs)
    algn = [algn '(' model.algs{i}.sn ')'];
end

if (isfield(model, 'Rtx'))
    Rtxs = sprintf('%02.0f', 100 * abs(model.Rtx(1, 2)));
else
    Rtxs = '00';
end

if (isfield(model, 'Rrx'))
    Rrxs = sprintf('%02.0f', 100 * abs(model.Rtx(1, 2)));
else
    Rrxs = '00';
end



sname = 'JTxRx';

fn = sprintf('%s_Nt%d_Nr%d_B%d_Rtx%s_Rrx%s_%s_QAM%d_n%s_c%s_%s', sname, ...
    model.N_t, model.N_r, model.B, Rtxs, Rrxs, algn, model.mod.M, get_sn(model.sim_n), get_sn(model.chn_n), ...
    datestr(now, 'yymmddTHHMM'));
end

function [sn] = get_sn(d)
v = num2str(d);
sn = [v(1) num2str(length(v) - 1)];
end


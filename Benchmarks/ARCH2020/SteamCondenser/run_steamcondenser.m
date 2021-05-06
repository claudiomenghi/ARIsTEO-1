% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function [tout, yout] = run_steamcondenser(u, T)
    ts = u(:,1);
    us = u(:,2:end);

    tin = 0:0.01:T;
    xin = interp1(ts, us, tin, 'previous');
    u = [tin' xin];

    assignin('base','u',u);
    assignin('base','T',T);

    result = sim('steamcondense_RNN_22_BR', ...
        'StopTime', 'T', ...
        'LoadExternalInput', 'on', 'ExternalInput', 'u', ...
        'SaveTime', 'on', 'TimeSaveName', 'tout', ...
        'SaveOutput', 'on', 'OutputSaveName', 'yout', ...
        'SaveFormat', 'Array');
    tout = result.tout;
    yout = result.yout;
end

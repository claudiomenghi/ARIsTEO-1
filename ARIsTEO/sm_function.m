% Copyright 2021 University of Luxembourg
%
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function [T, XT, YT, LT, CLG, Guards] = sm_function(X0,simT,TU,U)
% SM_FUNCTION encapsulates the surrogatemodel in a function handle
% the surrogate model format is not supported in staliro
% To run staliro on the surrogate model, we call the model as a function handle
% called sm_function.m which simulates the surrogate model
% this function replaces the old function @blackbox_exec_identified_model
% INPUTS:
% - X0: the initial conditions values if any.
% - simT: teh simulation time.
% - TU: the array of time.
% - U: the input signals.
% OUTPUTS:

%  - T      = col' vector of time stamps of simulated trajectory,
%  - XT     = internal states at the time stamps,
%  - YT     = outputs at the time stamps,
%  - LT     = locations at the time stamps
%  - Guards    = guard structure,
%  - CLG = Control Location Graph.
%         Some of these will be empty, depending on the type of system being simulated. E.g. if it's a simulink model,
%         at least CLG, LT and GRD will be empty.
global abstractedmodel;
XT=[];
LT=[];
CLG=[];
Guards=[];
if( ~isempty(X0))
    opt = simOptions('InitialCondition', X0);
    simOut=sim(abstractedmodel,U, opt);
else
    simOut=sim(abstractedmodel,U);
end

v=ver('Matlab');
if(isequal(v.Version,'7.11'))
    YT=simOut.OutputData;
    T=simOut.SamplingInstants;
else
    YT=simOut;
    T=TU;
end


end

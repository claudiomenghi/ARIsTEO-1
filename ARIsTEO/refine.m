% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function [data, abstractedmodel]=refine(data, abstractedmodel, aristeo_options, idOptions)
% REFINE performs the refinement step to build a surrogate model using the
% previously generated data and SM
% INPUTS:
% - data: the previously generated data
% - abstractedmodel: the surrogate model built in the previous iteration.
% - aristeo_options: aristeo_options options. opt should be of type "aristeo_options".
%       If the default options are going to be used, then this input may be
%       omitted. For instructions on how to change Aristeo options,
%       see the aristeo_options help file for each desired property.
% - idOptions: the system identification setting options
% OUTPUTs:
% - data: the data after reconstruction inputs and outputs data.
% - abstractedmodel: the refined surrogate model
    global absreftime;
    absreftimetic=tic;
    data = misdata(data);
    try
        if (isempty(idOptions))
            abstractedmodel=pem(iddata(data.y{size(data.y,2)},data.u{size(data.u,2)},aristeo_options.SampTime),abstractedmodel);
        else
            abstractedmodel=pem(iddata(data.y{size(data.y,2)},data.u{size(data.u,2)},aristeo_options.SampTime),abstractedmodel,idOptions);
        end
    catch e
        warning(e.message);
    end
    absreftime=toc(absreftimetic);
end

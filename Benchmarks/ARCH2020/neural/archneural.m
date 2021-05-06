% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

% speed->1 RPM->2 gear->3

function []=archneural(modelsstructure,currentOrder)

    clearvars -except modelsstructure currentOrder

    opt = aristeo_options();
    opt.runs = 50;
    opt.optim_params.n_tests = 300; %ntest will specifies the rest of tests

    assignin('base','u_ts',0.001)
    NN   = '[]_[1,18] (!(pred)->(<>_[0,2]([]_[0,1] (pred))))';

    preds(1).str='pred';
    preds(1).A = 1;
    preds(1).b = 0;

    properties={NN};

    propertiesString={'NN'};


    for i=1:size(properties,2)
        disp(properties{i})

        phi=properties{i};

        opt.abstraction_algorithm=modelsstructure;
        opt.n_refinement_rounds=300;
        opt.dispinfo=0;

        opt.nx=str2double(currentOrder);
        ninputs=1;
        noutputs=1;

        opt.na=ones(noutputs,noutputs)*str2double(currentOrder);
        opt.nb=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nf=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nk=ones(noutputs,ninputs)*0;
        opt.nc=ones(noutputs,1)*str2double(currentOrder);
        opt.nd=ones(noutputs,1)*str2double(currentOrder);
        opt.np=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nz=ones(noutputs,ninputs)*str2double(currentOrder);

        opt.interpolationtype={'pconst'};

        opt.TotSimTime=40;
        init_cond = [];

        input_range = [1 3];
        cp_array = 3;
        model = 'nn';
        warning off

        [results,history] = aristeo(model,init_cond,input_range,cp_array,phi,preds,opt.TotSimTime,opt);

        save(strcat('ARIsTEO_Arch19Bench','_','NN','_',propertiesString{i},'_',opt.abstraction_algorithm,currentOrder),'results','history');

    end
    exit();
end

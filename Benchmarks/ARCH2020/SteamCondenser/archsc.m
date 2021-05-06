% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

% speed->1 RPM->2 gear->3

function []=archsc(modelsstructure,currentOrder)


    opt = aristeo_options();
    opt.runs = 50;
    opt.optim_params.n_tests = 300; %ntest will specifies the rest of tests



    SC = '[]_[30,35] (p1/\p2)';

    preds(1).str = 'p1';
    preds(1).A =  [0 0 0 1];
    preds(1).b =  87.5;

    preds(2).str = 'p2';
    preds(2).A =  [0 0 0 -1];
    preds(2).b =  -87;

    properties={SC};

    propertiesString={'SC'};


    for i=1:size(properties,2)
        disp(properties{i})

        phi=properties{i};

        opt.abstraction_algorithm=modelsstructure;
        opt.n_refinement_rounds=300;
        opt.dispinfo=0;

        opt.nx=str2double(currentOrder);
        ninputs=1;
        noutputs=4;

        opt.na=ones(noutputs,noutputs)*str2double(currentOrder);
        opt.nb=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nf=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nk=ones(noutputs,ninputs)*0;
        opt.nc=ones(noutputs,1)*str2double(currentOrder);
        opt.nd=ones(noutputs,1)*str2double(currentOrder);
        opt.np=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nz=ones(noutputs,ninputs)*str2double(currentOrder);

        opt.interpolationtype={'pconst'};


        opt.TotSimTime=100;
        init_cond = [];

        input_range = [3.99 ,4.01];
        cp_array = 12;
        model = 'stc';
        warning off

	opt.runs = 50;
        [results,history] = aristeo(model,init_cond,input_range,cp_array,phi,preds,opt.TotSimTime,opt);

        save(strcat('ARIsTEO_Arch19Bench','_','SC','_',propertiesString{i},'_',opt.abstraction_algorithm,currentOrder),'results','history');

    end
    exit();
end

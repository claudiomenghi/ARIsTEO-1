% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

% speed->1 RPM->2 gear->3

function []=archcarscc4(modelsstructure,currentOrder)

    clearvars -except modelsstructure currentOrder

    opt = aristeo_options();
    opt.runs = 50;
    opt.optim_params.n_tests = 300; %ntest will specifies the rest of tests
    cp_array = [7,3];

    cc1   = '[]_[0,100] (cc1pred)';
    cc2   = '[]_[0,70] (<>_[0,30] (cc2pred))';
    cc3   = '[]_[0,80] (([]_[0,20] (cc3pred1))  \/ (<>_[0,20] (cc3pred1)))';
    cc4   = '[]_[0,65] (<>_[0,30] ([]_[0,5] (cc4pred1)))';
    cc5   = '[]_[0,72] (<>_[0,8] (([]_[0,5] (cc5pred1)) -> ([]_(5,20) (cc5pred1))))';

    preds(1).str='cc1pred';
    preds(1).A = [0 0 0 1 -1];
    preds(1).b = 40;

    preds(2).str='cc2pred';
    preds(2).A = [0 0 0 -1 1];
    preds(2).b = -15;

    preds(3).str='cc3pred1';
    preds(3).A = [1 -1 0 0 0];
    preds(3).b = 20;

    preds(4).str='cc3pred2';
    preds(4).A = [0 0 0 1 -1];
    preds(4).b = 40;

    preds(5).str='cc4pred1';
    preds(5).A = [0 0 0 -1 1];
    preds(5).b = -8;

    preds(6).str='cc5pred1';
    preds(6).A = [-1 1 0 0 0];
    preds(6).b = -9;

    preds(7).str='cc5pred2';
    preds(7).A = [0 0 0 -1 1];
    preds(7).b = -9;


    %properties={cc1,cc2,cc3,cc4,cc5};

    properties={cc4};
    %propertiesString={'cc1','cc2','cc3','cc4','cc5'};
    propertiesString={'cc4'};

    for i=1:size(properties,2)
        disp(properties{i})

        phi=properties{i};

        opt.abstraction_algorithm=modelsstructure;
        opt.n_refinement_rounds=300;
        opt.dispinfo=0;

        opt.nx=str2double(currentOrder);
        ninputs=2;
        noutputs=5;

        opt.na=ones(noutputs,noutputs)*str2double(currentOrder);
        opt.nb=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nf=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nk=ones(noutputs,ninputs)*0;
        opt.nc=ones(noutputs,1)*str2double(currentOrder);
        opt.nd=ones(noutputs,1)*str2double(currentOrder);
        opt.np=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nz=ones(noutputs,ninputs)*str2double(currentOrder);

        opt.interpolationtype={'pconst','pconst'};

        init_cond = [];

        input_range = [0.0 1.0; 0.0 1.0];
        model = 'carsmodel';
        opt.TotSimTime= 30;
        warning off

        [results,history] = aristeo(model,init_cond,input_range,cp_array,phi,preds,opt.TotSimTime,opt);
        %[results,history] = staliro(model,init_cond,input_range,cp_array,phi,preds,tf,opt);

        save(strcat('ARIsTEO_Arch19Bench','_','CC','_',propertiesString{i},'_',opt.abstraction_algorithm,currentOrder),'results','history');

    end
end

% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

% speed->1 RPM->2 gear->3

function []=archpt(modelsstructure,currentOrder)
    disp(class(currentOrder));
    clearvars -except modelsstructure currentOrder

    opt = aristeo_options();
    opt.runs = 50;
    opt.optim_params.n_tests = 300; %ntest will specifies the rest of tests


    cp_array = [1,10];

    i=0;
    eta = 1;
    % parameter h used for event definition
    h = 0.05;
    % parameter related to the period of the pulse signal
    zeta_min = 5;

    Ut = 0.008;

    low=8.8;
    high=40;

    taus = 10 + eta;
    i = i+1;
    preds(i).str = 'low'; % for the pedal input signal
    preds(i).A =  [0 0 1] ;
    preds(i).b =  low ;
    i = i+1;

    preds(i).str = 'high'; % for the pedal input signal
    preds(i).A =  [0 0 -1] ;
    preds(i).b =  -high ;
    i = i+1;
    % rise event is represented as low/\<>_(0,h)high
    % fall event is represented as high/\<>_(0,h)low
    preds(i).str = 'norm'; % mode < 0.5 (normal mode = 0)
    preds(i).A =  [0 1 0] ;
    preds(i).b =  0.5 ;
    i = i+1;
    preds(i).str = 'pwr'; % mode >0.5 (power mode = 1)
    preds(i).A =  [0 -1 0] ;
    preds(i).b =  -0.5 ;
    i = i+1;
    preds(i).str = 'utr'; % u<=Ut
    preds(i).A =  [1 0 0] ;
    preds(i).b =  Ut ;
    i = i+1;
    preds(i).str = 'utl'; % u>=-Ut
    preds(i).A =  [-1 0 0] ;
    preds(i).b =  Ut ;
    i = i+1;

    AFC27 = ['[]_[11,50](!((low/\<>_(0,' num2str(h) ')high) \/ (high/\<>_(0,' num2str(h) ')low)))'];
    AFC29 = '[]_[11,50](utr /\ utl)';
    AFC33 = '[]_[11,50](utr /\ utl)';

    properties={AFC27,AFC29,AFC33};

    propertiesString={'AFC27','AFC29','AFC33'};

    global simTime
    simTime=50;

    assignin('base','simTime',50)
    assignin('base','en_speed',1000)
    assignin('base','measureTime',1)
    assignin('base','spec_num',1)
    assignin('base','fuel_inj_tol',1)
    assignin('base','MAF_sensor_tol',1)
    assignin('base','AF_sensor_tol',1)
    assignin('base','sim_time',50)

    assignin('base','en_speed',1000)
    assignin('base','measureTime',1)
    assignin('base','fault_time',60)
    assignin('base','fault_time',60)


    for i=1:size(properties,2)
        disp(properties{i})

        phi=properties{i};

        opt.abstraction_algorithm=modelsstructure;
        opt.n_refinement_rounds=300;
        opt.dispinfo=0;

        opt.nx=str2double(currentOrder);
        ninputs=2;
        noutputs=3;

        opt.na=ones(noutputs,noutputs)*str2double(currentOrder);
        opt.nb=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nf=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nk=ones(noutputs,ninputs)*0;
        opt.nc=ones(noutputs,1)*str2double(currentOrder);
        opt.nd=ones(noutputs,1)*str2double(currentOrder);
        opt.np=ones(noutputs,ninputs)*str2double(currentOrder);
        opt.nz=ones(noutputs,ninputs)*str2double(currentOrder);

        opt.interpolationtype={'pconst','pconst'};
        opt.TotSimTime=50;
        init_cond = [];

        input_range = [900  1100; 0 61.1];
        phi
        if (strcmp(phi,'AFC33')==1)
            input_range = [900  1100; 61.2 81.2];
        end
        model = 'AbstractFuelControl_M1_2018a';
       	opt.TotSimTime = 30;
        warning off

        [results,history] = aristeo(model,init_cond,input_range,cp_array,phi,preds,opt.TotSimTime,opt);
        %[results,history] = staliro(model,init_cond,input_range,cp_array,phi,preds,tf,opt);

        save(strcat('ARIsTEO_Arch19Bench','_','PT','_',propertiesString{i},'_',opt.abstraction_algorithm,currentOrder),'results','history');

    end
    exit();
end

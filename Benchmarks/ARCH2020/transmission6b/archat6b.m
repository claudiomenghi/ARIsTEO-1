% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

% speed->1 RPM->2 gear->3

function []=archat6b(modelsstructure,currentOrder)

    clearvars -except modelsstructure currentOrder

    opt = aristeo_options();
    opt.runs = 50;
    opt.optim_params.n_tests = 300; %ntest will specifies the rest of tests
    cp_array = [10,10];



    at1 = '[]_[0,20] speed120';
    at2 = '[]_[0,10] rpm4750';
    at51 = '[]_[0, 30] ((!(gearleq1 /\ geargeq1) /\ <>_[0.001,0.1] (gearleq1 /\ geargeq1))-> <>_[0.001,0.1] []_[0.0,2.5] (gearleq1 /\ geargeq1))';
    at52 = '[]_[0, 30] ((!(gearleq2 /\ geargeq2) /\ <>_[0.001,0.1] (gearleq2 /\ geargeq2))-> <>_[0.001,0.1] []_[0.0,2.5] (gearleq2 /\ geargeq2))';
    at53 = '[]_[0, 30] ((!(gearleq3 /\ geargeq3) /\ <>_[0.001,0.1] (gearleq3 /\ geargeq3))-> <>_[0.001,0.1] []_[0.0,2.5] (gearleq3 /\ geargeq3))';
    at54 = '[]_[0, 30] ((!(gearleq4 /\ geargeq4) /\ <>_[0.001,0.1] (gearleq4 /\ geargeq4))-> <>_[0.001,0.1] []_[0.0,2.5] (gearleq4 /\ geargeq4))';


    at6a= ('[]_[0.0, 30.0] (rpmlt3000) -> []_[0.0, 4.0] (speed35))');
    at6b= ('[]_[0.0, 30.0] (rpmlt3000) -> []_[0.0, 8.0] (speed50))');
    at6c= ('[]_[0.0, 30.0] (rpmlt3000) -> []_[0.0, 20.0] (speed65))');


    preds(1).str='speed120';
    preds(1).A = [1 0 0];
    preds(1).b = 120;

    preds(2).str='rpm4750';
    preds(2).A = [0 1 0];
    preds(2).b = 4750;


    preds(3).str='gearleq1';
    preds(3).A = [0 0 1];
    preds(3).b = 1;

    preds(4).str='geargeq1';
    preds(4).A = [0 0 -1];
    preds(4).b = -1;

    preds(5).str='gearleq2';
    preds(5).A = [0 0 1];
    preds(5).b = 2;

    preds(6).str='geargeq2';
    preds(6).A = [0 0 -1];
    preds(6).b = -2;

    preds(7).str='gearleq3';
    preds(7).A = [0 0 1];
    preds(7).b = 3;

    preds(8).str='geargeq3';
    preds(8).A = [0 0 -1];
    preds(8).b = -3;


    preds(9).str='gearleq4';
    preds(9).A = [0 0 1];
    preds(9).b = 4;

    preds(10).str='geargeq4';
    preds(10).A = [0 0 -1];
    preds(10).b = -4;

    preds(11).str='speed35';
    preds(11).A = [1 0 0];
    preds(11).b = 35;

    preds(12).str='speed35';
    preds(12).A = [1 0 0];
    preds(12).b = 35;

    preds(13).str='speed50';
    preds(13).A = [1 0 0];
    preds(13).b = 35;

    preds(14).str='speed65';
    preds(14).A = [0 1 0];
    preds(14).b = 65;

    preds(15).str='rpmlt3000';
    preds(15).A = [0 1 0];
    preds(15).b = 3000;

    %properties={at1,at2,at6a,at6b,at6c};

properties={at6b};
    %propertiesString={'at1','at2','at6a','at6b','at6c'};
propertiesString={'at6b'};

    for i=1:size(properties,2)
        disp(properties{i})

        phi=properties{i};

        opt.abstraction_algorithm=modelsstructure;
        opt.n_refinement_rounds=300;
        opt.dispinfo=0;
        opt.SampTime=0.04;

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

        opt.interpolationtype={'pchip','pchip'};

        init_cond = [];

        input_range = [0 100;0 350];
        model = 'Autotrans_shift';
        opt.TotSimTime=50;
        warning off

        [results,history] = aristeo(model,init_cond,input_range,cp_array,phi,preds,opt.TotSimTime,opt);
        %[results,history] = staliro(model,init_cond,input_range,cp_array,phi,preds,tf,opt);

        save(strcat('ARIsTEO_Arch19Bench','_','Autotrans_shift','_',propertiesString{i},'_',opt.abstraction_algorithm,currentOrder),'results','history');

    end
end

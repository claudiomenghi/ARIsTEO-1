% Copyright 2021 University of Luxembourg
%
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

% defining a variable that contains the name of the model
model='pendulum';

% As we aim to consider the initial conditions of the model, we set init_cond to an empty set
init_cond = [];

% input_range we assume that the user can apply a momentun in the range [-0.5 0.5]
input_range = [-2 2];

% we want inputs with 100 control points
cp_array = 100;

% we want the pendulum to remain below the orizontal line
phi='[]_[0,10000] (a/\b)';

preds(1).str = 'a';
preds(1).A = [1];
preds(1).b = [1.5];

preds(2).str = 'b';
preds(2).A = [-1];
preds(2).b = [1.5];


TotSimTime=10;

opt=aristeo_options();
opt.n_refinement_rounds=9;

opt.fals_at_zero=0;
opt.abstraction_algorithm='ss';
opt.nx=2;
opt.optim_params.n_tests=10;
opt.dispinfo=1;

% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

t__ = [0; 10; 20; 30];
u__ = [50 0; 100 0; 0 325; 0 325];

u = [t__, u__];
T = 30;

[tout, yout] = run_transmission(u, T);

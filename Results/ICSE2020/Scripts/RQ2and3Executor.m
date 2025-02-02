% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

curPath=fileparts(which('RQ1.m'));
mainpath = strrep(curPath,'RQs','');
addpath(genpath(mainpath));

for i=[7 9 12 14 16 19]
    RQ2andRQ3('staliro',i);
end
for i=[5 7 9 11 13 15]
	RQ2andRQ3('aristeo',i);
end

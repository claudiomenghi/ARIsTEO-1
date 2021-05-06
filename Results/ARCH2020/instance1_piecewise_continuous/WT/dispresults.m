% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

k='2';
j={'arx'};
prefix='ARIsTEO_Arch19Bench_pchip_WT_';

for p={'WT1','WT2','WT3','WT4'}%p={'WT1','WT2','WT3'}

    config=strcat(j,num2str(k));
    file=strcat(prefix,p,'_',config,'.mat');
    load(file{1})
    disp('FR     mean(S)     median(S)');
    disp(strcat(num2str(sum([results.run.bestRob]<0)),{'     '},num2str(mean([results.run.nTests].*[[results.run.bestRob]<0])),...
         {'     '},num2str(median([results.run.nTests].*[[results.run.bestRob]<0]))));

end

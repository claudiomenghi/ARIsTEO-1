% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

k='2';
j={'arx'};
prefix='ARIsTEO_Arch19Bench_NN_';

for p={'NN'}

    config=strcat(j,num2str(k));
    file=strcat(prefix,p,'_',config,'.mat');
    load(file{1})
    disp('FR     mean(S)     median(S)');
    disp(strcat(num2str(sum([results.run.bestRob]<0)),{'     '},num2str(mean([results.run.nTests].*[[results.run.bestRob]<0])),...
         {'     '},num2str(median([results.run.nTests].*[[results.run.bestRob]<0]))));

end

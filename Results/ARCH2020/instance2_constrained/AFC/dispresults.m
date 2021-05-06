% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

j={'arx2','arx2','arx2'};
prefix='ARIsTEO_Arch19Bench_PT_';

i=1;
for  p={'AFC27','AFC29','AFC33'}

    config=strcat(j{i});
    file=strcat(prefix,p,'_',config,'.mat');
    load(file{1})
    disp('FR     mean(S)     median(S)');
    disp(strcat(num2str(sum([results.run.bestRob]<0)),{'     '},num2str(mean([results.run.nTests].*[[results.run.bestRob]<0])),...
         {'     '},num2str(median([results.run.nTests].*[[results.run.bestRob]<0]))));
    i=i+1;
end

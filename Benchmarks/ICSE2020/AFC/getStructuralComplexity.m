% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

close all
clc
load_system('AbstractFuelControl_M1')
a=find_system;

structuralcomplexity=0;
for i=2:1:size(a,1)
   t=a(i);
   t{1}
  bdroot

   if(strfind(t{1},'/'))
       type=get_param(a(i),'BlockType');
       if(~isequal(type{1},'Outport'))

           ports=get_param(a(i),'PortHandles');


           structuralcomplexity=structuralcomplexity+size(ports{1}.Outport,2)*size(ports{1}.Outport,2);
       end
   end
end
structuralcomplexity=structuralcomplexity./(size(a,1)-1);
disp(structuralcomplexity);

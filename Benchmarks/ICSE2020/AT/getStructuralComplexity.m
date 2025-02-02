% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt
 
load_system('sldemo_autotrans_mod01')
a=find_system;

structuralcomplexity=0;
for i=2:1:size(a,1)


   type=get_param(a(i),'BlockType');
   if(~isequal(type{1},'Outport'))

       ports=get_param(a(i),'PortHandles');


       structuralcomplexity=structuralcomplexity+size(ports{1}.Outport,2)*size(ports{1}.Outport,2);
   end

end
structuralcomplexity=structuralcomplexity./(size(a,1)-1);
disp(structuralcomplexity);

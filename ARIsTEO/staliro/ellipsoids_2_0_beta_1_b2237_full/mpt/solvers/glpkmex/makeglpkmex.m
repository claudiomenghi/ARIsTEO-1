% Use this function to build your own MEX interface
% GLPKMEX MEX Interface for the GLPK Callable Library
%
% Copyright (C) 2004-2007, Nicolo' Giorgetti, All right reserved.
% E-mail: <giorgetti __ at __ ieee.org>.
% 
% This file is part of GLPKMEX.
% 
% GLPKMEX is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% This part of code is distributed with the FURTHER condition that it is 
% possible to link it to the Matlab libraries and/or use it inside the Matlab 
% environment.
%
% GLPKMEX is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
% License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GLPKMEX; see the file COPYING. If not, write to the Free
% Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
% 02111-1307, USA.

disp('GLPKMEX - A Matlab interface for GLPK. Script installer. ');
disp('Version 2.4 compatible with GLPK 4.20 (or higher)');
disp('(C) 2001-2007, Nicolo'' Giorgetti.');
disp(' ');

usegraph=0;
reply = input('Do you want to use graphic installer? Y/N [Y]: ', 's');
if isempty(reply) || strcmp(reply,'Y') || strcmp(reply,'y')
   usegraph=1;
else
   usegraph=0;
end

% Provide GLPK path
if usegraph
   GLPKpath=uigetdir(pwd,'Please specify where GLPK is located');
else
   GLPKpath = input('GLPK path: ', 's');
end
if ~ischar(GLPKpath)
   fprintf('In order to compile GLPKMEX you have to specify a correct path of GLPK\n');
   return;
end
fprintf('GLPK path... %s\n',GLPKpath);

%% Check if all the data needed are available
include=[GLPKpath '/include'];
while 1
   d=dir(include);
   if isempty(d)
      if usegraph
         uiwait(msgbox('Specified directory does not contain GLPK include directory!','Include files','error'));
         include=uigetdir(pwd,'Location where GLPK include files are located');
      else
         fprintf('\t Warning: Specified directory does not contain GLPK include directory!\n');
         include = input('GLPK include directory: ', 's');
      end
   else
      break;
   end
end
fprintf('GLPK include files...');
d=dir([include, '/glpk.h']);
if isempty(d)
   fprintf('NO\n');
   fprintf('\tError: glpk.h not found\n');
   return;
else
   fprintf('OK\n');
end

%% Check library file
library=[GLPKpath '/lib/libglpk.a'];
while 1
   d=dir(library);
   if isempty(d)
      if usegraph
         uiwait(msgbox('Specified directory does not contain GLPK libglpk.a library!','Lib file','error'));
         [filename, pathname] = uigetfile('*.a','Choose libglpk.a file');
         library=[pathname,filename];
      else
         fprintf('\t Warning: Specified directory does not contain GLPK libglpk.a library!\n');
         library = input('Specify full GLPK libglpk.a path: ', 's');
      end
   else
      break;
   end
end
fprintf('GLPK library file...');
fprintf('OK\n');

filename='glpkcc.cpp';

cmd=['-I' include ' ' filename ' ' library];

if ispc
   reply = input('Do you compile GLPKMEX with CYGWIN? Y/N [Y]: ', 's');
   if isempty(reply) || strcmp(reply,'Y') || strcmp(reply,'y')
      [filename, pathname] = uigetfile('*.bat', 'Pick a MEXOPTS.BAT file');
      cmd=['-f ', pathname, filename,' ',cmd];      
   end
end

eval(['mex ' cmd]);




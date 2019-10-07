function [sys,x0,str,ts,simStateCompliance] = sldemo_tanksfun(t,x,u,flag,inflow,outflow, ...
                                           lolim,hilim,area)
%SLDEMO_TANKSFUN example animation S-function for model sldemo_tank
%
%   SLDEMO_TANKSFUN updates the graphical representation of the tank for
%   the Simulink example sldemo_tank.
%
%   See also sldemo_tank, sldemo_tankgui.
    
%   Copyright 1990-2012 The MathWorks, Inc.

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u,inflow,outflow,lolim,hilim,area);
  
  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 1, 3, 4, 9 }
    sys=[];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
    otherwise
    error(message('simdemos:sldemo_tanksfun:UnhandledFlag', flag));

end

% end sldemo_tanksfun

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 0;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% specify that the simState for this s-function is same as the default
 simStateCompliance = 'DefaultSimState';
% end mdlInitializeSizes

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u,inflow,outflow,lolim,hilim,area) %#ok

%
% initialize the return arg
%
sys=[];

invalve=u(1);
tankheight=u(2);
outvalve=u(3);

figs = sldemo_tankgui('figures');
if ~ishghandle(figs(1),'figure'),
  sldemo_tankgui
end

h=get(figs(1),'UserData');

ax=findobj(figs(1),'type','axes','Tag','picture');
if isempty(ax)
  sldemo_tankgui
  ax=findobj(figs(1),'type','axes','Tag','picture');
end

g=get(ax,'UserData');

%
% update water height patches
%
edge=sqrt(area);
outwidth = sqrt(outflow)/pi;
inwidth = sqrt(inflow)/pi;
len=max(1,max(edge,hilim)/5);
top=1.25*max(hilim, lolim+outwidth);
tubeheight=min(tankheight,lolim+outwidth);

ref.a = (edge/2-inwidth/2);
ref.b = ref.a+inwidth;
ref.c = ref.a-len/2;
ref.d = top+2*len+inwidth;
ref.e = lolim+outwidth;
in.vert.y=[top+len top+len top+2*len+inwidth top+2*len+inwidth top+len];
bottomedge = [-len -len -len -len -len];

set([g.water.main, g.water.outh, g.water.outv g.grey.outh ... 
      g.grey.outv g.grey.main], ...
   {'Ydata'},...
   {[0 0 tankheight tankheight 0];...
    [lolim lolim tubeheight tubeheight lolim];...
    [-len -len tubeheight tubeheight -len];...
    [tubeheight tubeheight ref.e ref.e tubeheight];...
    [tubeheight tubeheight ref.e ref.e  tubeheight];...
    [tankheight tankheight top top tankheight]})

on = {'on';'on';'off';'on';'on';'on';'on';'on'};
off = {'off';'off';'on';'off';'off';'off';'off';'off'};

%
% if in-valve is on....
%
if invalve==1;
  set([g.water.inv, g.valve.in.on, g.valve.in.off, g.drip.in], ...
      {'Visible'},on)
  set(g.water.inh,'Xdata',[-len ref.a ref.a -len -len])
  set(g.grey.inv,'Ydata',bottomedge)
  set(g.water.inv,'Ydata',in.vert.y)    
else
  set([g.water.inv, g.valve.in.on, g.valve.in.off, g.drip.in],{'Visible'},off)
  set(g.water.inh,'Xdata',[-len ref.c ref.c -len -len])
  set(g.grey.inv,'Ydata',in.vert.y)  
  set(g.water.inv,'Ydata',bottomedge)
end

%
% if out valve is on.....
%
if outvalve==1;
  set([g.water.outv,g.valve.out.on,g.valve.out.off,g.drip.out],{'Visible'},on)
  set(g.water.outh,'Xdata',[edge edge+len edge+len edge edge])
  set(g.grey.outh,'Xdata',[edge edge+len edge+len edge edge])    
else
  set([g.water.outv,g.valve.out.on,g.valve.out.off,g.drip.out],{'Visible'},off)
  set(g.water.outh,'Xdata',[edge edge+len/2 edge+len/2 edge edge])
  set(g.grey.outh,'Xdata',[edge edge+len/2 edge+len/2 edge edge])
  set(g.grey.outv,'Ydata',[-len -len ref.e ref.e -len])
end

%
% update comment text
%
if outvalve==1
  set(h.text(6),'String','Flusssh!');
elseif invalve==1
  set(h.text(6),'String','Filling');
else
  set(h.text(6),'String','Full');
end

%
% that's all folks
%
drawnow

% end mdlUpdate

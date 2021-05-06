% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function [data, abstractedmodel, X0, idOptions] = abstract(model,init_cond,input_range,cp_array,TotSimTime,sigData,hs,aristeo_options)
% ABSTRACT
%
% SYNOPSYS
%
%     [data, abstractedmodel] = abstract(model,init_cond,input_range,cp_array,TotSimTime,aristeo_options)
%
% DESCRIPTION
%
%     Simulate the input inputModel, starting at the provided initial conditions, and with the provided control input and computes an abstraction of the model.
%
% INPUTS :
%
%   - model : can be of type: function handle, an object of class hautomaton, or a simulink .mdl model.
%         In the last case, inputModel is of data type string.
%
%       * function handle :
%         It represents a pointer to a function which will be numerically
%         integrated using an ode solver (the default solver is ode45).
%         The solver can be changed through the option
%                   staliro_options.ode_solver
%         See documentation: <a href="matlab: doc staliro_options.ode_solver">staliro_options.ode_solver</a>
%
%       * Blackbox class object :
%         The user provides a function which returns the system behavior
%         based on given inputs and initial conditions. For example, this
%         option can be used when the system simulator is external to
%         Matlab. Please refer tp the staliro_blackbox help file.
%         See documentation: <a href="matlab: doc staliro_blackbox">staliro_blackbox</a>
%
%       * string :
%         It should be the name of the Simulink model to be simulated.
%
%       * hautomaton :
%         A hybrid automaton of the class hautomaton.
%         See documentation: <a href="matlab: doc hautomaton">hautomaton</a>
%
%       Examples:
%
%           % Providing directly a function that depends on state and time
%           model = @(t,x) [x(1) - x(2) + 0.1*t; ...
%                   x(2) * cos(2*pi*x(2)) - x(1)*sin(2*pi*x(1)) + 0.1 * t];
%
%           % Just an empty Blackbox object
%           model = staliro_blackbox;
%
%           % For other blackbox examples see the demos in demos folder:
%           staliro_demo_sa_simpleODE_param.m
%           staliro_demo_autotrans_02.m
%           staliro_demo_autotrans_03.m
%
%           % Simulink model under demos\SystemModelsAndData
%           model = 'SimpleODEwithInp';
%
%           % Hybrid automaton example (demo staliro_navbench_demo.m)
%           model = navbench_hautomaton(1,init,A);%

%       * Blackbox class object :
%         The user provides a function which returns the system behavior
%         based on given inputs and initial conditions. For example, this
%         option can be used when the system simulator is external to
%         Matlab. Please refer tp the staliro_blackbox help file.
%         See documentation: <a href="matlab: doc staliro_blackbox">staliro_blackbox</a>
%
%       * ss or dss :
%         A (descriptor) state-space model (see help file of ss or dss).
%         If the ss or dss models are discrete time models, then the
%         sampling time should match the sampling time for the input
%         signals (see staliro_options.SampTime). If they are not the same,
%         then an error will be issued.
%         See documentation: <a href="matlab: doc ss">ss</a>, <a href="matlab: doc dss">dss</a>, <a href="matlab: doc staliro_options.SampTime">staliro_options.SampTime</a>
%   - init_cond : a hyper-rectangle that holds the range of the initial
%       conditions (or more generally, constant parameters) and it should be a
%       Matlab n x 2 array, where
%			n is the size of the vector of initial conditions.
%		In the case of a Simulink model or a Blackbox model:
%			The array can be empty indicating no search over initial conditions
%			or constant parameters. For Simulink models in particular, an empty
%			array for initial conditions implies that the initial conditions in
%			the Simulink model will be used.
%
%       Format: [LowerBound_1 UpperBound_1; ...
%                          ...
%                LowerBound_n UpperBound_n];
%
%       Examples:
%        % A set of initial conditions for a 3D system
%        init_cond = [3 6; 7 8; 9 12];
%        % An empty set in case the initial conditions in the model should be
%        % used
%        init_cond = [];
%
%       Additional constraints on the initial condition search space can be defined
%       using the staliro option <a href="matlab: doc staliro_options.search_space_constrained">staliro_options.search_space_constrained</a>.
%
%   - input_range :
%       The constraints for the parameterization of the input signal space.
%       The following options are supported:
%
%          * an empty array : no input signals.
%              % Example when no input signals are present
%              input_range = [];
%
%          * a hyper-rectangle that holds the range of possible values for
%            the input signals. This is a Matlab m x 2 array, where m is the
%            number of inputs to the model. Format:
%               [LowerBound_1 UpperBound_1; ...
%                          ...
%                LowerBound_m UpperBound_m];
%            Examples:
%              % Example for two input signals (for example for a Simulink model
%              % with two input ports)
%              input_range = [5.6 7.8; 8 12];
%
%          * a cell vector. This is a more advanced option. Each input signal is
%            parameterized using a number of parameters. Each parameter can
%            range within a specific interval. The cell vector contains the
%            ranges of the parameters for each input signal. That is,
%                { [p_11_min p_11_max; ...; p_1n1_min p_1n1_max];
%                                    ...
%                  [p_m1_min p_m1_max; ...; p_1nm_min p_1nm_max]}
%            where m is the number of input signals and n1 ... nm is the number
%                  of parameters (control points) for each input signal.
%            Example:
%               See staliro_demo_constraint_input_signal_space_01.m
%       Additional constraints on the input signal search space can be defined
%       using the staliro option <a href="matlab: doc staliro_options.search_space_constrained">staliro_options.search_space_constrained</a>.
%            Example:
%               See staliro_demo_constraint_input_signal_space_01.m
%
%   - cp_array : contains the control points that parameterize each input
%       signal. It should be a vector (1 x m array) and its length must be equal
%       to the number of inputs to the system. Each element in the vector
%       indicates how many control points each signal will have.
%
%       Specific cases:
%
%       * If the signals generated using interpolation between the control
%         points, e.g., piece-wise linear or splines (for more options see
%         <a href="matlab: doc staliro_options.interpolationtype">staliro_options.interpolationtype</a>):
%
%         Initially, the control points are equally distributed over
%         the time duration of the simulation. The time coordinate of the
%         control points will remain constant unless the option
%
%					<a href="matlab: doc staliro_options.varying_cp_times">staliro_options.varying_cp_times</a>
%
%         is set (see the staliro_options help file for further instructions and
%         restrictions). The time coordinate of the first and last control
%         points always remains fixed.
%
%         Example:
%           cp_array = [1];
%               indicates 1 control point for only 1 input signal to the model.
%               One control point can only be used with piecewise constant
%               signals. If we assume that the total simulation time is 6 time
%               units and the input range is [0 2], then the input signal will
%               be:
%                  for all time t in [0,6] u(t) = const with const in [0,2]
%
%           cp_array = [4];
%               indicates 4 control points for only 1 input signal to the model.
%               If we assume that the total simulation time is 6 time units,
%               then the initial distribution of the control points will be:
%                            0   2   4   6
%
%           cp_array = [10 14];
%               indicates 10 control points for the 1st input signal and
%               14 for the second input.
%
%      * If the input_range is a cell vector, then the input range for each
%        control point variable is explicitly set. Therefore, we can
%        extract the number of control points from the number of
%        constraints. In this case, the cp_array should be set to emptyset.
%
%           cp_array = [];
%
%
%
%   - TotSimTime : total simulation time.
%
% - sigData: Return the input signals to the model. This is an array where the
%         first column contains the timestamps and the rest of the columns
%         the input signals. This is not used in staliro and it is provided
%         as an option for stand alone calls.
%
% - hs: Output trajectory is a struct with keys:
%               T      = col' vector of time stamps of simulated trajectory,
%               XT     = internal states at the time stamps,
%               YT     = outputs at the time stamps,
%               LT     = locations at the time stamps
%               locHis = sequence of visited locations
%               STraj  = XT or YT depending on staliro_options.spec_space,
%               GRD    = guard structure,
%               CLG = Control Location Graph.
%         Some of these will be empty, depending on the type of system being simulated. E.g. if it's a simulink model,
%         at least CLG, LT and GRD will be empty.
%
%   - aristeo_options : aristeo_options options. aristeo_options should be of type "aristeo_options".
%       If the default options are going to be used, then this input may be
%       omitted. For instructions on how to change Aristeo options,
%       see the aristeo_options help file for each desired property.
%

%
%     OUTPUTS
% data: an object of type iddata used to build the surrogate models
% abstractedmodel: the surrogate model
% X0: the initial conditions
% idOptions: the system identification set of options.

    %% checks that the parameters are correctly set
    global staliro_opt;
    global staliro_InputModel
    staliro_opt=aristeo_options;
    if(strcmp(aristeo_options.abstraction_algorithm,'arx')==0 && strcmp(aristeo_options.abstraction_algorithm,'armax')==0 && strcmp(aristeo_options.abstraction_algorithm,'bj')==0 && strcmp(aristeo_options.abstraction_algorithm,'tf')==0 && strcmp(aristeo_options.abstraction_algorithm,'ss')==0 && strcmp(aristeo_options.abstraction_algorithm,'nlarx')==0 && strcmp(aristeo_options.abstraction_algorithm,'hw')==0)
        error(strcat('The model structure ',aristeo_options.abstraction_algorithm,' is currently not supported by the abstraction procedure'));
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'arx') && (isequal(aristeo_options.na,-1) || isequal(aristeo_options.nb,-1) || isequal(aristeo_options.nk,-1)))
         error('If you are using the arx abstraction algoritm, the value of aristeo_options.na, aristeo_options.nb, aristeo_options.nc must be set');
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'armax') && (isequal(aristeo_options.na,-1) || isequal(aristeo_options.nb,-1) || isequal(aristeo_options.nc,-1) || isequal(aristeo_options.nk,-1)))
         error('If you are using the arx abstraction algoritm, the value of aristeo_options.na, aristeo_options.nb, aristeo_options.nc, aristeo_options.nk must be set');
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'bj') && (isequal(aristeo_options.nb,-1) || isequal(aristeo_options.nc,-1) || isequal(aristeo_options.nf,-1) || isequal(aristeo_options.nd,-1) || isequal(aristeo_options.nk,-1)))
         error('If you are using the arx abstraction algoritm, the value of aristeo_options.nb, aristeo_options.nc, aristeo_options.nf, aristeo_options.nd, aristeo_options.nk must be set');
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'tf') && (isequal(aristeo_options.np,-1) || isequal(aristeo_options.nz,-1)))
         error('If you are using the tf abstraction algoritm, the value of aristeo_options.np, aristeo_options.nz, must be set');
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'ss') && isequal(aristeo_options.nx,-1))
         error('If you are using the ss abstraction algoritm, the value of aristeo_options.nx must be set');
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'nlarx') && (isequal(aristeo_options.na,-1) || isequal(aristeo_options.nb,-1) || isequal(aristeo_options.nk,-1)))
         error('If you are using the ss abstraction algoritm, the value of aristeo_options.na, aristeo_options.nb, aristeo_options.nk must be set');
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'hw') && (isequal(aristeo_options.nb,-1) || isequal(aristeo_options.nf,-1) || isequal(aristeo_options.nk,-1)))
         error('If you are using the ss abstraction algoritm, the value of aristeo_options.nb, aristeo_options.nf, aristeo_options.nk  must be set');
    end

    % Determine the model type and perform error checks
    inputModelType = determine_model_type(model);
    % For back_box model types
    % for backward compatability of Blackbox option.
    % if user still using a Blackbox option then type cast model to Blackbox
    % class by copying the "inputModel" variable to "model_fcnptr" field
    if strcmp(inputModelType, 'function_handle')
        if staliro_opt.black_box == 1
            temp_model_ptr = model;
            model = staliro_blackbox();
            model.model_fcnptr = temp_model_ptr;
            staliro_InputModel = model;
        end
    end
    if isempty(sigData)|| isempty(hs)

    %% Generates the first input of the simulation

        if ~isempty(init_cond)

            InitialState= rand(size(init_cond,1),1).*(init_cond(:,2)-init_cond(:,1))+init_cond(:,1);

        else

            InitialState=[];

        end



        global simtime;
        UPoint=generateUPoint(cp_array,input_range);
        disp('Abstract');

        if(aristeo_options.dispinfo==1)

             disp(datestr(now));

        end

        simtimetic=tic;

        [hs, ~, sigData] = systemsimulator(model, InitialState, UPoint, TotSimTime, input_range, cp_array, aristeo_options);

        simtime=toc(simtimetic);
    end

    if(aristeo_options.dispinfo==1)
    disp('Resampling');
     disp(datestr(now));
    end
    v=ver('Matlab');
    if(isequal(v.Version,'7.11'))
        if(hs.T(2,1)>staliro_opt.SampTime)
            val= ceil(hs.T(2,1)/staliro_opt.SampTime);
            YT=hs.YT;
            YT=resample(YT,val,1);
        else
            val= ceil(staliro_opt.SampTime/hs.T(2,1));
            YT=hs.YT;
            YT=resample(YT,1,val);
        end
        %YT(:,preds(1).proj);
        if(sigData(2,1)>staliro_opt.SampTime)
            val= ceil(sigData(2,1)/staliro_opt.SampTime);
            U=sigData(:,2:1:size(sigData,2));
            U=resample(U,val,1);
        else
            val= ceil(staliro_opt.SampTime/sigData(2,1));
            U=sigData(:,2:1:size(sigData,2));
            U=resample(U,1,val);
        end
        minrow=min(size(YT,1),size(U,1));

        % If the developer specifies a specific portion of the trace to be used for system ideantification [aristeo_options.init_identification_time,aristeo_options.end_identification_time] and a sample step
        if aristeo_options.init_identification_time ~= 0 && aristeo_options.SampTime~=0 && aristeo_options.end_identification_time ~= 0
            data=iddata(YT((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),U((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),aristeo_options.SampTime);
        else
            data=iddata(YT(1:1:minrow,:),U(1:1:minrow,:),aristeo_options.SampTime);
        end
    else
        % when the model is hybrid, the output includes the output of the systme over time (YT) and the
        % current state (location) of the state machine over time (LT).
        % After the simulation, the location assumed by the state machine over time is in LT and the output in YT.
        % In order to evaluate the test inputs, staliro sets taliro_metric option with 'hybrid_inf' for this type of models.
        % Therefore, the abstraction procedure involves LT and YT for this type of
        % models
        if strcmp(aristeo_options.taliro_metric,'hybrid_inf')
            YT=[hs.YT hs.LT];
        else
            YT=hs.YT;
        end
        YT=resample(YT,hs.T,1/staliro_opt.SampTime);

        U=sigData(:,2:1:size(sigData,2));
        U=resample(U,sigData(:,1),1/staliro_opt.SampTime);

        minrow=min(size(YT,1),size(U,1));

        % If the developer specifies a specific portion of the trace to be used for system ideantification [aristeo_options.init_identification_time,aristeo_options.end_identification_time] and a sample step
        if aristeo_options.init_identification_time ~= 0 && aristeo_options.SampTime~=0 && aristeo_options.end_identification_time ~= 0
             data=iddata(YT((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),U((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),aristeo_options.SampTime);
        else
            data=iddata(YT(1:1:minrow,:),U(1:1:minrow,:),aristeo_options.SampTime);
        end
    end

    data = misdata(data);

    global absreftime;

    absreftimetic=tic;
    X0=[];

        disp('Abstracting');
         if(aristeo_options.dispinfo==1)
         disp(datestr(now));
        end

    try
    if(strcmp(aristeo_options.abstraction_algorithm,'arx'))

      v=ver('Matlab');
        if(isequal(v.Version,'7.11'))
          idOptions=[];
          abstractedmodel=arx(data,[aristeo_options.na aristeo_options.nb aristeo_options.nk],'Focus','simulation');
        else
          if(size(aristeo_options.na,1)==1)
              idOptions = arxOptions('EnforceStability',true,'InitialCondition','estimate','Focus','simulation');
          else
              idOptions = arxOptions('InitialCondition','estimate','Focus','simulation');
          end
          abstractedmodel=arx(data,[aristeo_options.na aristeo_options.nb aristeo_options.nk],idOptions);
        end
    end
    if(strcmp(aristeo_options.abstraction_algorithm,'armax'))
      idOptions = armaxOptions('EnforceStability',true,'InitialCondition','estimate','Focus','simulation');
      abstractedmodel=armax(data,[aristeo_options.na aristeo_options.nb aristeo_options.nc aristeo_options.nk],idOptions);
    end
    if(strcmp(aristeo_options.abstraction_algorithm,'bj'))
      try

        if(isequal(v.Version,'7.11'))
            idOptions = [];
             abstractedmodel=bj(data,[aristeo_options.nb aristeo_options.nc aristeo_options.nd aristeo_options.nf aristeo_options.nk],'InitialCondition','estimate','Focus','Stability');
        else
             idOptions = bjOptions('EnforceStability',true,'InitialCondition','estimate','Focus','simulation');
           abstractedmodel=bj(data,[aristeo_options.nb aristeo_options.nc aristeo_options.nd aristeo_options.nf aristeo_options.nk],idOptions);
         end
      catch exception
          idOptions = bjOptions('EnforceStability',true,'InitialCondition','estimate','Focus','simulation');
          abstractedmodel=NaN;
      end
    end
    if(strcmp(aristeo_options.abstraction_algorithm,'tf'))
      idOptions = tfestOptions('EnforceStability',true,'InitialCondition','estimate','Focus','simulation');
      abstractedmodel=tfest(data,aristeo_options.np,aristeo_options.nz,idOptions,'Ts',data.Ts);
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'ss'))
        if(isequal(v.Version,'7.11'))
          idOptions=[];
          abstractedmodel=pem(data,'nx',aristeo_options.nx,'Focus','stability');
        else
        idOptions = ssestOptions('EnforceStability',true,'InitialState','estimate');
        [abstractedmodel, X0]=ssest(data,aristeo_options.nx,idOptions);
        X0=[X0 X0];
        end
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'nlarx'))
       idOptions = nlarxOptions('Focus','simulation');
       idOptions.SearchOptions.MaxIter = 2;
       idOptions.Focus = 'simulation';
       idOptions.SearchMethod = 'lm';
       idOptions.SearchOptions.Advanced.LMStep=2;
              idOptions.SearchOptions.Advanced.MaxFunctionEvaluations=2;
      [abstractedmodel]=nlarx(data,[aristeo_options.na aristeo_options.nb aristeo_options.nk],wavenet('num',1),idOptions);
    end
    if (strcmp(aristeo_options.abstraction_algorithm,'hw'))
       idOptions = nlhwOptions();
       idOptions.SearchOptions.MaxIterations = 2;
        idOptions.SearchMethod = 'lm';
           idOptions.SearchOptions.Advanced.LMStep=2;
      [abstractedmodel]=nlhw(data,[aristeo_options.nb aristeo_options.nf aristeo_options.nk]);
    end
    absreftime=toc(absreftimetic);
    catch e
        absreftime=toc(absreftimetic);
        warning(e.message);
        abstractedmodel=[];

    end
end

% Copyright 2021 University of Luxembourg
%
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function [data,input,robustness]=check(model,cp_array,input_range,init_cond,X0,input,UPoint,phi,preds,data,TotSimTime,aristeo_options)
%CHECK Simulates the original model with the new falsification input
%signals returned by the falsification step.
% INPUTS:
%   - model : can be of type:
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
%       * ss or dss :
%         A (descriptor) state-space model (see help file of ss or dss).
%         If the ss or dss models are discrete time models, then the
%         sampling time should match the sampling time for the input
%         signals (see staliro_options.SampTime). If they are not the same,
%         then an error will be issued.
%         See documentation: <a href="matlab: doc ss">ss</a>, <a href="matlab: doc dss">dss</a>, <a href="matlab: doc staliro_options.SampTime">staliro_options.SampTime</a>
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
%           model = navbench_hautomaton(1,init,A);
%
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
%   - X0: A vector of initial conditions appropriate for the inputModel.
%
%   - UPoint:
%         A vector of control points. E.g. if the input signal is to be a linear interpolation between 3 input values,
%         UPoint contains the 3 values. If the system has more than one
%         input signal, then their control point values are concatenated in
%         UPoint: e.g. f inputSignal1 has 2 control points, and inputSignal2 has 3
%         control points, then UPoint will be a vector of the form
%         [cp1_1 cp1_2   cp2_1 cp2_2 cp2_3]
%         ------------   -----------------
%         Cntrl Pnts      Cntrl Pnts for
%         for inpSig1     inpSig2
%
%   - data: the data generated in previous iterations
%
%   - input: the input signals to check.
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
%   - phi : The formula to falsify. It should be a string. For the syntax of MTL
%       formulas type "help dp_taliro" (or see staliro_options.taliro for other
%       supported options depending on the temporal logic robustness toolbox
%       that you will be using).
%
%       Example:
%           phi = '!<>_[3.5,4.0] b)'
%
%       Note: phi can be empty in case the model is a hybrid automaton
%       object. In this case, an unsafe set must be provided in the hybrid
%       automaton.
%
%   - preds : contains the mapping of the atomic propositions in the formula to
%       predicates over the state space or the output space of the model. For
%       help defining predicate mappings type "help dp_taliro" (or see
%       staliro_options.taliro for other supported options depending on the
%       temporal logic robustness toolbox that you will be using).
%
%       In case of parameter mining:
%           If staliro is run for specification parameter mining, then set the
%           staliro option parameterEstimation to 1 (the default value is 0):
%               opt.parameterEstimation = 1;
%           and read the instructions under staliro_options.parameterEstimation
%           on how to define the mapping of the atomic propositions.
%
%   - TotSimTime : total simulation time.
%
%   - opt : s-taliro options. opt should be of type "staliro_options".
%       If the default options are going to be used, then this input may be
%       omitted. For instructions on how to change S-Taliro options,
%       see the staliro_options help file for each desired property.

    global aristeo_options_backup;
    global simtime;
    global RUNSTATS;
    RUNSTATS.new_run();
    LT=[];
    inputModelType = determine_model_type(model);
    % For back_box model types
    % for backward compatability of Blackbox option.
    % if user still using a Blackbox option then type cast model to Blackbox
    % class by copying the "model" variable to "model_fcnptr" field
    if strcmp(inputModelType, 'function_handle')
        if aristeo_options.black_box == 1
            temp_model_ptr = model;
            model = staliro_blackbox();
            model.model_fcnptr = temp_model_ptr;
            % update model type % For back_box model types
            inputModelType = determine_model_type(model);
        end
    end

    InitialState= X0;
    v=ver('Matlab');
    % the simulation depends on the model types
    if strcmp(inputModelType, 'simulink')
        simopt = simget(model);
        simopt = simset(simopt, 'InitialState', InitialState);

         if ~strcmp(aristeo_options.ode_solver, 'default')
            simopt = simset(simopt, 'Solver', aristeo_options.ode_solver);
        end

        load_system(model)
    end
    % exclude the case where the model has no input signals
    if (~isempty(input_range)) && (isempty(input) || sum(sum(isnan(input(:,:))))>0 || sum(sum(isinf(input(:,:))))>0 )
        disp('auto input generation')
        UPoint=generateUPoint(cp_array,input_range);
        opt=aristeo_options_backup;
        simtimetic=tic;
        [hs, ~, input] = systemsimulator(model, InitialState, UPoint, TotSimTime, input_range, cp_array, opt);
         simtime=toc(simtimetic);
         YT=hs.YT;
         T=hs.T;
    else
        % the simulation of the function handle is added
        if strcmp(inputModelType, 'function_handle')
            opt=aristeo_options_backup;
            % Choose ODE solver and simulate the system
            if strcmp(opt.ode_solver, 'default')
                ode_solver = 'ode45';
            else
                ode_solver = opt.ode_solver;
            end

            simtimetic=tic;
            %% Simulate system
            nb_ControlPoints=cp_array;
            for ii = 2:length(cp_array)
                nb_ControlPoints(ii) = nb_ControlPoints(ii)+nb_ControlPoints(ii-1);
            end
            [T, XT] = SimFcnPtr(ode_solver, model, TotSimTime, X0, UPoint, input_range, nb_ControlPoints, opt.interpolationtype, opt.varying_cp_times);
            simtime=toc(simtimetic);
            YT = XT; % Since there is no output signal explicitly defined
    % For back_box model types
        elseif strcmp(inputModelType, 'blackbox_model')
            steptime = (0:aristeo_options.SampTime:TotSimTime)';
             if ~isa(model.model_fcnptr,'function_handle')
                 error('S-Taliro: the balckbox model must be a function pointer')
             end

             [T, XT, YT, LT, CLG, GRD] = model.model_fcnptr(X0, TotSimTime, steptime, input(:,2:end));
             if (isempty(CLG) && isempty(GRD))
                 CLG = model.CLG ;
                 GRD = model.Grd ;
             else
                 if ~( isempty(model.CLG) && isempty(model.Grd))
                     if ~(isequal(CLG,model.CLG)&& isequal(GRD, model.Grd))
                         error('S-Taliro: Mismatch between staliro_blackbox and blackbox model in setting CLG/Grd. Please use either staliro_blackbox or blackbox model to set CLG and Grd')
                     end
                 end
             end
    elseif strcmp(inputModelType, 'hautomaton')

        error('Model not supported at the moment');
        % If you're advertising yourself as a hautomaton, then your h0 will contain the initial location and simulation time.
        % If multiple initial locations, then find the first that does not satisfy any switching conditions.
        % We always assume non-deterministic hybrid automata with ASAP transitions.
        % At this point, it is assumed that x0 is inside the invariant of at least one location. The simulator has its own error checking.
        if isempty(X0)
            if ~isfield(model.init,'h0') || isempty(model.init.h0)
                error(' systemsimulator: if a search space for the initial conditions is not provided to staliro, then the initial conditions must be specified in the field init.h0 of the hautomaton.')
            end
            if size(model.init.h0,1)==1
                X0 = model.init.h0';
            else
                X0 = model.init.h0;
            end
        end
        if length(model.init.loc)>1
            for cloc = 1:length(model.init.loc)
                l0 = cloc;
                actGuards = model.adjList{cloc};
                noag = length(actGuards);
                for jLoc = 1:noag
                    notbreak = 0;
                    nloc = actGuards(jLoc);
                    if staliro_opt.hasim_params(1) == 1
                        if isPointInSet(X0,model.guards(cloc,nloc).A,model.guards(cloc,nloc).b,'<')
                            notbreak = 1;
                            break
                        end
                    else
                        if isPointInSet(X0,model.guards(cloc,nloc).A,model.guards(cloc,nloc).b)
                            notbreak = 1;
                            break
                        end
                    end
                end
                if ~notbreak
                    break
                end
            end
        else
            l0 = model.init.loc(1);
        end
        if isempty(UPoint)
            h0 = [l0 0 XPoint'];
        else
            h0.h0 = [l0 0 XPoint'];
            h0.u = sigData;
        end
        if isfield(model, 'simulator')
            mysimulator = model.simulator;
            % A custom simulator's interface is assumed the same as the interface of systemsimulator
            [ht, ~, rc] = mysimulator(model, h0, UPoint, staliro_SimulationTime, staliro_InputBounds, nb_ControlPoints);
        else
            if strcmp(staliro_opt.ode_solver, 'default')
                ode_solver = 'ode45';
            else
                ode_solver = staliro_opt.ode_solver;
            end
            [ht, ~, rc] = hasimulator(model, h0, staliro_SimulationTime, ode_solver, staliro_opt.hasim_params);
        end
        T =  ht(:, 2); % get time
        XT = ht(:, 3:end); % get continuous state trajecotry
        LT = ht(:, 1); % get location trace
        YT = []; % no output signals (GF: change later to allow observations?)
        CLG = model.adjList; % location graph
        GRD = model.guards; % the transition guards
        else
            %% Similink models: Simulate system
            % The sim command is updated after the Matlab version 2010
            simtimetic=tic;
            if(isequal(v.Version,'7.11'))
                [T, ~, YT] = sim(model, [0 TotSimTime], simopt, input);
            else
                nb_ControlPoints=cp_array;
                for ii = 2:length(cp_array)
                    nb_ControlPoints(ii) = nb_ControlPoints(ii)+nb_ControlPoints(ii-1);
                end
                [hs, rc] = systemsimulator(model, X0, UPoint, TotSimTime, input_range, nb_ControlPoints);

                T = hs.T;
                XT = hs.XT;
                YT = hs.YT;
                LT = hs.LT;
                CLG = hs.CLG;
                GRD = hs.GRD;
            end
            simtime=toc(simtimetic);
        end
    end


    if(isequal(v.Version,'7.11'))
        if(input(2,1)>aristeo_options.SampTime)
            val= ceil(input(2,1)/aristeo_options.SampTime);
            input=input(:,2:1:size(input,2));
            input=resample(input,val,1);
        else
            val= ceil(aristeo_options.SampTime/input(2,1));
            input=input(:,2:1:size(input,2));
            input=resample(input,1,val);
        end
        robustness = dp_taliro(phi, preds, YT, T);

       if(YT(2,1)>aristeo_options.SampTime)
            val= ceil(YT(2,1)/aristeo_options.SampTime);
            YT=hs.YT;
            YT=resample(YT,val,1);
       else
           if(isequal(v.Version,'7.11'))
               val= ceil(aristeo_options.SampTime/T(2,1));
           else
                val= ceil(aristeo_options.SampTime/hs.T(2,1));
           end
            YT=hs.YT;
            YT=resample(YT,1,val);
       end
       minrow=min(size(YT,1),size(input,1));
       d=iddata(YT((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),input((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),aristeo_options.SampTime);
       data = merge(data,d);

    else

        robustness = dp_taliro(phi, preds, YT, T);
        input=resample(input(:,2:size(input,2)),input(:,1),1/aristeo_options.SampTime);
        if ((aristeo_options.fals_at_zero==0 && robustness>=0) || (aristeo_options.fals_at_zero==1 && robustness>0))
               YT=resample(YT,T,1/aristeo_options.SampTime);
               minrow=min(size(YT,1),size(input,1));
               d=iddata(YT((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),input((aristeo_options.init_identification_time/aristeo_options.SampTime+1):1:minrow,:),aristeo_options.SampTime);
               data = merge(data,d);
        end
    end
end

% the function wrapper is added to support the model type : function handle
function [T, XT] = SimFcnPtr(odesolver, fcn_ptr, simTime, Xpt, Upt, inpBound, tmpCP, IntType, CPType)

[T, XT] = feval(odesolver, @ourmodel_wrapper, [0 simTime], Xpt);

%% Nested functions
    function DX = ourmodel_wrapper(t, X)

        if isempty(Upt)

            DX = feval(fcn_ptr, t, X);

        else

            Uin = ComputeInputSignals(t, Upt, IntType, tmpCP, inpBound, simTime, CPType);
            if isempty(Uin)
                error('S-Taliro: the bounds of the input signals have been violated. The interpolation function does not respect the input signal bounds.')
            end
            DX = feval(fcn_ptr, t, X, Uin);

        end
    end
%% End of nested functions

end

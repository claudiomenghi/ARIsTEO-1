% Copyright 2021 University of Luxembourg
%
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

 function  [input,robustness,X0,UPoint]=falsify(abstractedmodel,init_cond, input_range, cp_array, phi, preds, TotSimTime, opt)
% FALSIFY performs the falsification of a given surrogate model
% INPUTS:
% - abstractedmodel: the surrogate model.
% - init_cond : a hyper-rectangle that holds the range of the initial
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
%   - opt : aristeo_options options. opt should be of type "aristeo_options".
%       If the default options are going to be used, then this input may be
%       omitted. For instructions on how to change Aristeo options,
%       see the aristeo_options help file for each desired property.
% OUTPUTS:
%   - input: the falsifying input signals.
%   - surrogaterobustness: the best robustness value returned by the surrogate model falsification.
%   - X0: the initial conditions values
%   - UPoint: A vector of control points. E.g. if the input signal is to be a linear interpolation between 3 input values,
%         UPoint contains the 3 values. If the system has more than one
%         input signal, then their control point values are concatenated in
%         UPoint: e.g. f inputSignal1 has 2 control points, and inputSignal2 has 3
%         control points, then UPoint will be a vector of the form
%         [cp1_1 cp1_2   cp2_1 cp2_2 cp2_3]
%         ------------   -----------------
%         Cntrl Pnts      Cntrl Pnts for
%         for inpSig1     inpSig2
%   - done: keeps true if the surrogate model exists.
%           Otherwise, done takes false.

    input =[];
    disp('Falsifying');
    if(opt.dispinfo==1)

     disp(datestr(now));
    end
    opt.black_box = 1;
    opt.runs=1;
    % aristeo-v2 changes
    % X0 and UPoint added
    X0=[];
    UPoint=[];
    if(isnan(abstractedmodel))
         warning('No faulty input found by the falsification procedure');
         input=[];
         robustness=100;
    else
        % encapsulate the surrogatemodel in a function handle
        sm=@sm_function;
        try
            % When you run SI on a blackbox model, and the spec predicate on X, the SI procedure consider the values of X as outputs.
            % Therefore, the surrogate model consideres the Xs as outputs.
            % Consequently, during the falsification it is necessary to consider Y as signals upon which the spec predicates.
            % a boolean variable indicating that the specification state of the MUT is X
            x_spec_state_flag=0;
            if strcmp(opt.spec_space,'X')
                opt.spec_space='Y';
                x_spec_state_flag=1;
            end
            [res, ~, ~] = staliro(sm,init_cond, input_range, cp_array, phi, preds, TotSimTime,opt);
            if x_spec_state_flag==1
                opt.spec_space='Y';
            end
            robustness=res.run.bestRob;
            %% Compute input signals - if any (this code is copied from systemsimulator function to remove sigData from params)
            % All this part is added
            % X0 and UPoint are returned, needed to run the simulation in Check.m
            if isempty(init_cond)
                UPoint=res.run.bestSample;
            else
                X0=res.run.bestSample(1:size(init_cond,1));
                UPoint=res.run.bestSample(size(init_cond,1)+1:end);
            end

            % temporal control points are used to compute the input
            % signals (as done in staliro)
            nb_ControlPoints=cp_array;
            for ii = 2:length(cp_array)
                nb_ControlPoints(ii) = nb_ControlPoints(ii)+nb_ControlPoints(ii-1);
            end
            steptime = (0:opt.SampTime:TotSimTime)';

            InpSignal = ComputeInputSignals(steptime, UPoint, opt.interpolationtype, nb_ControlPoints, input_range, TotSimTime, opt.varying_cp_times);
            % if inpSignal is empty then that means that the interpolation
            % function calculated a value outside of the predefined bounds
            if isempty(InpSignal)
                if opt.dispinfo == 1
                    warning('falsification issue:');
                    txtmsg =          'An input signal does not satisfy the input constraints and          ';
                    txtmsg = [txtmsg; 'the optimization function is not checking for this. You may want to '];
                    txtmsg = [txtmsg; 'use a different interpolating function.                             '];
                    disp(txtmsg)
                end
                input =[];
            else
                % input is returned
                input = [steptime InpSignal];
            end


        catch exception
            warning(exception.message);
            warning('No faulty input found by the falsification procedure');
            input=[];
            robustness=100;
        end
    end

end

% Copyright 2021 University of Luxembourg
%
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function UPoint = generateUPoint(cp_array,input_range)
%GENERATEUPOINT generate an array of the control points values
%associated to all the input signals
% INPUTS:
% - cp_array: contains the control points that parameterize each input
%       signal. It should be a vector (1 x m array) and its length must be equal
%       to the number of inputs to the system. Each element in the vector
%       indicates how many control points each signal will have.
% - input_range: The constraints for the parameterization of the input signal space.
% OUTPUTS:
% - UPoint: A vector of control points. E.g. if the input signal is to be a linear interpolation between 3 input values,
%         UPoint contains the 3 values. If the system has more than one
%         input signal, then their control point values are concatenated in
%         UPoint: e.g. f inputSignal1 has 2 control points, and inputSignal2 has 3
%         control points, then UPoint will be a vector of the form
%         [cp1_1 cp1_2   cp2_1 cp2_2 cp2_3]
%         ------------   -----------------
%         Cntrl Pnts      Cntrl Pnts for
%         for inpSig1     inpSig2

     UPoint=zeros(1,sum(cp_array));

    index=1;
    tmp_CP=zeros(size(cp_array,2));
    for n=1:1:size(cp_array,2)
        if n==1
         tmp_CP(n)=cp_array(n);
        else
         tmp_CP(n)=tmp_CP(n-1)+cp_array(n);
        end
        maxvalue=input_range(n,2);
        minvalue=input_range(n,1);
        diff=maxvalue-minvalue;

        for current_cp=1:1:cp_array(n)

            UPoint(index)=minvalue+rand(1)*diff;
            index=index+1;
        end
    end
end

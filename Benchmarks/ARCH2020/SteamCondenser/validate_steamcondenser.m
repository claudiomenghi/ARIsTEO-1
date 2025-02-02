% Copyright 2021 University of Luxembourg
 
% SPDX-FileCopyrightText: 2021 University of Luxembourg
% SPDX-License-Identifier: GPL-2.0-or-later
% Authors: see Authors.txt

function r = validate_steamcondenser(phi,rin,u)
    javaaddpath('../../../falstar/scala-library.jar');
    javaaddpath('../../../falstar/falstar.jar');

    in = {'Fs'};
    out = {'T', 'Fcw', 'Q', 'pressure'};

    T = u(end, 1);

    [tout, yout] = run_steamcondenser(u, T);

    y = [tout yout];
    r = falstar.Main.robustness(in, out, u, y, phi);

    disp(['accuracy: ' num2str(abs(rin - r))]);
    assert((r <  0) == (rin <  0));
    assert((r <= 0) == (rin <= 0));
end

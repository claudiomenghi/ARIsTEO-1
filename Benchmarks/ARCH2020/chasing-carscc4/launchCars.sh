#!/bin/bash
source /etc/profile

module load swenv/default-env/v1.1-20180716-production

module load base/MATLAB/2018a

matlab -nodisplay -nosplash -r "addpath(genpath('/home/users/cmenghi/ARCH/ARIsTEO/')); archcarscc4('$1','$2'); exit();"

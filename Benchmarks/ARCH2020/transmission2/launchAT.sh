#!/bin/bash
source /etc/profile

echo $1
echo $2

module load swenv/default-env/v1.1-20180716-production

module load base/MATLAB/2018a

matlab -nodisplay -nosplash -r "addpath(genpath('/home/users/cmenghi/ARCH/ARIsTEO/'));archat2('$1','$2'); exit();"

#!/bin/bash  
source /etc/profile

modelstructures=(arx)

for k in 2
do
    for j in "${modelstructures[@]}";
    do   
      sbatch -J "at$k$j" -n 10 -N 4 --priority=TOP -t 4-00:00:00  /home/users/cmenghi/ARCH/ARIsTEO/Benchmarks/ARCH2019/transmission2/launchAT.sh "$j" "$k" ;
    done
done

exit

#!/bin/bash  
source /etc/profile
rm *.out

modelstructures=(arx armax tf ss)

for k in 2 3 4 5 
do
    for j in "${modelstructures[@]}";
    do   
      sbatch -J "at$k$j" -n 10 -N 4 --priority=TOP -t 4-00:00:00  /home/users/cmenghi/ARCH/ARIsTEO/Benchmarks/ARCH2019/transmission/launchAT1.sh "$j" "$k" ;
    done
done

exit

#!/bin/bash  
source /etc/profile

#modelstructures=(arx armax tf ss)

modelstructures=(arx)

for k in 2 #3 4 5 
do
    for j in "${modelstructures[@]}";
    do   
      sbatch -J "cc$k$j" -n 10 -N 4 --priority=TOP -t 4-00:00:00  /home/users/cmenghi/ARCH/ARIsTEO/Benchmarks/ARCH2019/chasing-cars/launchCars.sh "$j" "$k" ;
    done
done

exit

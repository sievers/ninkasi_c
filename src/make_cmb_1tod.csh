#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
set DIR=/cita/d/raid-sievers/sievers/act/mapper/


#  ./mercator  -s -t 1e-6 -m 100 -i ${DIR}/src/simple_cmb_map_seed_1.dat -o equatorial_strip_badcommon_seed_${SEED}.map -d ${DIR}/simdata2/simple_equatorial_strip  -T equatorial_strip_badcommon_ -S ${SEED} -a -P -M  -B 

   ../bin/ninkasi <<EOF 
   -t 
1e-8
 -m 75 
-o ${DIR}/src/cmb_alldet.map 
  -d ${DIR}/real_data/cmb_test_alldet_  
-T cmb_alldet_   --precon   --pixsize=0.5   
-M
#-Q 
EOF
#  mpirun N  ../bin/ninkasi  -t 1e-8 -m 75 -o ${DIR}/src/mars_alldet.map   -d ${DIR}/real_data/mars_1det_  -T mars_alldet_   --precon   --pixsize=0.5  



#Iteration took  279.090 seconds.
#residual is    5.90199e+02 at iteration 1.
#Iteration took  280.499 seconds.
#residual is    4.04736e+01 at iteration 2.
#Iteration took  279.482 seconds.
#residual is    6.10053e+00 at iteration 3.
#Iteration took  279.031 seconds.
#residual is    2.72405e+00 at iteration 4.
#Iteration took  278.897 seconds.
#residual is    3.38232e+00 at iteration 5.
#Iteration took  280.262 seconds.
#residual is    7.56068e+00 at iteration 6.
#Iteration took  280.600 seconds.
#residual is    1.02430e+01 at iteration 7.
#Iteration took  279.719 seconds.
#residual is    5.37215e+00 at iteration 8.
#Iteration took  280.100 seconds.
#residual is    3.03332e+00 at iteration 9.
#Iteration took  279.527 seconds.
#residual is    1.29086e+00 at iteration 10.

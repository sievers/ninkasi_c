#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
set DIR=/cita/d/raid-sievers/sievers/act/mapper/


#  ./mercator  -s -t 1e-6 -m 100 -i ${DIR}/src/simple_cmb_map_seed_1.dat -o equatorial_strip_badcommon_seed_${SEED}.map -d ${DIR}/simdata2/simple_equatorial_strip  -T equatorial_strip_badcommon_ -S ${SEED} -a -P -M  -B 


  ../bin/ninkasi <<EOF
@sim @tol 1e-6
@maxiter 100 
@input ${DIR}/src/simple_cmb_map_seed_1.dat 
@output equatorial_weightmap_wnoise.map 
@data ${DIR}/simdata2/equatorial_hits_ 
@temp equatorial_weightmap_temp_  
#@precon
@blank
#@mean
@maxtod 2
@add_noise
#@quit
 
EOF

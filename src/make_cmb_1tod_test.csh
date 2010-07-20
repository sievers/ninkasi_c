#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
set DIR=/cita/d/raid-sievers/sievers/act/mapper/

./ninkasi << EOF
#@rows [ 5 6 : 16 11:12 13:13 1 2]

@tol 1e-4
@maxiter 100 
#@input ${DIR}/src/simple_cmb_map_seed_1.dat 
@output ${DIR}/src/mars_alldet.map  
#@data ${DIR}/real_data/mars_test_alldet_
@data  ${DIR}/real_data/cmb_test_alldet_
@temp mars_alldet
@pointing_offsets /cita/d/raid-sievers/sievers/act/ninkasi/from_ishmael/pointing_offset_mbac145.txt
@deglitch

@precon
@pixsize 0.5

#@use_cols [8 9 11 12 13 15 16 18 19 20 21 23]

#@use_rows [19] # - bad detector
use_rows [0:18 20:33]


#19,8 crappy


#@no_noise
#@detsize 0.69
#@detsize 0.5
#@blank
@mean
@maxtod 2
#@add_noise
#@quit
#finished 
EOF



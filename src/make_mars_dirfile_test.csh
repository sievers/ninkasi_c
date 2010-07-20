#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
#set DIR=/cita/d/raid-nolta/act/data/season1/
set DIR=/cita/d/raid-sievers/sievers/act/mapper/

./ninkasi << EOF

@tol 1e-4
@maxiter 100 
#@input ${DIR}/src/simple_cmb_map_seed_1.dat 
@output ./mars_alldet.map  
@data [${DIR}/real_data/1193729422.1193729537]
    
@temp mars_alldet_2map
@pointing_offsets /cita/d/raid-sievers/sievers/act/ninkasi/from_ishmael/pointing_offset_mbac145.txt
#@deglitch

@precon
@pixsize 0.15

#@use_cols [8 9 11 12 13 15 16 18 19 20 21 23]
#@use_rows [19] # - bad detector



@no_noise
#@blank
@mean
#@maxtod 2
#@add_noise
#@quit
#finished 
EOF



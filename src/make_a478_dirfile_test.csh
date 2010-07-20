#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
set DIR=/cita/d/raid-nolta/act/data/season1/

if ( 0 ) then 
  echo "running mpi"
  set exec="mpirun N"
  set TAG="mpi"
else
  echo "running serial"
  set exec=" "
  set TAG="serial"
endif


${exec} ./ninkasi << EOF

@tol 1e-6
@maxiter 100 
#@input ${DIR}/src/simple_cmb_map_seed_1.dat 
@output ./a478_alldet_${TAG}.map  
@data [${DIR}/1197424216.1197424241 ${DIR}/1197434953.1197434974]
#@data [${DIR}/1197434953.1197434974 ]
#@data [${DIR}/1197424216.1197424241 ${DIR}/1197424216.1197424241]
    
@rawmap blob_${TAG}.map

@temp a478_alldet_${TAG}
@pointing_offsets /cita/d/raid-sievers/sievers/act/ninkasi/from_ishmael/pointing_offset_mbac145.txt
@deglitch

@precon
@pixsize 0.5

#@use_cols [8 9 11 12 13 15 16 18 19 20 21 23]
#@use_rows [19] # - bad detector


@write_pointing


#@no_noise
#@mean
#@blank
#@rawonly
#@maxtod 2
#@add_noise
#@quit
#finished 
EOF



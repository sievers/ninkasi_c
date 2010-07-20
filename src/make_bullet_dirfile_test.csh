#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
set DIR=/cita/d/raid-nolta/act/data/season1/bullet/
set TAG=bullet
mpirun N  ./ninkasi << EOF

@tol 1e-4
@maxiter 100 
#@input ${DIR}/src/simple_cmb_map_seed_1.dat 
@output ./${TAG}_alldet.map  
@data [
${DIR}/1197254913.1197254962 
${DIR}/1197253693.1197253731
${DIR}/1197254303.1197254351]
#${DIR}/1197255523.1197255573
#${DIR}/1197256133.1197256204
#${DIR}/1197256743.1197256791
#${DIR}/1197257353.1197257401
#]
    
@rawmap ${TAG}_raw.map
    
@temp ${TAG}_alldet
@pointing_offsets /cita/d/raid-sievers/sievers/act/ninkasi/from_ishmael/pointing_offset_mbac145.txt
@deglitch
@mean

@precon
@pixsize 0.5


#@no_noise
#@blank



#@maxtod 2
#@add_noise
#@quit
#finished 
EOF



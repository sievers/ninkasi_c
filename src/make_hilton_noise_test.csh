#!/bin/csh
#set DIR=/Users/sievers/act/mapper/
#set DIR=/mnt/scratch-3week/sievers/act/mapper/
set DIR=/cita/d/raid-cita/sievers/act/hilton_sims/


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
@maxiter 5 
#@input ${DIR}/src/simple_cmb_map_seed_1.dat 
@output ./hilton_noise_test_alldet_nocommon_${TAG}.map  
#@data [${DIR}/1197424216.1197424241 ${DIR}/1197434953.1197434974]
@data [
${DIR}/1221390748.1221390761.ar1/
#${DIR}/1221993788.1221993903.ar1/
#${DIR}/1222252186.1222252283.ar1/
#${DIR}/1223630848.1223630964.ar1/
#${DIR}/1223700209.1234326370.ar1/
#${DIR}/1223717251.1234327922.ar1/
#${DIR}/1223975543.1223975689.ar1/
#${DIR}/1224130828.1224130940.ar1/
#${DIR}/1224217241.1224217412.ar1/
#${DIR}/1224475437.1224475573.ar1/
#${DIR}/1224578536.1224578653.ar1/
#${DIR}/1224923283.1224923419.ar1/
#${DIR}/1225009631.1225009755.ar1/
#${DIR}/1225078521.1225078640.ar1/
#${DIR}/1225267927.1225268037.ar1/
#${DIR}/1225354329.1225354462.ar1/
#${DIR}/1225699022.1234639633.ar1/
#${DIR}/1225768004.1234649611.ar1/
#${DIR}/1226112714.1226112873.ar1/
#${DIR}/1226285106.1226357005.ar1/
#${DIR}/1226302013.1226348217.ar1/
#${DIR}/1226560306.1226619322.ar1/
#${DIR}/1226646707.1226646820.ar1/
#${DIR}/1226991401.1226991516.ar1/
#${DIR}/1227491393.1227491515.ar1/
#${DIR}/1227939086.1227939191.ar1/
#${DIR}/1228025486.1228025589.ar1/
#${DIR}/1228283781.1228283887.ar1/
#${DIR}/1228783817.1228783833.ar1/
]



#@data [${DIR}/1197434953.1197434974 ]
#@data [${DIR}/1197424216.1197424241 ${DIR}/1197424216.1197424241]
    
@rawmap hilton_noise_nocommon_${TAG}.map

@temp hilton_noise_alldet__nocommon_${TAG}
@pointing_offsets /cita/d/raid-sievers/sievers/act/ninkasi/from_ishmael/pointing_offset_mbac145.txt
#@deglitch

@precon
@pixsize 0.5

#@use_cols [8 9 11 12 13 15 16 18 19 20 21 23]
#@use_rows [19] # - bad detector


#@write_pointing

@no_noise
#@mean
#@blank
#@rawonly
#@maxtod 2
#@add_noise
#@quit
#finished 
EOF



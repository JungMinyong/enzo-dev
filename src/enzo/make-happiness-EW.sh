#!/bin/sh  

# This makes enzo.exe on Happiness@SNU

echo "Making enzo on"
pwd
#make clean
make default
cd ../../
./configure
cd src/enzo/

# node 14
# export PATH=$PATH:/home/vinicius/install/openmpi-4.0.5/bin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vinicius/install/sevn_enzo/build_n14/lib64/sevn
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sykim/local/hdf5-1.12.1/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sykim/local/grackle/lib

# node 15
export PATH=$PATH:/home/vinicius/install/mpich-3.3.2/mpich3_15/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vinicius/install/sevn_enzo/build_n15/lib64/sevn
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vinicius/install/hdf5-1.12.1/hdf5_n15/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vinicius/install/grackle15/lib


# module add icc/latest
# module add mpi/latest
# export PATH=/appl/intel/oneapi/compiler/2021.4.0/linux/bin/intel64:$PATH
# export PATH=/appl/intel/oneapi/mpi/2021.4.0/bin:$PATH

make machine-linux-mpich-EW
#make precision-64 integers-32 particle-id-32 max-baryons-30 opt-aggressive lcaperf-no max-tasks-per-node-36 grackle-yes individualstar-yes nbody-yes photon-no cuda-no uuid-no fewbody-yes sevn-no
make io-64 precision-64 integers-32 particle-id-128 max-baryons-70 max-particle-attr-30 lcaperf-no max-tasks-per-node-64 grackle-yes photon-yes opt-aggressive cuda-no uuid-no individualstar-yes new-yield-tables-yes sevn-no memorypool-yes nbody-yes fewbody-yes
make show-config
make show-flags
make -j8
#cp enzo.exe enzo_spare.exe
#cp enzo.exe /data1/wispedia/enzo_nfw/lessSFB
echo "Make done!"
pwd
date

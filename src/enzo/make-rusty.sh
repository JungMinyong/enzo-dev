#!/bin/sh  

# This makes enzo.exe on Rusty@Flatiron

echo "Making enzo on"
pwd
make clean
make default
cd ../../
./configure
cd src/enzo/

#module add openmpi4
#module add modules/2.1.1-20230405
#module add intel-oneapi-compilers
#module add intel-oneapi-mpi
#module add modules/2.1-20220630
#module add modules/2.2-20230808 /2.1.1-20230405
#2.1.1-20230405
#module add gcc/10.3.0
#module add openmpi

module add modules/2.0-20220630 
#module add modules/2.2-20230808
module add gcc/7.5.0
#module add gcc/10.5.0 #/7.5.0
#module add openmpi/1.10.7
module add cuda/12.1.1
module add openmpi #/cuda-4.0.7
#module add openmpi
module add hdf5/1.8.22
#module add hdf5/1.10.11
#module add ucx
#module add cuda/12.1.1
#module add openmpi/1.10.7
#module add intel-parallel-studio

module show cuda
module show openmpi



make machine-rusty
#make io-64 precision-64 integers-32 particle-id-128 max-baryons-70 max-particle-attr-60 lcaperf-no max-tasks-per-node-64 grackle-yes photon-yes opt-aggressive cuda-no uuid-no individualstar-yes new-yield-tables-yes sevn-no memorypool-yes nbody-yes fewbody-no
make io-64 precision-64 integers-32 particle-id-128 max-baryons-70 max-particle-attr-60 lcaperf-no max-tasks-per-node-64 grackle-yes photon-yes opt-debug cuda-no uuid-no individualstar-yes new-yield-tables-yes sevn-no memorypool-yes nbody-yes fewbody-no
#new-yield-tables-yes   new-problem-types-yes
#make photon-yes grackle-yes opt-high io-64 max-baryons-66 max-particle-attr-46 new-yield-tables-yes individualstar-yes nbody-yes integers-32 uuid-no fewbody-yes sevn-no max-tasks-per-node-72
make show-config
make show-flags
##make -j3
make -j8
#cp enzo.exe enzo_spare.exe
#cp enzo.exe enzo_spare_v100.exe
#cp enzo.exe enzo_test_openmpi_escape.exe
#cp enzo.exe enzo_test_openmpi.exe
#cp enzo.exe enzo_sf_test.exe
cp enzo.exe enzo_debug.exe
#cp enzo.exe enzo_cosmo.exe
#cp enzo.exe enzo_orbit.exe
#cp enzo.exe enzo_nbn.exe
#cp enzo.exe enzo_test_irr.exe
#cp enzo.exe enzo_const.exe
#cp enzo.exe enzo_1e6.exe
echo "Make done!"
pwd
date

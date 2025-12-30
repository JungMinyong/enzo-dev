make default
cd ../../
./configure
cd src/enzo/

module purge
source ~/intel/oneapi/setvars.sh

export LD_LIBRARY_PATH=/home/vjl4366/pkg/hdf5-1.10.8/build/lib:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

ldconfig -p | grep hdf5

make machine-hap2
make io-64 precision-64 integers-32 particle-id-128 max-baryons-60 max-particle-attr-40 lcaperf-no max-tasks-per-node-64 grackle-yes photon-yes opt-aggressive cuda-yes
#make uuid-no 
make individualstar-yes new-yield-tables-yes memorypool-yes 
make nbody-yes fewbody-yes sevn-no
#make show-config
#make show-flags
make -j8


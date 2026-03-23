module load intel21/openmpi-5.0.0 
module load intel21/compiler-21
module load intel21/hdf5-1.10.5

make machine-linux-gnu
make io-64 precision-64 integers-32 particle-id-128 max-baryons-60 max-particle-attr-30 lcaperf-no max-tasks-per-node-64 grackle-yes photon-yes opt-aggressive cuda-no
make individualstar-yes new-yield-tables-yes memorypool-yes

make show-config
make show-flags
make -j8

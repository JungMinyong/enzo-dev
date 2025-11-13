make default
cd ../../
./configure
cd src/enzo/

export PATH=/home/wispedia/local/hdf5-1.12.1/bin:$PATH
#export PATH=/home/wispedia/local/gsl-2.7.1/bin:/home/wispedia/local/fftw-2.1.5/bin:$PATH
#export LD_LIBRARY_PATH=/home/wispedia/local/hdf5-1.12.1/lib:/home/wispedia/local/gsl-2.7.1/lib/gsl:/home/wispedia/local/fftw-2.1.5/lib:/home/wispedia/local/gsl-2.7.1/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/wispedia/local/hdf5-1.12.1/lib

make machine-icc
make io-64 precision-64 integers-32 particle-id-128 max-baryons-60 max-particle-attr-30 lcaperf-no max-tasks-per-node-64 grackle-yes photon-yes opt-aggressive cuda-no
#make uuid-no 
make individualstar-yes new-yield-tables-yes memorypool-yes 
#nbody-yes fewbody-yes sevn-no
make show-config
make show-flags
make -j8


#!/bin/bash

#module load nvhpc/23.1
#module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8

module load /leonardo/home/userexternal/ecalore0/nvidia/hpc_sdk/modulefiles/nvhpc/20.9

OPENMPI_LIB=/leonardo/home/userexternal/ecalore0/nvidia/hpc_sdk/Linux_x86_64/20.9/comm_libs/openmpi/openmpi-3.1.5

export MPILIB=$OPENMPI_LIB/lib
export MPIINC=$OPENMPI_LIB/include

# LOCAL LIB PATH
my_lib_path=$( echo "${LIBRARY_PATH}" | cut -d ':' -f 1 ) 

# CONFIGURE PARAMETERS SET FOR OPTIMIZED COMPILATION ON M100 (GPU = Nvidia Volta (V100) -> --ta=tesla,cc70) 
C_comp=pgcc 
C_comp_flags="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc80 -DUSE_MPI_CUDA_AWARE -I${MPIINC} -U__GNUC__" 
linker_flags="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc80 -lmpi -L${MPILIB}" #-L${my_lib_path} -U__GNUC__" 
CXX_comp=pgc++ 
CXX_comp_flags="-O3 -U__GNUC__"
cur_dir="--prefix=${PWD}" 

# LOAD PYTHON2.6
#module unload python
#module load profile/global
#module load pgi

# CONFIGURE + COMPILE
# if autotools are not built, build them
if [ ! -f ../configure ]; then cd ../; autoreconf -f -i; cd -; fi
# clean everything, configure, then make
if [ -f Makefile ];then make clean;fi
../configure CC="${C_comp}" CFLAGS="${C_comp_flags}" LDFLAGS="${linker_flags}" CXX="${CXX_comp}" CXXFLAGS="${CXX_comp_flags}" ${cur_dir}
make -j 32

# COPY EXECS IN CUR DIR
cp src/main .
cp tools/rgen .

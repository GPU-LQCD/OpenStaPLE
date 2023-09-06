#!/bin/bash

# LOCAL LIB PATH
my_lib_path=$( echo "${LIBRARY_PATH}" | cut -d ':' -f 1 ) 

# CONFIGURE PARAMETERS SET FOR OPTIMIZED COMPILATION ON LEONARDO
C_comp=pgcc 
C_comp_flags="-acc=noautopar,gpu -v -O0 -Minfo=all -gpu=cc80 -DUSE_MPI_CUDA_AWARE -I${OPENMPI_INCLUDE} -U__GNUC__"
linker_flags="-acc=noautopar,gpu -v -O0 -Minfo=all -gpu=cc80 -lmpi -L${OPENMPI_LIB} -L${my_lib_path} -U__GNUC__" 
CXX_comp=pgc++ 
CXX_comp_flags="-O3 -U__GNUC__"
cur_dir="--prefix=${PWD}" 

# LOAD PYTHON2.6
module unload python

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

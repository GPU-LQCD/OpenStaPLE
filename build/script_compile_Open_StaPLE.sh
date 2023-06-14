#!/bin/bash

# SIMULATION PARAMETERS
ACTION_TYPE='TLSM'         # must be either 'WILSON' or 'TLSM'

# LOCAL LIB PATH
my_lib_path=$( echo "${LIBRARY_PATH}" | cut -d ':' -f 1 ) 

# CONFIGURE PARAMETERS SET FOR OPTIMIZED COMPILATION ON M100 (GPU = Nvidia Volta (V100) -> --ta=tesla,cc70) 
C_comp=pgcc 
C_comp_flags="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc70 -DUSE_MPI_CUDA_AWARE -I${MPIINC}" 
linker_flags="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc70 -lmpi -L${MPILIB} -L${my_lib_path}" 
CXX_comp=pgc++ 
CXX_comp_flags=-O3 
cur_dir="--prefix=${PWD}" 

# LOAD PYTHON2.6
module unload python
module load profile/global
module load pgi

# CHECKS
if [ "${ACTION_TYPE}" != 'TLSM' ] && [ "${ACTION_TYPE}" != 'WILSON' ]; then
	echo "ERROR! ACTION_TYPE set to ${ACTION_TYPE} but must be either WILSON or TLSM"
	exit 1
fi


# COMMENT/UNCOMMENT DESIRED SIMULATION PARAMS
commented=$( grep "#define GAUGE_ACT_${ACTION_TYPE}" ../src/Include/common_defines.h | cut -d '#' -f 1 )
if test -z "${commented}"; then # desired action is already selected
	echo "Desired action ${ACTION_TYPE} already uncommented"
else # desired action is uncommented
	if [ "${ACTION_TYPE}" == 'TLSM' ]; then OTHER_ACTION='WILSON'; else OTHER_ACTION='TLSM'; fi
	echo "${ACTION_TYPE} will be uncommented; ${OTHER_ACTION} will be commented"
	sed -i "s://#define GAUGE_ACT_${ACTION_TYPE}:#define GAUGE_ACT_${ACTION_TYPE}:g" ../src/Include/common_defines.h
	sed -i "s:#define GAUGE_ACT_${OTHER_ACTION}://#define GAUGE_ACT_${OTHER_ACTION}:g" ../src/Include/common_defines.h
fi

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

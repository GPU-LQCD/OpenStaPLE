#!/bin/bash
# This is an example of a possible "configure" command for OpenStaPLE. 

MPFR_AND_GMP_LOCATION=/home/s.michele.mesiti/software

../configure CC="mpicc"\
 CFLAGS="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc70,cuda10.0 -DUSE_MPI_CUDA_AWARE"\
 LDFLAGS="-L$MPFR_AND_GMP_LOCATION/lib -lgmp"\
 CPPFLAGS="-I$MPFR_AND_GMP_LOCATION/include"\
 CXX="mpicxx"\
 CXXFLAGS="-O3"\
  --prefix=$(pwd)

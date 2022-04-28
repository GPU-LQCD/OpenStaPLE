#!/bin/bash

# LINES CONTAINING '#' ARE IGNORED
BLOCKS_FILE=dirac_blocks.txt

if test ! -f $BLOCKS_FILE 
then
    echo "ERROR: tile information file $BLOCKS_FILE does not exist."
    exit
fi

SRCDIR=${1-./src}

OBJDIR=$SRCDIR/OpenAcc

if test ! -d $OBJDIR
then 
    echo "ERROR: '.o' directory $OBJDIR does not exist."
    exit
fi

while read DDT0 DDT1 DDT2 DDG3
do
   echo 
   echo Compiling dirac test for 
   echo DEODOETILE0=$DDT0 DEODOETILE1=$DDT1  DEODOETILE2=$DDT2 DEODOEGANG3=$DDG3
   echo 
   rm $OBJDIR/fermion_matrix.o $OBJDIR/sp_fermion_matrix.o
   cd $SRCDIR
   make DEODOETILE0=$DDT0 DEODOETILE1=$DDT1  DEODOETILE2=$DDT2 DEODOEGANG3=$DDG3 deo_doe_test
   echo mv deo_doe_test $OLDPWD/deo_doe_test_$DDT0\_$DDT1\_$DDT2\_$DDG3
   mv deo_doe_test $OLDPWD/deo_doe_test_$DDT0\_$DDT1\_$DDT2\_$DDG3
   echo Compiled deo_doe_test_$DDT0\_$DDT1\_$DDT2\_$DDG3
   cd -

done < <(grep -v '#' $BLOCKS_FILE)

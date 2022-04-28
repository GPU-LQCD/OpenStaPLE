#!/bin/bash

# LINES CONTAINING '#' ARE IGNORED
BLOCKS_FILE=puregauge_blocks.txt

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


while read IST0 IST1 IST2 ISG3 ST0 ST1 ST2 SG3
do
   echo 
   echo Compiling pure gauge test for 
   echo IMPSTAPTILE0=$IST0  IMPSTAPTILE1=$IST1  IMPSTAPTILE2=$IST2  IMPSTAPGANG3=$ISG3
   echo STAPTILE0=$ST0  STAPTILE1=$ST1  STAPTILE2=$ST2  STAPGANG3=$SG3

   echo 
   rm $OBJDIR/rettangoli.o $OBJDIR/ipdot_gauge.o $OBJDIR/plaquettes.o $OBJDIR/su3_utilities.o
   rm $OBJDIR/sp_rettangoli.o $OBJDIR/sp_ipdot_gauge.o 
   rm $OBJDIR/sp_plaquettes.o $OBJDIR/sp_su3_utilities.o
   cd $SRCDIR
   make IMPSTAPTILE0=$IST0  IMPSTAPTILE1=$IST1  IMPSTAPTILE2=$IST2  IMPSTAPGANG3=$ISG3\
       STAPTILE0=$ST0  STAPTILE1=$ST1  STAPTILE2=$ST2  STAPGANG3=$SG3 main
   echo mv main $OLDPWD/pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   mv main $OLDPWD/pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   echo Compiled pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   cd -
done < <(grep -v '#' $BLOCKS_FILE)

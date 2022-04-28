#!/bin/bash

MDDIAGFILE=$1

ACCEPTANCE=$(grep ACCEPTED $MDDIAGFILE |  awk '{a+=$2;b++;}END{print a/b "+-" sqrt(a*(b-a)/b)/b,"("b" trajs)" }')
echo ACCEPTANCE: $ACCEPTANCE

STEPS=$(grep 'STEPS:' $MDDIAGFILE | tail -n 1 | awk '{print $4}')
GAUGESTEPS=$(grep 'STEPS:' $MDDIAGFILE | tail -n 1 | awk '{print $7}')
echo STEPS/GAUGESUBSTEPS: $STEPS / $GAUGESTEPS


NMEAS=$(grep DGFHN $MDDIAGFILE | awk '($3!=0){b++}END{print b}')
GFHN=$(grep GFHN $MDDIAGFILE  |grep -v DGFHN | awk '($3!=0){a+=$3;c+=$3*$3;b++}END{print a/b,"(~"sqrt(c/b-a*a/(b*b))")"}')
DGFHN=$(grep DGFHN $MDDIAGFILE   | awk '($3!=0){a+=$3;c+=$3*$3;b++}END{print a/b,"(~"sqrt(c/b-a*a/(b*b))")"}') 
echo Gauge Force Norm / Norm of the Gauge Force Difference        $GFHN  /  $DGFHN   '['$NMEAS measurements ']'
NMEAS=$(grep DFFHN $MDDIAGFILE | awk '($3!=0){b++}END{print b}')
FFHN=$(grep FFHN $MDDIAGFILE |grep -v DFFHN  | awk '($3!=0){a+=$3;c+=$3*$3;b++}END{print a/b,"(~"sqrt(c/b-a*a/(b*b))")"}')
DFFHN=$(grep DFFHN $MDDIAGFILE   | awk '($3!=0){a+=$3;c+=$3*$3;b++}END{print a/b,"(~"sqrt(c/b-a*a/(b*b))")"}') 
echo Fermion Force Norm / Norm of the Fermion Force Difference    $FFHN  /  $DFFHN '['$NMEAS measurements ']'
DELTA_ACTION=$(grep 'D ' $MDDIAGFILE | awk '{a+=$14;a2+=$14*$14}END{print "mean=",a/NR,"sigma=",sqrt(a2/NR - (a/NR)^2) }')
echo DELTA ACTION:  $DELTA_ACTION
echo Update time / MD setup time / MD time / Metro time: 
UPDATETIME=$(grep TOTUPDATETIME $MDDIAGFILE | awk '{a+=$NF;c+=$NF*$NF;b++}END{print a/b,"+-",sqrt(c/b-(a/b)^2)}')
MDSETUPTIME=$(grep MDSETUPTIME $MDDIAGFILE | awk '{a+=$2;c+=$2*$2;b++}END{print a/b,"+-",sqrt(c/b-(a/b)^2)}')
MDTIME=$(grep MDTIME $MDDIAGFILE | awk '{a+=$2;c+=$2*$2;b++}END{print a/b,"+-",sqrt(c/b-(a/b)^2)}')
METROTIME=$(grep METROTIME $MDDIAGFILE | awk '{a+=$2;c+=$2*$2;b++}END{print a/b,"+-",sqrt(c/b-(a/b)^2)}')
echo $UPDATETIME / $MDSETUPTIME / $MDTIME / $METROTIME
echo CG-M iterations '(averages)':
CGMTOT=$(grep 'CG-M-iters\[TOT]' $MDDIAGFILE | awk '{a+=$2}END{print a/NR}' )
CGMMD=$(grep 'CG-M-iters\[MD] ' $MDDIAGFILE | awk '{a+=$2}END{print a/NR}' ) 
CGMFI=$(grep 'CG-M-iters\[FI] ' $MDDIAGFILE | awk '{a+=$2}END{print a/NR}' ) 
CGMLI=$(grep 'CG-M-iters\[LI] ' $MDDIAGFILE | awk '{a+=$2}END{print a/NR}' ) 
echo Tot: $CGMTOT,MD: $CGMMD, FI:$CGMFI, LI: $CGMLI

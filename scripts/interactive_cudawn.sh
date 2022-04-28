MACHINE=$1
N=$2

if !  [ $N  ]
then
N=1
fi
echo Going on cudawn$MACHINE with $N slots...

echo bsub -q gpuday -n $N -m cudawn$MACHINE -Is /bin/bash -l
bsub -q gpuday -n $N -m cudawn$MACHINE -Is /bin/bash -l

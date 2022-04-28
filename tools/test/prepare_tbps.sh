
GEOMFILE="geom_defines.txt"
PREFIX="test"
PREPARESLURM=no
MODULES_TO_LOAD_CSV="auto"
NRESGPUS=
SLURMPARTITION=
MEMORY=12000

while (( "$#" ))
do
    case $1 in 
        -c|--compileinfo)
            GEOMFILE=$2
            shift
            ;;
        -p|--pref)
            PREFIX=$2
            shift;;
        -s|--slurm)
            PREPARESLURM=yes
            ;;
        -m|--modules)
            MODULES_TO_LOAD_CSV=$2
            shift
            ;;
        -a|--partition)
            SLURMPARTITION=$2
            shift
            ;;
        -n|-ngpus)
            NRESGPUS=$2 # 'none' is an available option
            shift
            ;;
    esac
    shift
done

echo "geomfile $GEOMFILE"
echo "prefix $PREFIX"
echo "PREPARESLURM = $PREPARESLURM"
echo "modules = $MODULES_TO_LOAD_CSV"
echo "partition = $SLURMPARTITION"
echo "ngpus = $NRESGPUS"

if test \
    $PREFIX != "test" -a \
    $PREFIX != "benchmark" -a \
    $PREFIX != "profiling"
then
    echo "usage: $0 -c compile_info_filename -p prefix -s -m modules,to,load"
    echo "       prefix = test OR benchmark OR profiling"
    exit
fi

if test $PREPARESLURM == "yes"
then
    if test -z $SLURMPARTITION
    then
        echo "$0 Error: Slurm Partition must be selected with -a or --partition."
        exit
    fi
else
    if [ ! -z "$NRESGPUS" ] && [ "$NRESGPUS" != "none" ] 
    then
        echo "$0 Error: incompatible options:"
        echo "   number of reqested gpus $NRESGPUS specified without -s flag."
        exit
    fi      
    if [ ! -z $SLURMPARTITION ]
    then
        echo "$0 Error: incompatible options:"
        echo "   slurm partition $SLURMPARTITION specified without -s flag."
        exit
    fi
fi


SCRIPTSDIR=$(dirname $BASH_SOURCE)
echo Setting SCRIPTSDIR to $SCRIPTSDIR

FILETEMPLATEDIR=$SCRIPTSDIR


if test ! -f $GEOMFILE
then
    echo File $GEOMFILE does not exist! 
    exit
fi

echo Reading compilation/geometry parameters from $GEOMFILE ...



#####################################
### reading basic geometry info #####
#####################################
L0=$( grep -E "^\s*LOC_N0\s+" $GEOMFILE | awk '{print $2}')
L1=$( grep -E "^\s*LOC_N1\s+" $GEOMFILE | awk '{print $2}')
L2=$( grep -E "^\s*LOC_N2\s+" $GEOMFILE | awk '{print $2}')
L3=$( grep -E "^\s*LOC_N3\s+" $GEOMFILE | awk '{print $2}')

TASKS=$( grep -E "^\s*NRANKS_D3\s+" $GEOMFILE | awk '{print $2}')

# size of the local gauge configuration in bytes.
# the '+4' is needed only on the direction which is cut (L3)
SIZE=$((576*L0*L1*L2*(L3+4)))

LABEL=$L0$L1$L2$L3\_$TASKS

#####################################
## preparing module load command ####
#####################################
if test $MODULES_TO_LOAD_CSV == "auto"
then 
    echo Reading modules to load in slurm scripts...
    MODULES_TO_LOAD=$(module list -t 2>&1 | tail -n+2)
    # creating a string with all the modules in one line
    for MODULE_TO_LOAD in $MODULES_TO_LOAD
    do
        MODULES_TO_LOAD_TEMP="$MODULES_TO_LOAD_TEMP $MODULE_TO_LOAD"
    done
    MODULES_TO_LOAD=$MODULES_TO_LOAD_TEMP
else
    MODULES_TO_LOAD=$(echo $MODULES_TO_LOAD_CSV | sed 's/,/ /g')
fi

if test -z "$MODULES_TO_LOAD" 
then
    echo "WARNING: No module loaded."
    MODULES_TO_LOAD_COMMAND=""
else
    MODULES_TO_LOAD_COMMAND="module load $MODULES_TO_LOAD"
fi

echo Modules to load:$MODULES_TO_LOAD....


#####################################
## preparing base file ##############
#####################################
# preparing base benchmark file
PGFILEBASE=$PREFIX.pg.DPSP.set.proto #PREFIX= test, benchmark, profiling
if test ! -f "$FILETEMPLATEDIR/$PGFILEBASE" -a \
    -f "$FILETEMPLATEDIR/fermion_parameters.set"
then
    echo Directory $SCRIPTSDIR does not exist or not contains benchmarks templates.
    exit
else
    echo Copying template benchmark input files from directory $SCRIPTSDIR
fi

for MODE in dp sp
do
    echo
    echo "MODE: $MODE"
	if test $MODE == sp
	then 
		SINGLEPRECISIONMD=1
	elif test $MODE == dp
	then 
		SINGLEPRECISIONMD=0
	fi

	OUTPROTOSETFILE=$(echo ${PGFILEBASE%.proto} |sed 's/DPSP/'$MODE'/;')
	echo $PGFILEBASE | sed 's/DPSP/'$MODE'/'
	echo "  Writing $OUTPROTOSETFILE"

    cat $FILETEMPLATEDIR/$PGFILEBASE |\
 		sed 's/SEDNRANKS/'$TASKS'/;s/SEDNX/'$L0'/;s/SEDNY/'$L1'/;s/SEDNZ/'$L2'/;
		    s/SEDNODEDIM/16/;s/SEDNT/'$((L3*TASKS))'/;s/DPSP/'$MODE'/;
		    s/SINGLEPRECISIONMD/'$SINGLEPRECISIONMD'/'  >  tmp_$OUTPROTOSETFILE

    sed 's/SEDREVTEST/0/'  tmp_$OUTPROTOSETFILE >  $OUTPROTOSETFILE

	echo "  Writing $PREFIX.main.$MODE.set"
	cat $FILETEMPLATEDIR/fermion_parameters.set $OUTPROTOSETFILE >\
		$PREFIX.main.$MODE.set


    if test $PREFIX == "test"
    then 
        OUTPROTOSETFILEREVT=$(echo ${PGFILEBASE%.proto} |sed 's/DPSP/'$MODE'.revt/;')
        sed 's/SEDREVTEST/1/'  tmp_$OUTPROTOSETFILE >  $OUTPROTOSETFILEREVT

        echo "  Writing $PREFIX.main.$MODE.revt.set"
        cat $FILETEMPLATEDIR/fermion_parameters.set $OUTPROTOSETFILEREVT >\
            $PREFIX.main.$MODE.revt.set
    fi
    rm tmp_$OUTPROTOSETFILE



done

echo
echo "PRECISION INDEPENDENT TESTS:"
# CREATING .set FILES FOR deo_doe_test and inverter_multishift_test
echo "  Writing $PREFIX.deo_doe_test.set"
cp $PREFIX.main.dp.set $PREFIX.deo_doe_test.dpsp.set
echo "  Writing $PREFIX.inverter_multishift_test.set"
cp $PREFIX.main.dp.set $PREFIX.inverter_multishift_test.dpsp.set


#####################################
## setting up directories and #######
## scripts                    #######
#####################################
EXECUTABLES="deo_doe_test inverter_multishift_test pg" # in all cases
if test $PREFIX == "test"
then                                # only with tests
    EXECUTABLES="main $EXECUTABLES" # adding 'main' to executables
fi

# executable 'pg' is actually 'main'
rm ./bin/pg
ln -s ./main ./bin/pg

echo
echo CREATING DIRECTORIES...

for EXECUTABLE in $EXECUTABLES
do
    echo "  Executable: $EXECUTABLE"
    if test "$EXECUTABLE" == 'main' -o "$EXECUTABLE" == 'pg'
    then 
        if test "$PREFIX" == "test"
        then
           MODES='sp dp sp.revt dp.revt'
           if test "$NRESGPUS" != "none"
           then 
               NRESGPUS=$TASKS
               echo "selecting minimum number of gpus per testing, $TASKS"
           fi
        else
           MODES='sp dp'
        fi
    else 
        MODES='dpsp'
    fi
    for MODE in $MODES
    do
        echo "    MODE: $MODE"

        INPUTFILENAME="$PREFIX.$EXECUTABLE.$MODE.set"

        if test "$PREFIX" == "profiling"
        then 
            PROFILINGSTR="nvprof --log-file \"$EXECUTABLE.%p"\"
        else
            PROFILINGSTR=""
        fi
        OUTPUTFILENAME="out.$EXECUTABLE.$MODE.txt"
        SLURMFILENAME="$PREFIX.$EXECUTABLE.$MODE.slurm"
        BASHFILENAME="$PREFIX.$EXECUTABLE.$MODE.sh"

        DIRNAME=${BASHFILENAME%.sh}
        echo "      mkdir -p  $DIRNAME"
        mkdir -p  $DIRNAME


        if test $EXECUTABLE != "pg"
        then
            cp $SCRIPTSDIR/ratapproxes/* ./$DIRNAME
        fi 

        if test $PREPARESLURM == yes
        then

            if test $NRESGPUS != "none"
            then
                SLURM_GPU_GRES='#SBATCH --gres=gpu:'$NRESGPUS
            else
                SLURM_GPU_GRES=''
            fi

        cat > $DIRNAME/$SLURMFILENAME << EOF
#!/bin/bash
#SBATCH --job-name=test.${EXECUTABLE}_$LABEL
#SBATCH --ntasks=$TASKS
#SBATCH --cpus-per-task=1
#SBATCH --error=test.${EXECUTABLE}.%J.err 
#SBATCH --output=test.${EXECUTABLE}.%J.out
$SLURM_GPU_GRES
#SBATCH --partition=$SLURMPARTITION
#SBATCH --mem-per-cpu=$MEMORY

module purge
$MODULES_TO_LOAD_COMMAND
export PGI_ACC_BUFFERSIZE=$SIZE

rm stop

srun --cpu_bind=v,sockets $PROFILINGSTR ../bin/$EXECUTABLE ./$INPUTFILENAME > $OUTPUTFILENAME

EOF
        fi

        cat > $DIRNAME/$BASHFILENAME << EOF
#!/bin/bash

mpirun -n $TASKS $PROFILINGSTR ../bin/$EXECUTABLE ./$INPUTFILENAME | tee $OUTPUTFILENAME

EOF
        echo "      cp $INPUTFILENAME ./$DIRNAME/"
        cp $INPUTFILENAME ./$DIRNAME/
        echo "      chmod +x $DIRNAME/$BASHFILENAME"
        chmod +x ./$DIRNAME/$BASHFILENAME

    done 


done

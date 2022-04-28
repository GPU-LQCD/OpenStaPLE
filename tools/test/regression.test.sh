#!/usr/bin/env bash
# This script runs the tests comparing the results between different commits.
# Options:
#   -s, --scheduler: The job scheduler that will be used to launch the tests.
#                    If 'bash' is chosen, everything will be run in serial on 
#                    the machine this scripts runs on.
#   -o, --oldcommit: The old commit to test against.
#   -n, --newcommit: The new commit to test.
#   -g, --geomfile : The file containing the geometry information for the 
#                    compilation.
#   -m, --modules  : Comma-separated-value list of modules to load (will be 
#                    propagated to the rest of the code.
#   -c, --config-options: Comma-separated-value list of options to pass to the
#                    configure-wrapper script for compilation.

# Default values
# Scheduler can be either 'bash' or a proper scheduling system. 
SCHEDULER=bash
V1COMMIT=$(git rev-parse HEAD)
V2COMMIT=
GEOMFILE=
GEOMFILE1=
GEOMFILE2=
MODULES_TO_LOAD_CSV=""
WORKDIR=$PWD
SCRIPTSDIR=$PWD/$(dirname $BASH_SOURCE)
REPODIR=$(dirname $(dirname $SCRIPTSDIR))
BUILDJOBS=4
SLURMPARTITION=
CONFIGOPTIONS_CSV=""
CONFIGOPTIONS_CSV1=""
CONFIGOPTIONS_CSV2=""


# Types of test
DOBUILDTESTONLY=no
DOMAINTEST=no
DOMAINREVTEST=no
DOPGTEST=no
DOPGREVTEST=no
DODIRACTEST=no
DOINVTEST=no

CURRENTCOMMIT=$(git rev-parse --abbrev-ref HEAD)


CLEANUP(){
    echo "Restoring repo state to $CURRENTCOMMIT..."
    git checkout -f $CURRENTCOMMIT
    exit
}

WAIT_FOR_RETURN(){
           echo "Next command: "
           echo $@
           echo "Press Return to continue"
           read
}

trap CLEANUP EXIT SIGINT
###################################################
####   READING COMMAND LINE OPTIONS   #############
###################################################
while (( "$#" ))
do
    case $1 in
    -s|--scheduler)
        SCHEDULER=$2
        shift
        ;;
    -v1|--version1) 
        V1COMMIT=$2
        shift
        ;;
    -v2|--version2)
        V2COMMIT=$2
        shift
        ;;
    -g|--geomfile)
        GEOMFILE=$2
        shift
        ;;
    -g1|--geomfile1)
        GEOMFILE1=$2
        shift
        ;;
    -g2|--geomfile2)
        GEOMFILE2=$2
        shift
        ;;
    -m|--modules)
        MODULES_TO_LOAD_CSV=$2
        shift
        ;;
    -c|--config-options)
        CONFIGOPTIONS_CSV=$2
        shift
        ;;
    -c1|--config-options1)
        CONFIGOPTIONS_CSV1=$2
        shift
        ;;
    -c2|--config-options2)
        CONFIGOPTIONS_CSV2=$2
        shift
        ;;
    -bj|--build-jobs)
        BUILDJOBS=$2
        shift
        ;;
    -t|--test)
        TESTSCSV=$2
        shift
        ;;
    -a|--partition) # -a is choosen for uniformity with prepare_tbps.sh
        SLURMPARTITION=$2
        shift
        ;;
    esac
    shift 
done

#########################################################
#######   LOGIC OPERATIONS ON OPTIONS   #################
#########################################################

# this function checks if an option has been set equal for both runs, e.g. GEOMFILE, or
# CONFIGOPTIONS_CSV, or MODULES_TO_LOAD_CSV. In that case sets the run-specific options 
# equal to the common value.
# If instead it has been set separately for the two runs (e.g. using GEOMFILE1 and GEOMFILE2)
# basically nothing is done, except some checking (you cannot specify GEOMFILE and GEOMFILE1 
# at the same time, for example).

OPTION12CHECK(){

    VARNAME=$1
    VARDEFAULT=$2
    CLOPTION=$3
    VARNAME1=${VARNAME}1
    VARNAME2=${VARNAME}2


    if [ ! -z ${!VARNAME1} ] || [ ! -z ${!VARNAME2} ] &&  [ ! -z ${!VARNAME} ]  
    then
        echo "$0 Error: Incompatible options: -$CLOPTION and -$CLOPTION\1 or -$CLOPTION\2"
        echo ${!VARNAME1}
        echo ${!VARNAME2}
        echo ${!VARNAME}
        exit
    fi

    if [ -z ${!VARNAME} ] && [ -z ${!VARNAME1} ] && [ -z ${!VARNAME2} ]
    then 
            eval $VARNAME=$VARDEFAULT
            if [ -z $VARDEFAULT ];then
                echo "No $VARNAME provided. Using default: (none)"
            else
                echo "No $VARNAME provided. Using default: \"${!VARNAME}\""
            fi
        fi

        if [ -z ${!VARNAME1} ] && [ ! -z ${!VARNAME} ] 
        then 
            eval $VARNAME1=${!VARNAME}
        fi
        if [ -z ${!VARNAME2} ] && [ ! -z ${!VARNAME} ] 
        then 
            eval $VARNAME2=${!VARNAME}
        fi

}

OPTION12CHECK GEOMFILE geom_defines.txt -g
OPTION12CHECK MODULES_TO_LOAD_CSV " " -m
OPTION12CHECK CONFIGOPTIONS_CSV "gcc" -c

# SLURM-related checks
if [ $SCHEDULER == "slurm" ]
then
    if [ -z $SLURMPARTITION ] 
    then
        echo "$0 Error: Slurm Partition must be selected with -a or --partition."
        exit
    fi
    SLURMFLAGS="-s -a $SLURMPARTITION"
else 
    SLURMFLAGS=""
fi

# Showing which options we got.
echo "Selecting scheduler \"$SCHEDULER\" (choose with -s)"
echo "Selecting commit \"$V1COMMIT\" for (1) (choose with -v1)"
echo "Selecting geom_defines file \"$GEOMFILE1\" for (1) (choose with -g1 or -g)" 
if [ ! -f "$GEOMFILE1" ]; then
    echo "Error: \"$GEOMFILE1\" does not exist. Exiting"
    exit 1
fi
echo "Selecting config wrapper options \"$CONFIGOPTIONS_CSV1\" for (1) (choose with -c1 or -c)" 
echo "Selecting modules \"$MODULES_TO_LOAD_CSV1\" for (1) (choose with -m1 or -m)" 

if [ ! -z "$V2COMMIT" ]
then
    echo "Selecting commit \"$V2COMMIT\" for (2) (choose with -v2)" 
    echo "Selecting geom_defines file \"$GEOMFILE2\" for (2) (choose with -g2 or -g)" 
    if [ ! -f "$GEOMFILE2" ]; then
        echo "Error: \"$GEOMFILE2\" does not exist. Exiting"
        exit 1
    fi
    echo "Selecting config wrapper options \"$CONFIGOPTIONS_CSV2\" for (2) (choose with -c2 or -c)" 
    echo "Selecting modules \"$MODULES_TO_LOAD_CSV2\" for (1) (choose with -m1 or -m)" 
else
    echo "Commit (2)  not selected (choose with -v2)" 
fi

# creating directory with all the test scripts, which will not be changed by git
# when we checkout to another commit.
echo "Using regression test scripts at commit \"$CURRENTCOMMIT\""

cp -r $SCRIPTSDIR test

TESTSLIST=(${TESTSCSV//,/ }) # using a bash array
for TEST in ${TESTSLIST[@]}
do
    case $TEST in
        build) 
            DOBUILDTESTONLY=yes
            echo "Build test will be performed for commit \"$V1COMMIT\""
            ;;
        main)
            if [ ! -z "$V2COMMIT" ]
            then 
                DOMAINTEST=yes
                DOBUILDTESTONLY=no # included automatically
                echo "main program regression test will be performed." 
	    else
                DOMAINTEST=no
		echo "ERROR: Main program regression test not possible."
	        echo "       Only one code version specified: select other commit with -v2"
		exit
	    fi
            ;;
        mainrev)
            DOMAINREVTEST=yes
            DOBUILDTESTONLY=no # included automatically
	        echo "reversibility test on main program will be performed on commit \"$V1COMMIT\"."
            ;;
        pg)
            if [ ! -z "$V2COMMIT" ]
            then 
                DOPGTEST=yes
                DOBUILDTESTONLY=no # included automatically
                echo "pure gauge molecular dynamics regression test will be performed." 
	    else
                DOPGTEST=no
                echo "ERROR: Pure gauge regression test not possible."
                echo "       Only one code version specified: select other commit with -v2"
		exit
	    fi
            ;;
        pgrev)
            DOPGREVTEST=yes
            DOBUILDTESTONLY=no # included automatically
            echo "reversibility test on pure gauge will be performed for commit $V1COMMIT."
            ;;
        dirac)
	    if [ ! -z "$V2COMMIT" ]
            then 
                DODIRACTEST=yes
                DOBUILDTESTONLY=no # included automatically
                echo "dirac operator regression test will be performed."
	    else
                DODIRACTEST=no
                echo "ERROR: Dirac regression test not possible."
                echo "       Only one code version specified: select other commit with -v2"
		exit
	    fi
            ;;
        inv|inverter)
            DOINVTEST=yes
            DOBUILDTESTONLY=no # included automatically
	        echo "inverter test will be performed on commit \"$V1COMMIT\"."
            ;;
        *)
            echo "ERROR: test \"$TEST\" not recognized."
            exit 1
            ;;
    esac
done

if [[ ${#TESTSLIST[@]} == 0 ]]
then 
    echo "No test has been selected to run."
    echo "Possible tests are: build,main,mainrev,pg,pgrev,dirac,inv (or inverter)."
    echo "Exiting."
    exit 1 
fi
echo "Press Return to continue..."
read




########################################################
########  COMPILING AND PREPARING RUNS   ###############
########################################################
# function for compilation and setting up of all directories
BUILD(){
    local COMMIT=$1
    local TEMPGEOMFILE=$2
    local TEMPCONFIGOPTIONS_CSV=$3 # compiler versions may change
    local TEMPMODULES_TO_LOAD_CSV=$4
 
    echo "Building at commit \"$COMMIT\"..."
    if [ ! -z $TEMPMODULES_TO_LOAD_CSV ]
    then
        echo 'Loading modules...'
        for MODULE in $(echo $TEMPMODULES_TO_LOAD_CSV | sed 's/,/ /g')
        do
            echo "module load $MODULE"
            module load $MODULE
        done
    fi

    cd $REPODIR
    autoreconf --install
    git checkout $COMMIT
    if [ $? -ne 0 ]
    then 
        echo "ERROR: could not checkout to commit $COMMIT"
        echo "Exiting."
        exit 1
    fi
    mkdir $WORKDIR/$COMMIT
    # configure_wrapper and configure need the geom_defines.txt file to be named
    # exactly geom_defines.txt
    cp $WORKDIR/$TEMPGEOMFILE $WORKDIR/$COMMIT/geom_defines.txt
    cd $WORKDIR/$COMMIT
    yes | $REPODIR/configure_wrapper $( echo $TEMPCONFIGOPTIONS_CSV | sed 's/,/ /g')
    make -j$BUILDJOBS &&  make install 
    if [ $? -ne 0 ]
    then 
        echo "Error: make failed, in $PWD"
        CLEANUP
    fi
}

PREPARE_RUN(){
    # variables that are expected to change form call to call are passed as arguments.
    local COMMIT=$1
    local TEMPGEOMFILE=$2
    local TEMPCONFIGOPTIONS_CSV=$3 # compiler versions may change
    local TEMPMODULES_TO_LOAD_CSV=$4
    # the other are taken as global.
    echo "Preparing run:"
    echo "Commit: \"$COMMIT\" "
    echo "compilation/geometry info file:  \"$TEMPGEOMFILE\""
    echo "Options for configure_wrapper:  \"$TEMPCONFIGOPTIONS_CSV\""
    echo "environment modules to load: \"$TEMPMODULES_TO_LOAD_CSV\""
    echo "Press return to continue..."
    read
 

    BUILD $COMMIT $TEMPGEOMFILE $TEMPCONFIGOPTIONS_CSV $TEMPMODULES_TO_LOAD_CSV

    if [ ! -z $TEMPMODULES_TO_LOAD_CSV ]
    then
        local TEMPMODULESFLAGS="-m $TEMPMODULES_TO_LOAD_CSV"
    else 
        local TEMPMODULESFLAGS=""
    fi

    $WORKDIR/test/prepare_tbps.sh -c geom_defines.txt -p test $SLURMFLAGS $TEMPMODULESFLAGS
    if [ $? -ne 0 ]
    then 
        echo "Error: prepare_tbps.sh failed, in $PWD"
        CLEANUP
    fi

    cd -
}

# compilation and setting up of all directories
if [ $DOBUILDTESTONLY == "yes" ] 
then 
    BUILD $V1COMMIT  $GEOMFILE1 $CONFIGOPTIONS_CSV1 $MODULES_TO_LOAD_CSV1
    exit 
else
    PREPARE_RUN $V1COMMIT $GEOMFILE1 $CONFIGOPTIONS_CSV1 $MODULES_TO_LOAD_CSV1
fi

if [ ! -z "$V2COMMIT" ]
then
      	PREPARE_RUN $V2COMMIT $GEOMFILE2 $CONFIGOPTIONS_CSV2 $MODULES_TO_LOAD_CSV2
fi


#####################################################
####### SETTING UP ALL THE TESTS ####################
#####################################################
# setting up all the tests, submitting all the jobs, linking the files
# that must be the same

echo "Runs have been prepared."

if [ $DOMAINTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    echo "Setting up main program test on both commits."
    echo "Press Return to continue..."
    read
    for DPSP in dp sp
    do 
        # 1. job on first commit
        cd $WORKDIR/$V1COMMIT/test.main.$DPSP
        if [  $SCHEDULER == "slurm" ]
        then 
            COMMAND="sbatch test.main.$DPSP.slurm"
            WAIT_FOR_RETURN $COMMAND
            CAJOB=$($COMMAND| cut -d' ' -f4)
        else
           echo $PWD
           COMMAND="bash ./test.main.$DPSP.sh"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
        fi
        # 2. "connection" job
        cd $WORKDIR
        cat > test.main.$DPSP\.connection.sh-slurm << EOF
#!/bin/bash
#SBATCH --job-name=main.$DPSP\.connection
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=main.$DPSP\.connection.%J.err 
#SBATCH --output=main.$DPSP\.connection.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in $WORKDIR/$V1COMMIT/test.main.$DPSP/*norndtest* 
do 
	echo "Linking file \$file ..."
	echo ln -s \$file $WORKDIR/$V2COMMIT/test.main.$DPSP/\$(basename \$file)
	ln -s \$file $WORKDIR/$V2COMMIT/test.main.$DPSP/\$(basename \$file)
done

EOF
        if [  $SCHEDULER == "slurm" ]
        then
            COMMAND="sbatch test.main.$DPSP.connection.sh-slurm --dependency=afterok:$CAJOB"
            WAIT_FOR_RETURN $COMMAND
            CONNJOB=$( $COMMAND | cut -d' ' -f4)
        else
            COMMAND="bash ./test.main.$DPSP.connection.sh-slurm"
            WAIT_FOR_RETURN $COMMAND
            $COMMAND
        fi
        # 3. job on second commit
        cd $WORKDIR/$V2COMMIT/test.main.$DPSP
        if [  $SCHEDULER == "slurm" ]
        then
            COMMAND="sbatch test.main.$DPSP.slurm --dependency=afterok:$CONNJOB"
            WAIT_FOR_RETURN $COMMAND
            CBJOB=$( $COMMAND | cut -d' ' -f4)
        else
            COMMAND="bash ./test.main.$DPSP.sh"
            WAIT_FOR_RETURN $COMMAND
            $COMMAND
        fi
        # 4. final-check job to compare the results
        cd $WORKDIR
        # WHAT EXACTLY NEED TO BE CHECKED? TO DEFINE
        cat > test.main.$DPSP\.finalcheck.sh-slurm << EOF
#!/bin/bash
#SBATCH --job-name=main.$DPSP\.finalcheck
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=main.$DPSP\.finalcheck.%J.err 
#SBATCH --output=main.$DPSP\.finalcheck.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in $WORKDIR/$V1COMMIT/test.main.$DPSP/global* 
do 
    echo "Checking differences for \$(basename \$file)" 
    ./test/diff_ascii_files.py \$file $WORKDIR/$V1COMMIT/test.main.$DPSP/\$(basename \$file)
done

EOF

       if [  $SCHEDULER == "slurm" ]
       then
           COMMAND="sbatch test.main.$DPSP.finalcheck.sh-slurm --dependency=afterok:$CBJOB"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
       else
           COMMAND="bash test.main.$DPSP.finalcheck.sh-slurm | tee test.main.$DPSP.check"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
       fi
    done
fi

if [ $DOMAINREVTEST == "yes" ] 
then
    # both single and double precision tests will be run on V1 commit.
    echo "Setting up main program reversibility test"
    echo "Commit: " $V1COMMIT
    echo "Press Return to continue..."
    read

    for DPSP in dp sp
    do
	cd $WORKDIR/$V1COMMIT/test.main.$DPSP.revt
        if [  $SCHEDULER == "slurm" ]
        then 
           COMMAND="sbatch test.main.$DPSP.revt.slurm"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
        else
           COMMAND="bash ./test.main.$DPSP.revt.sh"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
        fi
    done
fi

if [ $DOPGTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    echo "Setting up pure gauge molecular dynamics test on both commits."
    echo "Press Return to continue..."
    read
    for DPSP in dp sp
    do 
        # 1. job on first commit
        cd $WORKDIR/$V1COMMIT/test.pg.$DPSP
        if [  $SCHEDULER == "slurm" ]
        then 
            COMMAND="sbatch test.pg.$DPSP.slurm"
            WAIT_FOR_RETURN $COMMAND
            CAJOB=$($COMMAND | cut -d' ' -f4)
        else
           COMMAND="bash ./test.pg.$DPSP.sh"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
        fi
        # 2. "connection" job
        cd $WORKDIR
        cat > test.pg.$DPSP\.connection.sh-slurm << EOF
#!/bin/bash
#SBATCH --job-name=pg.$DPSP\.connection
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=pg.$DPSP\.connection.%J.err 
#SBATCH --output=pg.$DPSP\.connection.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in $WORKDIR/$V1COMMIT/test.pg.$DPSP/*norndtest* 
do 
	echo "Linking file \$file ..."
	echo ln -s \$file $WORKDIR/$V2COMMIT/test.pg.$DPSP/\$(basename \$file)
	ln -s \$file $WORKDIR/$V2COMMIT/test.pg.$DPSP/\$(basename \$file)
done

EOF
        if [  $SCHEDULER == "slurm" ]
        then
            COMMAND="sbatch test.pg.$DPSP.connection.sh-slurm --dependency=afterok:$CAJOB"
            WAIT_FOR_RETURN $COMMAND
            CONNJOB=$($COMMAND | cut -d' ' -f4)
        else
            COMMAND="bash ./test.pg.$DPSP.connection.sh-slurm"
            WAIT_FOR_RETURN $COMMAND
            $COMMAND
        fi
        # 3. job on second commit
        cd $WORKDIR/$V2COMMIT/test.pg.$DPSP
        if [  $SCHEDULER == "slurm" ]
        then
            COMMAND="sbatch test.pg.$DPSP.slurm --dependency=afterok:$CONNJOB"
            WAIT_FOR_RETURN $COMMAND
            CBJOB=$($COMMAND | cut -d' ' -f4)
        else
            COMMAND="bash ./test.pg.$DPSP.sh"
            WAIT_FOR_RETURN $COMMAND
            $COMMAND
        fi
        # 4. final-check job to compare the results
        cd $WORKDIR
        # WHAT EXACTLY NEED TO BE CHECKED? TO DEFINE
        cat > test.pg.$DPSP\.finalcheck.sh-slurm << EOF
#!/bin/bash
#SBATCH --job-name=pg.$DPSP\.finalcheck
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=pg.$DPSP\.finalcheck.%J.err 
#SBATCH --output=pg.$DPSP\.finalcheck.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in $WORKDIR/$V1COMMIT/test.pg.$DPSP/global* 
do 
    echo "Checking differences for file \$(basename \$file)"
    ./test/diff_ascii_files.py \$file $WORKDIR/$V1COMMIT/test.main.$DPSP/\$(basename \$file) \
    \$file $WORKDIR/$V1COMMIT/test.pg.$DPSP/\$(basename \$file)
done

EOF

       if [  $SCHEDULER == "slurm" ]
       then
           COMMAND="sbatch test.pg.$DPSP.finalcheck.sh-slurm --dependency=afterok:$CBJOB"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
       else
           COMMAND="bash test.pg.$DPSP.finalcheck.sh-slurm | tee test.pg.$DPSP.check"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
       fi
    done

fi

if [ $DOPGREVTEST == "yes" ] 
then
    # both single and double precision tests will be run on V1 commit.
    for DPSP in dp sp
    do
    	cd $WORKDIR/$V1COMMIT/test.pg.$DPSP.revt
        if [  $SCHEDULER == "slurm" ]
        then 
        COMMAND="sbatch test.pg.$DPSP.revt.slurm"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
        else
           COMMAND="bash ./test.pg.$DPSP.revt.sh"
           WAIT_FOR_RETURN $COMMAND
           $COMMAND
        fi
    done
fi

if [ $DODIRACTEST == "yes" ] 
then
    # 1. job on first commit
    cd $WORKDIR/$V1COMMIT/test.deo_doe_test.dpsp
    if [  $SCHEDULER == "slurm" ]
    then 
        COMMAND="sbatch test.deo_doe_test.dpsp.slurm"
        WAIT_FOR_RETURN $COMMAND
        CAJOB=$($COMMAND | cut -d' ' -f4)
    else
        COMMAND="bash ./test.deo_doe_test.dpsp.sh"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    fi
    # 2. "connection" job
    cd $WORKDIR
    cat > test.deo_doe_test.dpsp.connection.sh-slurm << EOF
#!/bin/bash
#SBATCH --job-name=deo_doe_test.dpsp.connection
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=deo_doe_test.dpsp.connection.%J.err 
#SBATCH --output=deo_doe_test.dpsp.connection.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in test.dp.gaugeconf_save test_fermion
do 
	echo "Linking file \$file ..."
    ln -s $WORKDIR/$V1COMMIT/test.deo_doe_test.dpsp/\$file \
    $WORKDIR/$V2COMMIT/test.deo_doe_test.dpsp/\$file  
done

EOF
    if [  $SCHEDULER == "slurm" ]
    then
        COMMAND="sbatch test.deo_doe_test.dpsp.connection.sh-slurm --dependency=afterok:$CAJOB"
        WAIT_FOR_RETURN $COMMAND
        CONNJOB=$($COMMAND | cut -d' ' -f4)
    else
        COMMAND="bash ./test.deo_doe_test.dpsp.connection.sh-slurm"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    fi
    # 3. job on second commit
    cd $WORKDIR/$V2COMMIT/test.deo_doe_test.dpsp
    if [  $SCHEDULER == "slurm" ]
    then
        COMMAND="sbatch test.deo_doe_test.dpsp.slurm --dependency=afterok:$CONNJOB"
        WAIT_FOR_RETURN $COMMAND
        CBJOB=$( $COMMAND | cut -d' ' -f4)
    else
        COMMAND="bash ./test.deo_doe_test.dpsp.sh"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    fi
    # 4. final-check job to compare the results
    cd $WORKDIR
    cat > test.deo_doe_test.dpsp.finalcheck.sh-slurm << EOF
#!/bin/bash
#SBATCH --job-name=deo_doe_test.dpsp.finalcheck
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=deo_doe_test.dpsp.finalcheck.%J.err 
#SBATCH --output=deo_doe_test.dpsp.finalcheck.%J.out
#SBATCH --partition=$SLURMPARTITION

FILES_TO_CHECK=( "test_fermion" \
    "test_fermion_result_doe2" \
    "test_fermion_result_deo2" \ 
    "test_fermion_result_fulldirac2" \
    "sp_test_fermion_result_doe2" \ 
    "sp_test_fermion_result_deo2" \
    "sp_test_fermion_result_fulldirac2" )


for file in \${FILES_TO_CHECK[@]}
do 
    echo "Checking differences for file \$file" 
    ./test/diff_ascii_files.py $WORKDIR/$V1COMMIT/test.deo_doe_test.dpsp/\$file \
      $WORKDIR/$V2COMMIT/test.deo_doe_test.dpsp/\$file 
done

EOF
    if [  $SCHEDULER == "slurm" ]
    then
        COMMAND="sbatch test.deo_doe_test.dpsp.finalcheck.sh-slurm --dependency=afterok:$CBJOB"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    else
        COMMAND="bash test.deo_doe_test.dpsp.finalcheck.sh-slurm | tee test.deo_doe_test.dpsp.check"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    fi
fi

if [ $DOINVTEST == "yes" ] 
then
    cd $WORKDIR/$V1COMMIT/test.inverter_multishift_test.dpsp
    if [  $SCHEDULER == "slurm" ]
    then
        COMMAND="sbatch test.inverter_multishift_test.dpsp.slurm"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    else
        COMMAND="bash test.inverter_multishift_test.dpsp.sh | tee test.inverter_multishift_test.dpsp.check"
        WAIT_FOR_RETURN $COMMAND
        $COMMAND
    fi
fi



git checkout $CURRENTCOMMIT

```
   ____                    _____ _        _____  _      ______ 
  / __ \                  / ____| |      |  __ \| |    |  ____|
 | |  | |_ __   ___ _ __ | (___ | |_ __ _| |__) | |    | |__   
 | |  | | '_ \ / _ \ '_ \ \___ \| __/ _` |  ___/| |    |  __|  
 | |__| | |_) |  __/ | | |____) | || (_| | |    | |____| |____ 
  \____/| .__/ \___|_| |_|_____/ \__\__,_|_|    |______|______|
        | |                                                    
        |_|                                                    

        OpenACC STAggered Parallel LatticeQCD Everywhere

By E. Calore, M. Mesiti, F. Negro, G. Silvi 
(and others, which can be added if they feel offended)
```
# An humorous intro to this code

## 0. Create the `configure` and `Makefile.in`s in the repo
This step is done with 
```
autoreconf --install
```
This will generate a number of files (among which the configure and the
Makefile.in files), starting from the Makefile.am in the various directories
and configure.ac.

## 1. Build - General ideas

### 1.0 Dependencies
You need 
- gmp and mpfr, for the rational approximation generation. You need to point 
  your `configure` to the include/lib locations of these libraries with `-L`
  and `-I` (you should not have to use `-l`).
- very likely the pgi compiler. For MPI, best thing is to use the mpicc and
  mpicxx (c++ is used only to compile the remez algorithm).

### 1.1 Build

At build time the geometry of the lattice must be known, and must be passed to 
the "configure" script (notice that there is also a "configure_wrapper" script, 
which is a python script, which makes easier to call "configure": more on that is written
below). At present, "configure" reads a file named "geom_defines.txt" 
which contains the lattice dimensions and the dimensions of the OpenAcc tiles.
A commented example of the file can be found in the 'doc' directory.

In case the system you are working on uses environment modules, you need to load them before 
invoking "configure".

Example:
```
mkdir build
cd build

<add the geom_defines.txt file>

../configure_wrapper pgi cc35 cuda8.0 # if you like the configure_wrapper script

make && make install  # you may want to use make -j32 to be faster.
```
### 1.2 Usage of the plain "configure" script

You must also pass to configure:
- The compiler : `CC=<the compiler you want to use>`
- prefix : The place where you want the code to be installed. Note that it defaults to 
  a system directory you're likely not to have access to, but even if you have access to it 
  you should not install this software in system directories.
- CFLAGS : compiler flags. Will be discussed later, depends on the compiler used.
- LDFLAGS : the linker flags (e.g. the `-L/path`)
- CPPFLAGS : the C-preprocessor flags (e.g. `-I/path`)

### 1.3 Shortcuts to configure

For some cases, the python script `configure_wrapper` may help.
This simple wrapper allows the user to call call configure without remembering all the
options. If you don't trust this wrapper, you can look into the code and look at the 
options you need.
It is called, e.g., in this way:
```
../configure_wrapper pgi cc35 cuda8.0
```
where you must change the options according to your goals. If invoked in ths way, it will 
produce the command
```
./configure CC=pgcc CFLAGS="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc35,cuda8.0 -DUSE_MPI_CUDA_AWARE -I${MPIINC}" LDFLAGS="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc35,cuda8.0 -lmpi -L${MPILIB}"  --prefix=$(pwd)
```
and you will be asked if you want to execute it or not. Note that, for simplicity, 
the install directory (`--prefix`) is set to the working directory (the programs will be
installed in `$PWD/bin`).  If you want to change that, just specify 
`--prefix=<your explicit installing directory>`

If you know what you are doing or the `configure_wrapper` script does not work
for any reason, you can have a look at `doc/configure.example.sh`. 
```
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
```
On some machines, you will have to specify manually the location of the `mpfr` 
and `gmp` libraries (that are needed only by the code that generates the 
rational approximations). Also, using the mpi compiler wrapper may or may not
be beneficial.

After configuring, you need to run 
```
make -jN # with N > 1 to make it faster, 
```
And
```
make install 
```
to finally install in the `--prefix` directory (see above).

## 2. Benchmarks

In order to run all the necessary benchmarks, it is advised to create an aptly named 
directory. In such directory, create a `geom_defines.txt` file (see an example in 'docs').
Then build the software, using the commands described in the previous section.
You can also produce all the necessary slurm scripts and setting files using the script
`prepare_tbps.sh` in `tools/test`. For example,
```
mkdir my_benchmark
cd my_benchmark
<create the geom_defines.txt file>
bash ./configure.sh # a wrapper like ../doc/configure.sh
make && make install 
../tools/test/prepare_tbps.sh -c geom_defines.txt -p benchmark -s -a gpu -n 1
```
This script prepares some benchmark/test cases with an associated slurm script,
reading many of the necessary information from the file passed with `-c`, in 
this case `geom_defines.txt` (*the generated script should be reviewed by a human*).
Notice that `tools/test/prepare_tpbs.sh` assumes the executables are in ./bin.

**NOTE**: the `tools/test/prepare_tpbs.sh` script has not been tested too much.
See the script for the meaning of the options and for other options.

## 3. Documentation
Do you want to know more? A 26-page latex document is waiting for you.
How to create it:
```
cd doc
pdflatex general_guide.tex
```

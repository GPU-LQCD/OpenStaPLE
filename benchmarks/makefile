
# DIMENSIONS OF THE LOCAL LATTICE, to be written in input files
N0:=$(shell grep LOC_N0 make_settings.set | awk '{print $$2}')
N1:=$(shell grep LOC_N1 make_settings.set | awk '{print $$2}')
N2:=$(shell grep LOC_N2 make_settings.set | awk '{print $$2}')
N3:=$(shell grep LOC_N3 make_settings.set | awk '{print $$2}')
# NUMBER OF RANKS to be written in input files
NR3:=$(shell grep NRANKS_D3 make_settings.set | awk '{print $$2}')
# NUMBER OF RANKS PER NODE to be written in input files
NODEDIM:=$(shell grep NODEDIM make_settings.set | awk '{print $$2}')

# COMPILER CHOICE: PGI / GNU
# ARCH CHOICE: CC35 / CC60 / AMD / HASWELL
COMPILER:=$(shell grep COMPILER make_settings.set | awk '{print $$2}')
ARCH:=$(shell grep ARCH make_settings.set | awk '{print $$2}')

all: run/main run/deo_doe_test run/inverter_multishift_test \
     run/benchmark_for_profiling.set run/benchmark_for_profiling.puregauge.set \
     run/benchmark.set run/benchmark.puregauge.set run/rgen make_settings.set \
     test_blocks_dirac_operator.sh test_blocks_puregauge.sh dirac_blocks.txt \
     puregauge_blocks.txt

make_settings.set:
	$(error The file make_settings.set must be present.)


geom_defines.txt: make_settings.set
	grep "LOC_N0\|LOC_N1\|LOC_N2\|LOC_N3\|NRANKS_D3\|DEODOE\|IMPSTAP\|STAP\|SIGMA\|COMPILER\|ARCH" make_settings.set > geom_defines.txt

compiler_settings_library.txt: $(ROOTDIR)/build/compiler_settings_library.txt
	cp $(ROOTDIR)/build/compiler_settings_library.txt ./


makefile: make_settings.set geom_defines.txt compiler_settings_library.txt
	$(warning 'REMEMBER TO SPECIFY ROOTDIR')
	cd $(ROOTDIR)/src && ./double_to_single_transformer.py autoMode && cd - && echo 'MAKE SETTINGS:'  &&  cat make_settings.set && \
	$(ROOTDIR)/build/generate_makefile.py $(ROOTDIR)/src/OpenAcc/main.c $(ROOTDIR)/src/tests_and_benchmarks/deo_doe_test.c $(ROOTDIR)/src/tests_and_benchmarks/inverter_multishift_test.c $(ROOTDIR)/src/tools/checkinput.c  $(ROOTDIR)/src/tools/confrestart.c  $(ROOTDIR)/src/tools/confconverter.c > ./makefile 

run:
	if [ ! -d run ] ; then mkdir run ; fi

rgen run/rgen : makefile run
	make COMMIT_HASH=$(COMMIT_HASH) rgen

main run/main: makefile run
	make COMMIT_HASH=$(COMMIT_HASH) main -j32

deo_doe_test run/deo_doe_test: makefile run
	make COMMIT_HASH=$(COMMIT_HASH) deo_doe_test -j32

inverter_multishift_test run/inverter_multishift_test: makefile run
	make COMMIT_HASH=$(COMMIT_HASH) inverter_multishift_test -j32

fermion_parameters.set :
	$(warning 'REMEMBER TO SPECIFY ROOTDIR and BENCHDIR')
	cp $(ROOTDIR)/$(BENCHDIR)/fermion_parameters.set ./

run/benchmark.set: run/benchmark.puregauge.set fermion_parameters.set run
	cat fermion_parameters.set > run/benchmark.set && cat run/benchmark.puregauge.set >> run/benchmark.set 

run/benchmark_for_profiling.set: run/benchmark_for_profiling.puregauge.set fermion_parameters.set run
	cat fermion_parameters.set > run/benchmark_for_profiling.set && cat run/benchmark_for_profiling.puregauge.set >> run/benchmark_for_profiling.set 

benchmark.set.proto : 
	cp $(ROOTDIR)/$(BENCHDIR)/benchmark.set.proto ./

benchmark_for_profiling.set.proto : 
	cp $(ROOTDIR)/$(BENCHDIR)/benchmark_for_profiling.set.proto ./

run/benchmark.puregauge.set: benchmark.set.proto run
	sed 's/SEDNRANKS/'$(NR3)'/' benchmark.set.proto | sed 's/SEDNX/'$(N0)'/' |\
	   	sed 's/SEDNY/'$(N1)'/' | sed 's/SEDNZ/'$(N2)'/' | sed 's/SEDNODEDIM/'$(NODEDIM)'/' |  sed 's/SEDNT/'$$(($(N3)*$(NR3)))'/' > run/benchmark.puregauge.set 


run/benchmark_for_profiling.puregauge.set: benchmark_for_profiling.set.proto run
	sed 's/SEDNRANKS/'$(NR3)'/' benchmark_for_profiling.set.proto | sed 's/SEDNX/'$(N0)'/' |\
	   	sed 's/SEDNY/'$(N1)'/' | sed 's/SEDNZ/'$(N2)'/' | sed 's/SEDNODEDIM/'$(NODEDIM)'/' |  sed 's/SEDNT/'$$(($(N3)*$(NR3)))'/' > run/benchmark_for_profiling.puregauge.set 


test_blocks_dirac_operator.sh:
	cp $(ROOTDIR)/$(BENCHDIR)/test_blocks_dirac_operator.sh ./ 

test_blocks_puregauge.sh:
	cp $(ROOTDIR)/$(BENCHDIR)/test_blocks_puregauge.sh ./

dirac_blocks.txt:
	cp $(ROOTDIR)/$(BENCHDIR)/dirac_blocks.txt.example ./dirac_blocks.txt

puregauge_blocks.txt:
	cp $(ROOTDIR)/$(BENCHDIR)/puregauge_blocks.txt.example ./puregauge_blocks.txt

test: testdeodoe testmultishift testmain testmainpg

testdeodoe: run/deo_doe_test run/benchmark.set
	cd run && mpirun -n $(NR3) ./deo_doe_test ./benchmark.set | tee -a  deo_doe_test_output.txt

testmultishift: run/inverter_multishift_test run/benchmark.set
	cd run && mpirun -n $(NR3) ./inverter_multishift_test ./benchmark.set | tee -a  inverter_multishift_test_output.txt

testmainpg: run/main run/benchmark.puregauge.set
	cd run && mpirun -n $(NR3) ./main ./benchmark.puregauge.set | tee -a  puregauge_output.txt

testmain: run/main run/benchmark.set run/rgen
	cd run && mpirun -n $(NR3) ./main ./benchmark.set | tee -a  full_md_output.txt

clean: 
	cd run && rm -f main rgen deo_doe_test inverter_multishift_test benchmark.set benchmark.puregauge.set 	&& cd ..

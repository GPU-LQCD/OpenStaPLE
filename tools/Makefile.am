
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = rgen rgenfo ratnormalizer eval_ratapprox

__common_sources = \
	./rationalapprox.h\
	./rationalapprox.cpp\
	../src/RationalApprox/Remez/alg_remez.h\
	../src/RationalApprox/Remez/bigfloat.h

rgen_SOURCES = ./rationalapprox_calc.cc\
	../src/RationalApprox/Remez/alg_remez.cc\
	$(__common_sources)

rgen_LDFLAGS = -lmpfr -lgmp

rgenfo_SOURCES = ./rationalapprox_calc_fixed_order.cc\
	../src/RationalApprox/Remez/alg_remez.cc\
	$(__common_sources)

rgenfo_LDFLAGS = -lmpfr -lgmp

ratnormalizer_SOURCES = normalize_ratapprox.cc\
	$(__common_sources)

eval_ratapprox_SOURCES = eval_ratapprox.cc \
	$(__common_sources)



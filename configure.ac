
AC_PREREQ([2.63])
AC_INIT([OpenStaPLE],[1.0],[michele.mesiti@gmail.com])


AM_PATH_PYTHON([2.6],\
               [AC_MSG_NOTICE([Python found.])],\
      [AC_MSG_ERROR([Python >2.6 is necessary to create single precision code.],1)])

AC_MSG_NOTICE([Moving back to $srcdir/src and creating/updating single precision code...])
cd $srcdir/src
./double_to_single_transformer.py autoMode silentMode
AC_MSG_NOTICE([Single prec code created/updated. Moving back to $OLDPWD])
cd -

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

AC_PROG_RANLIB # Requested by autoreconf on sunbird
AM_PROG_AR     # Requested by autoreconf on sunbird

AC_PROG_CC 
# CHECKS SUGGESTED BY AUTOSCAN
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P

# check for libraries
AC_CHECK_LIB([m], [sqrt],[],[AC_MSG_ERROR([math library not found.])])
AC_CHECK_LIB([mpi], [MPI_Init],[],[AC_MSG_ERROR([mpi library not found.])])

# Checks for header files.
AC_CHECK_HEADER([inttypes.h],[],[AC_MSG_ERROR([inttypes.h not found])])
AC_CHECK_HEADER([stdint.h],[],[AC_MSG_ERROR([stdint.h not found])])
AC_CHECK_HEADER([stdlib.h],[],[AC_MSG_ERROR([stdlib.h not found])])
AC_CHECK_HEADER([string.h],[],[AC_MSG_ERROR([string.h not found])])
AC_CHECK_HEADER([sys/time.h],[],[AC_MSG_ERROR([sys/time.h not found])])
AC_CHECK_HEADER([gmp.h],[],[AC_MSG_ERROR([gmp.h not found] )])
AC_CHECK_HEADER([mpfr.h],[],[AC_MSG_ERROR([mpfr.h not found])])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# THIS MACRO IS OBSOLETE BUT FOR SOME REASONS AUTORECONF NEEDS IT ON COKA.
AM_PROG_CC_C_O

# Checks for library functions.
AC_FUNC_FSEEKO
AC_CHECK_FUNCS([gettimeofday memset pow sqrt strstr])



COMPILE_INFO_FILENAME="geom_defines.txt"

# geometry - CHECKS CAN BE DEEPER 
if test -f $COMPILE_INFO_FILENAME
then 
    AC_MSG_NOTICE(["Reading $COMPILE_INFO_FILENAME"]) 
else :
    AC_MSG_ERROR([ "File $COMPILE_INFO_FILENAME must be present. Aborting." ],[1])
fi



LOC_N0=$( grep -E "^\s*LOC_N0\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
LOC_N1=$( grep -E "^\s*LOC_N1\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
LOC_N2=$( grep -E "^\s*LOC_N2\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
LOC_N3=$( grep -E "^\s*LOC_N3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')

# ranks in the 3rd direction
NRANKS_D3=$( grep -E "^\s*NRANKS_D3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')

if test $NRANKS_D3 -gt 1 
then 
AC_CHECK_HEADERS([mpi.h])
fi

#tile sizes

DEODOETILE0=$( grep -E "^\s*DEODOETILE0\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
DEODOETILE1=$( grep -E "^\s*DEODOETILE1\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
DEODOETILE2=$( grep -E "^\s*DEODOETILE2\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
DEODOEGANG3=$( grep -E "^\s*DEODOEGANG3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
IMPSTAPTILE2=$( grep -E "^\s*IMPSTAPTILE2\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
IMPSTAPTILE1=$( grep -E "^\s*IMPSTAPTILE1\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
IMPSTAPTILE0=$( grep -E "^\s*IMPSTAPTILE0\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
IMPSTAPGANG3=$( grep -E "^\s*IMPSTAPGANG3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
STAPTILE2=$( grep -E "^\s*STAPTILE2\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
STAPTILE1=$( grep -E "^\s*STAPTILE1\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
STAPTILE0=$( grep -E "^\s*STAPTILE0\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
STAPGANG3=$( grep -E "^\s*STAPGANG3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
SIGMATILE2=$( grep -E "^\s*SIGMATILE2\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
SIGMATILE1=$( grep -E "^\s*SIGMATILE1\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
SIGMATILE0=$( grep -E "^\s*SIGMATILE0\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
SIGMAGANG3=$( grep -E "^\s*SIGMAGANG3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')

COMMIT_HASH=$(git rev-parse HEAD)

AC_DEFINE_UNQUOTED(LOC_N0,      [[$LOC_N0]])
AC_DEFINE_UNQUOTED(LOC_N1,      [[$LOC_N1]])
AC_DEFINE_UNQUOTED(LOC_N2,      [[$LOC_N2]])
AC_DEFINE_UNQUOTED(LOC_N3,      [[$LOC_N3]])
AC_DEFINE_UNQUOTED(NRANKS_D3,   [[$NRANKS_D3]])
AC_DEFINE_UNQUOTED(DEODOETILE0, [[$DEODOETILE0]])
AC_DEFINE_UNQUOTED(DEODOETILE1, [[$DEODOETILE1]])
AC_DEFINE_UNQUOTED(DEODOETILE2, [[$DEODOETILE2]])
AC_DEFINE_UNQUOTED(DEODOEGANG3, [[$DEODOEGANG3]])
AC_DEFINE_UNQUOTED(IMPSTAPTILE2,[[$IMPSTAPTILE2]])
AC_DEFINE_UNQUOTED(IMPSTAPTILE1,[[$IMPSTAPTILE1]])
AC_DEFINE_UNQUOTED(IMPSTAPTILE0,[[$IMPSTAPTILE0]])
AC_DEFINE_UNQUOTED(IMPSTAPGANG3,[[$IMPSTAPGANG3]])
AC_DEFINE_UNQUOTED(STAPTILE2,   [[$STAPTILE2]])
AC_DEFINE_UNQUOTED(STAPTILE1,   [[$STAPTILE1]])
AC_DEFINE_UNQUOTED(STAPTILE0,   [[$STAPTILE0]])
AC_DEFINE_UNQUOTED(STAPGANG3,   [[$STAPGANG3]])
AC_DEFINE_UNQUOTED(SIGMATILE2,  [[$SIGMATILE2]])
AC_DEFINE_UNQUOTED(SIGMATILE1,  [[$SIGMATILE1]])
AC_DEFINE_UNQUOTED(SIGMATILE0,  [[$SIGMATILE0]])
AC_DEFINE_UNQUOTED(SIGMAGANG3,  [[$SIGMAGANG3]])
AC_DEFINE_UNQUOTED(COMMIT_HASH, [[$COMMIT_HASH]])

#AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 tools/Makefile
])



AC_OUTPUT

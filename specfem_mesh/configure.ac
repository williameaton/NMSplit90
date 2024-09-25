AC_INIT([nmsplit90], [1.0], [weaton@princeton.edu])

# Initialise automake
AM_INIT_AUTOMAKE([foreign])

# Set the f90 compiler 
AC_PROG_FC

# Determine some options for compilation 
# Use wither --enable-cuda or --disable-cuda
# Default is without cuda 
AC_ARG_ENABLE([cuda],
    [AS_HELP_STRING([--enable-cuda], [Enable CUDA support])],
    [with_cuda=$enableval], [with_cuda=no])
# Same for OpenMPI
AC_ARG_ENABLE([openmpi],
    [AS_HELP_STRING([--enable-openmpi], [Enable OpenMPI support])],
    [with_openmpi=$enableval], [with_openmpi=no])


if test "$with_cuda" = "yes"; then
    AC_MSG_NOTICE([CUDA support enabled])
    AM_CONDITIONAL([USE_CUDA], [true])
else
    AC_MSG_NOTICE([CUDA support disabled])
    AM_CONDITIONAL([USE_CUDA], [false])
fi

if test "$with_openmpi" = "yes"; then
    AC_MSG_NOTICE([OpenMPI support enabled])
    AC_MSG_NOTICE([Using mpifort compiler])
    FC=mpifort
    AM_CONDITIONAL([USE_OPENMPI], [true])
else
    AC_MSG_NOTICE([OpenMPI support disabled])
    AC_MSG_NOTICE([Using nvfortran compiler])
    AM_CONDITIONAL([USE_OPENMPI], [false])
    FC=nvfortran
fi




# Check for CUDA Fortran support
AC_CHECK_PROG([NVFORTRAN], [nvfortran], [yes], [no])
if test "$NVFORTRAN" = "no"; then
    AC_MSG_ERROR([nvfortran is required for CUDA Fortran support])
fi

# Check for CUDA Fortran support
AC_CHECK_PROG([MPIFORT], [mpifort], [yes], [no])
if test "$MPIFORT" = "no"; then
    AC_MSG_ERROR([mpifort is required for OpenMPI Fortran suppor. Have you connected the path?])
fi




# Tell autoconf to find the Makefile and populate the missing bits
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

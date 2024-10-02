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


# Check for cuda support
if test "$with_cuda" = "yes"; then
    AC_MSG_NOTICE([CUDA support enabled])
    AM_CONDITIONAL([USE_CUDA], [true])
else
    AC_MSG_NOTICE([CUDA support disabled])
    AM_CONDITIONAL([USE_CUDA], [false])
fi
# Check for MPI and change compiler based on this 
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


# NMODES
AC_ARG_VAR([nmodes], [Number of modes used in params.f90])
if test -z "$nmodes"; then
    n_modes=1  # Default value if --nmodes is not provided
else
    n_modes=$nmodes
fi
# Ensure n_modes is a valid integer
case "$n_modes" in
    ''|*[!0-9]*)
        AC_MSG_ERROR([Please provide an integer value for nmodes using --nmodes=VALUE])
        ;;
esac
AC_SUBST([N_MODES], [$n_modes])


# NSETS
AC_ARG_VAR([nsets], [Number of sets used in params.f90])
if test -z "$nsets"; then
    n_sets=8 
else
    n_sets=$nsets
fi
# Ensure n_sets is a valid integer
case "$n_sets" in
    ''|*[!0-9]*)
        AC_MSG_ERROR([Please provide an integer value for nmodes using --nsets=VALUE])
        ;;
esac
AC_SUBST([N_SETS], [$n_sets])



# For the KD tree we want to incorporate the coretran library: 
AC_ARG_WITH([coretranmod],
    [AS_HELP_STRING([--with-coretranmod], [Path to the directory containing coretran .mod files])],
    [coretranmod=$withval], [coretranmod=/scratch/gpfs/we3822/coretran/include])

AC_ARG_WITH([coretranlib],
    [AS_HELP_STRING([--with-coretranlib], [Path to the directory containing .so files])],
    [coretranlib=$withval], [coretranlib=/scratch/gpfs/we3822/coretran/lib])

# Ensure that both directories are valid
if test ! -d "$coretranmod"; then
    AC_MSG_ERROR([Module directory $coretranmod does not exist.])
else 
    AC_MSG_NOTICE([Coretran modules: $coretranmod])
    AC_SUBST([CORETRAN_MODDIR], [$coretranmod])

fi

if test ! -d "$coretranlib"; then
    AC_MSG_ERROR([Library directory $coretranlib does not exist.])
else 
    AC_MSG_NOTICE([Coretran library: $coretranlib])
    AC_SUBST([CORETRAN_LIBDIR], [$coretranlib])
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

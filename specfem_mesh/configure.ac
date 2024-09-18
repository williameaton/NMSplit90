AC_INIT([nmsplit90], [1.0], [weaton@princeton.edu])

# Initialise automake
AM_INIT_AUTOMAKE

# Set the f90 compiler 
AC_PROG_FC
FC=nvfortran

# Check for CUDA Fortran support
AC_CHECK_PROG([NVFORTRAN], [nvfortran], [yes], [no])
if test "$NVFORTRAN" = "no"; then
    AC_MSG_ERROR([nvfortran is required for CUDA Fortran support])
fi


# Tell autoconf to find the Makefile and populate the missing bits
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

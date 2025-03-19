AC_INIT([test_nmsplit90], [1.0], [weaton@princeton.edu])
AM_INIT_AUTOMAKE([subdir-objects foreign])

# Check for Fortran compiler
AC_PROG_FC
FC=nvfortran


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




# Generate Makefile
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

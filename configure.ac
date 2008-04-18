AC_INIT([src/flowClust.c])

# : ${R_HOME=`R RHOME`}
# if test -z "${R_HOME}"; then
#   echo "could not determine R_HOME"
#   exit 1
# fi

## BLAS.
AC_ARG_WITH([blas],
[AC_HELP_STRING([--with-blas],
                [use system gsl BLAS library (if available), or specify it @<:@no@:>@])],
[BLAS_LIBS=$with_blas],
[BLAS_LIBS=`-lgslcblas`])

# I did not find any BLAS libraries, (R is probably not well configured)
# if test "${BLAS_LIBS}" = ""; then
#   BLAS_LIBS="-lgslcblas"
# fi

AM_PATH_GSL(1.6,[LIBS="${LIBS} ${BLAS_LIBS} -lgsl"],[AC_MSG_ERROR([Cannot find Gnu Scientific Library >=1.6])])


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT()

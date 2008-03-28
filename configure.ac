AC_INIT([src/flowClust.c])

# ACX_BLAS([], [AC_MSG_ERROR([Did not found any BLAS libraries])])

AM_PATH_GSL(1.8,[LIBS="${LIBS} -lgsl"],[AC_MSG_ERROR([Cannot find Gnu Scientific Library >=1.8])])

ACX_BLAS([], [BLAS_LIBS="${BLAS_LIBS} -lgslcblas"])

# : ${R_HOME=`R RHOME`}
# if test -z "${R_HOME}"; then
# 	echo "could not determine R_HOME"
# 	exit 1
# fi
# BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
# 
# # I did not find any of the BLAS libraries, then use gsl
# if test "${BLAS_LIBS}" = ""; then
# 	LIBS="${LIBS} -lgslcblas"
# 	else
# 	LIBS="-lgsl $BLAS_LIBS"
# fi
LIBS="-lgsl $BLAS_LIBS"
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT()

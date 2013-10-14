AC_PREREQ(2.13)

AC_INIT([flowClust], 1.x, [gfinak@fhcrc.org], [flowClust])

# We want the following logic
# Per default, use the included gsl and set GSL_DIR, GSL_VERSTR and do compilation of gsl
# Otherwise, specify location of gsl using --with-gsl

 AC_ARG_ENABLE(bundled-gsl,
           [AS_HELP_STRING([--enable-bundled-gsl=yes/no],
                           [wether to use bundled gsl library.]
                           )
              ],
           [GSL_BUNDLED=$enableval],
	   [GSL_BUNDLED=yes])


AC_ARG_WITH(gsl,
           [AS_HELP_STRING([--with-gsl=DIR],
                           [root directory of GSL installation (defaults to /usr/local)])],
           [GSL_DIR="${with_gsl}"],
	   [GSL_DIR=yes])

if test "${GSL_BUNDLED}" = "yes" ; then
    AC_MSG_NOTICE([Using bundled GSL.])
    BUILD_DIR="`pwd`"
    AC_MSG_NOTICE([Searching for the compilers specified by R.])
    : ${R_HOME=`R RHOME`}
    if test -z "${R_HOME}"; then
      echo "could not determine R_HOME"
      exit 1
    fi
    CC=`"${R_HOME}/bin/R" CMD config CC`
    CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
    CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
    LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
    MAKE=`"${R_HOME}/bin/R" CMD config MAKE`
    AC_PROG_CC
    AC_PROG_CPP
    
    BASENAME="gsl-1.13"
	TGZNAME=trim-${BASENAME}.tar.gz
	cd src
	if test -d ./${BASENAME};
		then 
			echo 'found ' $BASENAME ' sources and tar archive;using what is there.'
		else
			echo "untarring gsl...";
			gunzip -dc ${TGZNAME} | tar xf -;
	fi;
	cd ../
    AC_MSG_NOTICE([Preparing to compile GSL.])
    cd src/${BASENAME}
    touch configure 
    BUILD_INSTALL_DIR="${BUILD_DIR}/src/gslbuild"
    rm -rf ${BUILD_INSTALL_DIR}
	
    AC_MSG_NOTICE([Configuring bundled GSL.])
     ./configure CC="${CC}" \
     CPP="${CPP}" \
     CFLAGS="${CFLAGS}" \
     CPPFLAGS="${CPPFLAGS}" \
     LDFLAGS="${LDFLAGS}" \
     CPICFLAGS="${CPICFLAGS}" \
     --prefix="${BUILD_INSTALL_DIR}" \
     --libdir="${BUILD_INSTALL_DIR}/lib" \
     --includedir="${BUILD_INSTALL_DIR}/include"\
     --enable-shared\
     --disable-static
      
    ${MAKE} install
    
    cd "${BUILD_DIR}"
    BUILD_INCLUDE_DIR="${BUILD_INSTALL_DIR}/include/gsl"
    #BUILD_LIB_DIR="${BUILD_INSTALL_DIR}/lib"
    #cp "${BUILD_LIB_DIR}"/* "${BUILD_DIR}"/src
    GSL_VERSION="1.13.0"
    PKG_CFLAGS="${PKG_CFLAGS} -I${BUILD_INCLUDE_DIR}"
    FLOWCLUST_DEFS="-DGSL_MAJOR=1 -DGSL_MINOR=13 -DGSL_STATIC"
else
    AC_MSG_NOTICE([Using external GSL.])
    

    # The logic is as follows:
    # Either the user specifically points us to gsl using --with-gsl='/something'
    # Or the intention to use external gsl is indicated by --with-gsl  
    if test "${GSL_DIR}" = "yes"; then
        dnl Since the user did not specify a directory, we first try pkg-config and then we will use /usr/local
        AC_MSG_NOTICE([No directory was specified for --with-gsl. Trying to find GSL using other methods.])
        if test -z "${PKG_CONFIG}" ; then
            AC_PATH_PROG(PKG_CONFIG, pkg-config)
        fi
        if ! test -z "${PKG_CONFIG}" ; then
            GSL_CONFIG="${PKG_CONFIG} libgsl"
            GSL_VERSION="`${GSL_CONFIG} --modversion || echo ''`"
            if test -z "${GSL_VERSION}" ; then
            	AC_MSG_NOTICE([pkg-config was not able to find the GSL library libgsl.])
                AC_MSG_NOTICE([Trying gsl-config."])
            	GSL_CONFIG="gsl-config"
            	GSL_VERSION="`${GSL_CONFIG} --version || echo ''`"
            	if test -z "${GSL_VERSION}" ; then
                	AC_MSG_NOTICE([gsl-config was not able to find the GSL library. This either indicates that GSL is old or that something is wrong. Verify GSL is installed and that PKG_CONFIG_PATH is correct.])
                	AC_MSG_NOTICE([Trying with /usr/local anyway."])
                	GSL_DIR="/usr/local"
                else
                	PKG_CFLAGS="`${GSL_CONFIG} --cflags`"
                	PKG_LIBS="`${GSL_CONFIG} --libs`"
                fi
            else
                PKG_CFLAGS="`${GSL_CONFIG} --cflags`"
                PKG_LIBS="`${GSL_CONFIG} --libs`"
            fi
        fi
    fi

    # At this point, either GSL_DIR is yes (using pkg-config/or gsl-config worked) or
    # it is set to something else, in which case we use dot to get the version number
    if ! test "${GSL_DIR}" = "yes"; then
        AC_MSG_NOTICE([Using GSL dir '${GSL_DIR}'])
		gsl_config="${GSL_DIR}/bin/gsl-config"
        if ! test -x "${gsl_config}" ; then
            AC_MSG_ERROR([$gsl_config not found.  Check GSL installation.])
            exit 1
        fi
        GSL_VERSION="`${gsl_config} --version`"
        PKG_CFLAGS="`${gsl_config} --cflags`"
        PKG_LIBS="`${gsl_config} --libs`"
    fi

    # At this point GSL_VERSION ought to exist, we do computations on it to infer MAJOR and MINOR version
    if test -z ${GSL_VERSION} ; then 
       AC_MSG_ERROR([Unable to determine GSL version. Report to flowClust maintainer.])
       exit 1
    fi
    AC_MSG_NOTICE([Found GSL version '${GSL_VERSION}'.])
    MAJOR=`echo ${GSL_VERSION} | cut -f1 -d"."`
    MINOR=`echo ${GSL_VERSION} | cut -f2 -d"."`
    if test -z ${MAJOR} || test -z ${MINOR} ; then
        AC_MSG_ERROR([Unable to infer GSL major and minor version, report to flowClust maintainer])
        exit 1
    fi
    AC_MSG_NOTICE([GSL major version is '${MAJOR}' and minor version is '${MINOR}'.])
    if (test ${MAJOR} -lt "1") || (test ${MAJOR} -eq "1" && test ${MINOR} -lt "6") ; then
        AC_MSG_ERROR([Need GSL >= 1.6])
        exit 1
    fi
    PKG_CFLAGS="${PKG_CFLAGS} -DGSL_MAJOR=${MAJOR} -DGSL_MINOR=${MINOR}"
fi

AC_MSG_NOTICE([Using the following compilation and linking flags for flowClust])
AC_MSG_NOTICE([   PKG_CFLAGS=${PKG_CFLAGS}])
AC_SUBST(PKG_CFLAGS)
AC_MSG_NOTICE([   PKG_LIBS=${PKG_LIBS}])
AC_SUBST(PKG_LIBS)
AC_MSG_NOTICE([   FLOWCLUST_DEFS=${FLOWCLUST_DEFS}])
AC_SUBST(FLOWCLUST_DEFS)
AC_MSG_NOTICE([Compiling using GSL version: '${GSL_VERSION}'.])
AC_MSG_NOTICE([Using bundled GSL: '${GSL_BUNDLED}'.])
AC_SUBST(GSL_VERSION)
AC_SUBST(GSL_BUNDLED)
AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
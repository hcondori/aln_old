dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(aln, 0.0)
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_SYSTEM

AC_PROG_CC

#AC_PROG_YACC
#AM_PROG_LEX
AM_INIT_AUTOMAKE()
LT_INIT

AC_CONFIG_FILES(Makefile \
	src/Makefile \
	src/common/Makefile
	src/structs/Makefile
	src/avx_aligners/Makefile \
	src/avx2_aligners/Makefile)
AC_OUTPUT


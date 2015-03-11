noinst_LTLIBRARIES = libavx2_aligners.la
libavx2_aligners_la_CFLAGS=-mavx2 -std=gnu11 -I$(top_srcdir)/src -I$(top_srcdir)/include
libavx2_aligners_la_LIBADD= $(top_builddir)/src/common/libcommon.la
libavx2_aligners_la_SOURCES = avx2_sw.h avx2_swi.h avx2_sw.c avx2_swi.c
	

noinst_LTLIBRARIES = libavx_aligners.la
libavx_aligners_la_CFLAGS=-mavx -std=gnu11 -I$(top_srcdir)/src -I$(top_srcdir)/include
libavx_aligners_la_LIBADD= $(top_builddir)/src/common/libcommon.la
libavx_aligners_la_SOURCES = avx_sw.h avx_sw.c

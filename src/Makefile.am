#BUILT_SOURCES=fasta_parser.h
#AM_YFLAGS = -d
bin_PROGRAMS=aln
aln_CFLAGS=-std=gnu11 -I$(top_srcdir)/include
#aln_SOURCES= alignment.h queue.h sw.h utils.h\
aln_SOURCES=aln.c
aln_LDADD=-lpthread \
	avx_aligners/libavx_aligners.la \
	avx2_aligners/libavx2_aligners.la \
	structs/libstructs.la

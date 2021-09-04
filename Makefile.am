# Makefile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libTeuk.la
bin_PROGRAMS = MAIN

libTeuk_la_SOURCES = myAlloc.c myFileIO.c myLog.c myOptparser.c myUtils.c \
etUtils.c etCore.c etExample.c mySpherical.c

libTeuk_la_LIBADD = -lgsl -lgslcblas -lgobject-2.0 -lglib-2.0 -lgthread-2.0 -lhdf5
libTeuk_la_LDFLAGS = -module -avoid-version

MAIN_SOURCES = main.c
MAIN_LDADD = libTeuk.la

all:$(bin_PROGRAMS)
$(bin_PROGRAMS): $(lib_LTLIBRARIES) 	

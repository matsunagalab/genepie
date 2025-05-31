#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8
if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana .

noinst_LTLIBRARIES = libkmeans_clustering.la
bin_PROGRAMS = kmeans_clustering$(EXEEXT)

libkmeans_clustering_la_SOURCES = \
        kc_option_str.fpp \
        kc_option.fpp \
        kc_control.fpp \
        kc_setup.fpp \
        kc_analyze.fpp

kmeans_clustering_SOURCES = kc_main.fpp

kc_main.o: kc_main.fpp $(noinst_LTLIBRARIES)

libkmeans_clustering_la_LIBADD = \
	$(top_srcdir)/src/lib/lib.la \
	$(top_srcdir)/src/analysis/libana/libana.la \
	$(INT_LAPACK_LIB)

kmeans_clustering_LDADD = \
	libkmeans_clustering.la

AM_FCFLAGS = \
	-I$(top_srcdir)/src/lib \
	-I$(top_srcdir)/src/analysis/libana

EXTRA_DIST = \
	Makefile.depends

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) *.f90 *~ *.mod

.fpp.f90:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
endif

.f90.lo:
if USEKCOMP
	mv $*.cpp.f90 $*.f90
endif
	$(LIBTOOL) --mode=compile $(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		$(AM_FCFLAGS) $(FCFLAGS) -c $*.f90

.f90.o:
if USEKCOMP
	mv $*.cpp.f90 $*.f90
endif
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		$(AM_FCFLAGS) $(FCFLAGS) -c $*.f90

depend: clean_depend
	python $(top_srcdir)/fortdep.py --objext lo *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)

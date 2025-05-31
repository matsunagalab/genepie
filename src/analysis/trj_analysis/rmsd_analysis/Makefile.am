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

noinst_LTLIBRARIES = librmsd_analysis.la
bin_PROGRAMS = rmsd_analysis$(EXEEXT)

librmsd_analysis_la_SOURCES = \
	ra_option_str.fpp \
	ra_option.fpp \
	ra_control.fpp \
	ra_setup.fpp \
	ra_analyze.fpp

rmsd_analysis_SOURCES = ra_main.fpp

ra_main.o: ra_main.fpp $(noinst_LTLIBRARIES)

librmsd_analysis_la_LIBADD = \
	$(top_srcdir)/src/lib/lib.la \
	$(top_srcdir)/src/analysis/libana/libana.la \
	$(INT_LAPACK_LIB)

rmsd_analysis_LDADD = \
	librmsd_analysis.la

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

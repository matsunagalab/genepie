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

noinst_LTLIBRARIES = librg_analysis.la
bin_PROGRAMS = rg_analysis$(EXEEXT)

librg_analysis_la_SOURCES = \
	rg_option_str.fpp \
	rg_option.fpp \
	rg_control.fpp \
	rg_setup.fpp \
	rg_analyze.fpp

rg_analysis_SOURCES = rg_main.fpp

rg_main.o: rg_main.fpp $(noinst_LTLIBRARIES)

librg_analysis_la_LIBADD = \
	$(top_srcdir)/src/lib/lib.la \
	$(top_srcdir)/src/analysis/libana/libana.la \
	$(INT_LAPACK_LIB)

rg_analysis_LDADD = \
	librg_analysis.la

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

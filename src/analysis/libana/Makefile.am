#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS := ../../lib .

lib_LTLIBRARIES = libana.la

install-libLTLIBRARIES:
	@# Not install libana.la

libana_la_SOURCES = \
	input_str.fpp \
	output_str.fpp \
	trajectory_str.fpp \
	measure.fpp \
	fileio_trj.fpp \
	fileio_atp.fpp \
	fileio_rtp.fpp \
	fileio_mt.fpp \
	input.fpp \
	output.fpp \
	trajectory.fpp \
	pbc_correct.fpp \
	analysis_timers.fpp

EXTRA_DIST = \
	Makefile.depends

libana_ADDLIB = ../../lib/lib.la

SUFFIXES = .fpp
DEPENDS = Makefile.depends

generated_f90 = $(foreach src,$(libana_la_SOURCES),$(if $(filter %.fpp,$(src)),$(src:.fpp=.f90)))

CLEANFILES = $(generated_f90) *~ *.mod

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
	$(LIBTOOL) --mode=compile $(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../../lib $(FCFLAGS) -c $*.f90
else
	$(LIBTOOL) --mode=compile $(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../../lib $(FCFLAGS) -c $*.f90
endif

depend: clean_depend
	python $(top_srcdir)/fortdep.py --objext lo *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)

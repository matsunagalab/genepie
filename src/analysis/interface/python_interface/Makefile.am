SUBDIRS := \
        ../../libana \
        ../../converter/crd_convert \
        ../../trj_analysis/trj_analysis \
        ../../trj_analysis/diffusion_analysis \
        ../../trj_analysis/drms_analysis \
        ../../trj_analysis/hb_analysis \
        ../../trj_analysis/msd_analysis \
        ../../trj_analysis/rg_analysis \
        ../../trj_analysis/rmsd_analysis \
        ../../clustering/kmeans_clustering \
        ../mbar_analysis \
        ../../free_energy/wham_analysis \
        ../../mode_analysis/avecrd_analysis \
        .

#noinst_LTLIBRARIES = libpython_interface.la
lib_LTLIBRARIES = libpython_interface.la

libpython_interface_la_SOURCES = \
        conv_f_c_util.fpp \
        s_molecule_c_mod.fpp \
        ctrl_c_mod.fpp \
        define_molecule.fpp \
        s_trajectories_c_mod.fpp \
        crd_convert_convert.fpp \
        crd_convert_c_mod.fpp \
        trj_analysis_c_mod.fpp \
        trj_analysis_analysis.fpp \
        rg_analysis_c_mod.fpp \
        rg_analysis_analysis.fpp \
        ra_analysis_c_mod.fpp \
        ra_analysis_analysis.fpp \
        dr_analysis_c_mod.fpp \
        dr_analysis_analysis.fpp \
        ma_analysis_c_mod.fpp \
        ma_analysis_analysis.fpp \
        hb_analysis_c_mod.fpp \
        hb_analysis_analysis.fpp \
        diffusion_analysis_main_c_mod.fpp \
        diffusion_analysis_analyze.fpp \
        aa_analysis_c_mod.fpp \
        aa_analysis_analysis.fpp \
        wa_analysis_c_mod.fpp \
        wa_analysis_analysis.fpp \
        mbar_analysis_c_mod.fpp \
        mbar_analysis_analysis.fpp \
        kc_analysis_c_mod.fpp \
        kc_analysis_analysis.fpp \
        dynamic_string_mod.fpp \
        internal_file_type_mod.fpp \
        export_molecule_c_mod.fpp

EXTRA_DIST = Makefile.depends

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = *.f90 *~ *.mod

# Ensure Position-Independent Code is generated
AM_FCFLAGS = -fPIC

libpython_interface_la_LDFLAGS = -module -avoid-version -export-dynamic -shared -no-undefined -Wl,-undefined,dynamic_lookup -L$(libdir) 
libpython_interface_la_LDFLAGS += -Wl,-undefined,dynamic_lookup
libpython_interface_la_LDFLAGS += $(shell mpif90 --showme:link)
libpython_interface_la_LDFLAGS += -lgfortran -lmpi -lmpi_mpifh
libpython_interface_la_LIBADD = \
        $(top_srcdir)/src/lib/lib.a \
        $(top_srcdir)/src/analysis/libana/libana.a
#        -llapack -lblas -lgfortran -lmpi -lmpi_mpifh
libpython_interface_la_FCFLAGS = \
	$(FCFLAGS) $(PYTHON_CPPFLAGS) \
	-I$(top_srcdir)/src/ \
	-I$(top_srcdir)/src/lib \
	-I$(top_srcdir)/src/analysis/libana \
	-I$(top_srcdir)/src/analysis/converter/crd_convert \
	-I$(top_srcdir)/src/analysis/trj_analysis/trj_analysis \
	-I$(top_srcdir)/src/analysis/trj_analysis/rg_analysis \
	-I$(top_srcdir)/src/analysis/trj_analysis/rmsd_analysis \
	-I$(top_srcdir)/src/analysis/trj_analysis/drms_analysis \
	-I$(top_srcdir)/src/analysis/trj_analysis/msd_analysis \
	-I$(top_srcdir)/src/analysis/trj_analysis/hb_analysis \
	-I$(top_srcdir)/src/analysis/trj_analysis/diffusion_analysis \
	-I$(top_srcdir)/src/analysis/mode_analysis/avecrd_analysis \
	-I$(top_srcdir)/src/analysis/clustering/kmeans_clustering \
	-I$(top_srcdir)/src/analysis/free_energy/wham_analysis \
	-I$(top_srcdir)/src/analysis/interface/mbar_analysis


# Add these lines
CCLD = $(FC)
LINK = $(LIBTOOL) --tag=FC --mode=link $(FC) $(AM_CFLAGS) $(CFLAGS) \
       $(AM_FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

-include $(DEPENDS)

.fpp.f90:
if USEKCOMP
	cp $< $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $< $*.f90
endif

.f90.lo:
if USEKCOMP
	mv $*.cpp.f90 $*.f90
endif
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		$(libpython_interface_la_FCFLAGS) -c $*.f90

depend: clean_depend
	python $(top_srcdir)/fortdep.py --objext lo --objprefix libpython_interface_la- *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

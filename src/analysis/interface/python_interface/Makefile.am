if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB1)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana \
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

EXTRA_DIST = \
	Makefile.depends

SUFFIXES = .fpp
DEPENDS = Makefile.depends

# Ensure Position-Independent Code is generated
AM_FCFLAGS = -fPIC

libpython_interface_la_LDFLAGS = -module -avoid-version -export-dynamic -shared -no-undefined
libpython_interface_la_LDFLAGS += -Wl,-undefined,dynamic_lookup
libpython_interface_la_LDFLAGS += $(shell mpif90 --showme:link)
libpython_interface_la_LIBADD = \
		../../converter/crd_convert/libcrd_convert.a \
		../../trj_analysis/trj_analysis/libtrj_analysis.a \
		../../trj_analysis/rg_analysis/librg_analysis.a \
		../../trj_analysis/rmsd_analysis/librmsd_analysis.a \
		../../trj_analysis/drms_analysis/libdrms_analysis.a \
		../../trj_analysis/msd_analysis/libmsd_analysis.a \
		../../trj_analysis/hb_analysis/libhb_analysis.a \
		../../trj_analysis/diffusion_analysis/libdiffusion_analysis.a \
		../../mode_analysis/avecrd_analysis/libavecrd_analysis.a \
		../../clustering/kmeans_clustering/libkmeans_clustering.a \
		../../free_energy/wham_analysis/libwham_analysis.a \
		../mbar_analysis/libmbar_analysis.a \
		../../libana/libana.a ../../../lib/lib.a ../../../ \
		$(INT_LAPACK_LIB) $(LDFLAGS) \
		-lgfortran 

libpython_interface_la_FCFLAGS = $(FCFLAGS) $(PYTHON_CPPFLAGS)

# Add these lines
CCLD = $(FC)
LINK = $(LIBTOOL) --tag=FC --mode=link $(FC) $(AM_CFLAGS) $(CFLAGS) \
       $(AM_FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

# FC = mpif90
# CC = mpicc

.fpp.lo:
if USEKCOMP
	cp $< $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana \
		-I../../converter/crd_convert \
		-I ../../trj_analysis/trj_analysis \
		-I ../../trj_analysis/rg_analysis \
		-I ../../trj_analysis/rmsd_analysis \
		-I ../../trj_analysis/drms_analysis \
		-I ../../trj_analysis/msd_analysis \
		-I ../../trj_analysis/hb_analysis \
		-I ../../trj_analysis/diffusion_analysis \
		-I ../../mode_analysis/avecrd_analysis \
		-I ../../clustering/kmeans_clustering \
		-I ../../free_energy/wham_analysis \
		-I ../mbar_analysis \
		$(libpython_interface_la_FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $< $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana \
		-I../../converter/crd_convert \
		-I ../../trj_analysis/trj_analysis \
		-I ../../trj_analysis/rg_analysis \
		-I ../../trj_analysis/rmsd_analysis \
		-I ../../trj_analysis/drms_analysis \
		-I ../../trj_analysis/msd_analysis \
		-I ../../trj_analysis/hb_analysis \
		-I ../../trj_analysis/diffusion_analysis \
		-I ../../mode_analysis/avecrd_analysis \
		-I ../../clustering/kmeans_clustering \
		-I ../../free_energy/wham_analysis \
		-I ../mbar_analysis \
		$(libpython_interface_la_FCFLAGS) -c $*.f90
endif

depend: clean_depend
	python ../../../../fortdep.py *.fpp | sed -e 's/\.o/.lo/g' > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)

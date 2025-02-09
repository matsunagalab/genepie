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
		   ../../free_energy/mbar_analysis \
		   ../../free_energy/wham_analysis \
		   ../../mode_analysis/avecrd_analysis \
		   .

lib_LTLIBRARIES = libpython_interface.la

libpython_interface_la_SOURCES = \
        your_module.fpp \
        conv_f_c_util.fpp \
        s_molecule_c_mod.fpp \
        define_molecule.fpp \
        s_trajectories_c_mod.fpp \
        crd_convert_convert.fpp \
        crd_convert_c_mod.fpp \
        trj_analysis_c_mod.fpp \
        trj_analysis_analysis.fpp\
        rg_analysis_c_mod.fpp \
        rg_analysis_analysis.fpp

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
		../../libana/libana.a ../../../lib/lib.a ../../../ \
		-llapack -lblas -lgfortran -lmpi -lmpi_mpifh
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
		$(libpython_interface_la_FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $< $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana \
		-I../../converter/crd_convert \
		-I ../../trj_analysis/trj_analysis \
		-I ../../trj_analysis/rg_analysis \
		$(libpython_interface_la_FCFLAGS) -c $*.f90
endif

depend: clean_depend
	python ../../../../fortdep.py *.fpp | sed -e 's/\.o/.lo/g' > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)

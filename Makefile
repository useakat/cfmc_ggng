CC            = cc

CFLAGS        = -O

FC            = g77

INCLUDES      = -I./ -I./inc -I./helas3

FFLAGS	      = -O3 $(INCLUDES) -ffixed-line-length-132

INSTALL	      = install

LIBDIR = -L./lib 

LDFLAGS	      = $(LIBDIR)

#LIBS_SM = -lsmlib
LIBS_HELAS = -ldhelas3 -ldhelas5 -lmylib -lcfqcdv5
LIBS_BASES = -lbases50_xhsave
#LIBS_BASES = -lbases50_xhsave_max

LIBS          = $(LIBS_HELAS) \
		$(LIBS_BASES)
#		$(LIBS_SM)

MAKEFILE      = Makefile

OBJS          = bfunc.o userin.o bfunc_ini.o \
		rambo.o ps2bd.o ps2bdt.o ps3bd.o ps4bd.o \
		pshk2kai.o pshk_divps_dely.o pshk_v7.o pshk_funcs.o \
		ps_rambo.o wrap_ps.o wrap_pshk.o \
		get_chjacob_diagram_6ch.o get_chjacob_diagram_1tch.o \
		get_chjacob_diagram_tch.o get_chjacob_diagram_sch.o \
		get_chjacob_diagram_2ch.o \
		matrix_co1.f flow1.o flow2.o matrix_mg.o \
		matrix_cfqcd.o matrix_mg_per.o \
                wavsort.o mom_perm.o calc_jacob.o conjugate_cf.o \
		momconv.o randperm_finmom.o get_momenta.o impose_cuts.o \
		GenMom.o GetCF.o GenWgt.o GetChWgt.o get_amp2.o fill_histo.o \
		set_histo.o save_histo.o gen_basiccf.o check_cf.o \
		alfas_functions.o hmffrd.o bsffrd.o tolow.o comsep.o \
		leng.o pdf.o ctq6pdf.o opendata.o get_amp2_cfz.o

OBJS_BASES = mainb.o $(OBJS)
OBJS_SPRING = mains.o spevnt.o $(OBJS)
OBJS_MAKE = make_files.o make_cparam.o
OBJS_TEST = test.o pshk_v7.o pshk_funcs.o ps2bd.o

.f.o:
	$(FC) -c $(FFLAGS) $<

bases:  $(OBJS_BASES)
	$(FC) $(LDFLAGS) $(OBJS_BASES) $(LIBS) -o $@

spring:  $(OBJS_SPRING)
	$(FC) $(LDFLAGS) $(OBJS_SPRING) $(LIBS) -o $@

get_numch: get_numch.o
	 $(FC) $(LDFLAGS) get_numch.o $(LIBS) -o $@

get_nevent: get_nevent.o
	 $(FC) $(LDFLAGS) get_nevent.o $(LIBS) -o $@

combine_dat: combine_dat.o
	 $(FC) $(LDFLAGS) combine_dat.o $(LIBS) -o $@

makefiles: $(OBJS_MAKE)
	 $(FC) $(LDFLAGS) $(OBJS_MAKE) $(LIBS) -o $@


test: $(OBJS_TEST)
	$(FC) $(LDFLAGS) $(OBJS_TEST) $(LIBS) -o $@

clean:
	@rm *.o *~ *# core* xsec_* spring nevents_total.dat fort.* bases nevt_* *ajob* bases_err* combine_dat get_numch get_nevent

clean_all:
	@rm -r -f *.o *~ *# *.lhe log_* *.bases *.tdr fort* jobstatus.dat bases spring \
	get_numch get_nevent test nohup.out Events/* a.out *.dat *ajob* plots/* bases_plots/*

###

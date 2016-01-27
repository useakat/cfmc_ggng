CC            = cc

CFLAGS        = -O

FC            = g77

INCLUDES      = -I./ -I../hml/inc -I../helas3

FFLAGS	      = -O3 $(INCLUDES) -ffixed-line-length-132

# EXTHDRS	      = ../hml/inc/bsfile.inc \
# 		../hml/inc/bsxdim.inc \
# 		../hml/inc/crndm.inc \
# 		../hml/inc/fvevnt.inc \
# 		../hml/inc/gnbeam.inc \
# 		../hml/inc/gudat1.inc \
# 		../hml/inc/gujets.inc \
# 		../hml/inc/hmparm.inc \
# 		../hml/inc/hmunit.inc \
# 		../hml/inc/smcupl.inc \
# 		../hml/inc/smptab.inc

INSTALL	      = install

LIBDIR = -L../lib 

LDFLAGS	      = $(LIBDIR)

#LIBS_SM = -lsmlib
#LIBS_HML = -lhml -lhmlutil
LIBS_HELAS = -ldhelas3 -ldhelas5 -lmylib -lcfqcdv5
LIBS_BASES = -lbases50_xhsave
#LIBS_LIB = -llibrary

# LIBS          = $(LIBS_SM) \
# 		$(LIBS_HML) \
# 		$(LIBS_HELAS) \
# 		$(LIBS_LIB) \
# 		$(LIBS_BASES)
LIBS          = $(LIBS_HELAS) \
		$(LIBS_BASES)

MAKEFILE      = Makefile

OBJS          = find_chcf.o get_cflow.o program_per.o shift.o sort_cflow.o sort.o \
		check_cf.o check_atm.o analyse_atm.o gen_basiccf_all.o analyse_optch.o \
		analyse_perch.o get_chmom.o

OBJS1          = find_chcf.o get_cflow.o find_perch.o shift.o sort_cflow.o sort.o \
		check_cf.o check_atm.o analyse_atm.o gen_basiccf.o analyse_perch.o \
		write_chinfo.o analyse_optch.o get_chmom.o

.f.o:
	$(FC) -c $(FFLAGS) $<

find_cf:  $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

find_perch:  $(OBJS1)
	$(FC) $(LDFLAGS) $(OBJS1) $(LIBS) -o $@

#test: test.o mom_perm.o
#	 $(FC) $(LDFLAGS) test.o mom_perm.o $(LIBS) -o $@

clean:
	@rm *.o *~ *#

#clean_all:
#	@rm -r -f *.o *~ *# *.lhe log_* *.bases *.tdr fort* jobstatus.dat bases spring \
	get_numch get_nevent test nohup.out Events/* a.out *.dat *ajob* plots/* bases_plots/*

###

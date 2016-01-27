FC            = g77
FFLAGS        = -O
DEST	      = ../lib
INSTALL	      = install
LIBRARY	      = libmylib.a
OBJS  = fact.o ipnext.o simpson.o ps2b.o etaorder.o \
	order.o rsort.o write.o cut_chk2.o fvdot.o pkinvars.o ptorder_chk.o \
	runge2.o rambo.o get_momenta.o switch_mom.o kin_func.o set_pdf.o get_std.o \
	dacosh.o es4dot.o es4sq.o esbeta.o ordering.o data_array.o mom_check.o

install:	$(LIBRARY)
	        @echo Installing $(LIBRARY) in $(DEST)
	        @if [ $(DEST) != . ]; then \
	        (rm -f $(DEST)/$(LIBRARY); $(INSTALL) -c $(LIBRARY) $(DEST)); fi

$(LIBRARY):	$(OBJS)
		@echo  "Loading $(LIBRARY) ... "
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:
	rm -f *.o *.a *~

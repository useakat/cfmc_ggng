FC            = g77
FFLAGS        = -O
DEST	      = ../lib
INSTALL	      = install
LIBRARY	      = libcfqcdv5.a
OBJS  = FFV1_0.o FFV1_1.o FFV1_2.o FFV1_3.o VVV1_0.o VVV1_1.o VVV1_2.o VVV1_3.o VVVV1_0.o \
	VVVV1_1.o VVVV3_0.o VVVV3_1.o VVVV4_0.o VVVV4_1.o \
	cfsort.o coldel.o colfac.o colfact2.o ggggcf.o GLUON4_0.o GLUON4_1.o \
	GLUON4_4.o GLUON5_0.o GLUON5_1.o GLUON5_5.o GLUON6_0.o GLUON6_1.o GLUON6_6.o \
	GLUON7_0.o GLUON7_1.o GLUON7_7.o GLUON8_0.o jgggcf.o sumw.o 

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

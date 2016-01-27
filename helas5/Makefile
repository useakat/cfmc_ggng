FC = g77
FFLAGS        = -O
DEST	      = ../lib
LIBRARY	      = libdhelas5.a
OBJS = FFV1_0.o FFV1_1.o FFV1_2.o FFV1_3.o VVV1_0.o VVV1_1.o VVV1_2.o VVV1_3.o VVVV1_0.o VVVV1_1.o \
	   VVVV3_0.o VVVV3_1.o VVVV4_0.o VVVV4_1.o sumF2.o sumV2.o FFV2_5_0.o FFV2_5_3.o FFV5_0.o \
	FFV2_0.o FFV5_3.o FFV2_3.o sumV4.o sumV5.o sumV3.o sumF3.o

install:	$(LIBRARY)
	        @echo Installing $(LIBRARY) in $(DEST)
	        @if [ $(DEST) != . ]; then \
	        (rm -f $(DEST)/$(LIBRARY); install -c $(LIBRARY) $(DEST)); fi

$(LIBRARY):	$(OBJS)
		@echo  "Loading $(LIBRARY) ... "
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:
	rm -f *.o *.a *~
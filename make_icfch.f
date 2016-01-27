      subroutine make_icfch_1(nfin)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
      include 'opt_info.inc'
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nfin
C     LOCAL VARIABLES 
      integer i
      integer perpos1,nperch1
      integer icf_perch1(1000)
      integer optpos1,noptch1
      integer icf_optch1(1000)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      perpos1 = 1
      do i = 1,nch_per
         if (perch_type(i).eq.0) then
            icf_perch1(perpos1) = i
            perpos1 = perpos1 +1
         endif
      enddo

      optpos1 = 1
      do i = 1,nch_opt
         if (optch_type(i).eq.0) then
            icf_optch1(optpos1) = i
            optpos1 = optpos1 +1
         endif
      enddo

      nperch1 = ntch_per

      noptch1 = ntch_opt

      open(1,file='icf_ch.inc',status='replace')
      write(1,*) "     integer icf_perch1(",nperch1,")"
      call data_array(icf_perch1,nperch1,"icf_perch1",10,1)
      write(1,*)
      write(1,*) "     integer icf_optch1(",noptch1,")"
      call data_array(icf_optch1,noptch1,"icf_optch1",10,1)
      close(1)

      return
      end


      subroutine make_icfch_2(nfin)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
      include 'opt_info.inc'
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nfin
C     LOCAL VARIABLES 
      integer i
      integer perpos1,perpos2,nperch1,nperch2
      integer icf_perch1(1000),icf_perch2(1000)
      integer optpos1,optpos2,noptch1,noptch2
      integer icf_optch1(1000),icf_optch2(1000)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      perpos1 = 1
      perpos2 = 1
      do i = 1,nch_per
         if (perch_type(i).eq.0) then
            icf_perch1(perpos1) = i
            perpos1 = perpos1 +1
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
         elseif (perch_type(i).eq.1) then
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
         endif
      enddo

      optpos1 = 1
      optpos2 = 1
      do i = 1,nch_opt
         if (optch_type(i).eq.0) then
            icf_optch1(optpos1) = i
            optpos1 = optpos1 +1
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
         elseif (optch_type(i).eq.1) then
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
         endif
      enddo

      nperch1 = ntch_per
      nperch2 = ntch_per +nsch_per

      noptch1 = ntch_opt
      noptch2 = ntch_opt +nsch_opt

      open(1,file='icf_ch.inc',status='replace')
      write(1,*) "     integer icf_perch1(",nperch1,")"
      write(1,*) "     integer icf_perch2(",nperch2,")"
      call data_array(icf_perch1,nperch1,"icf_perch1",10,1)
      call data_array(icf_perch2,nperch2,"icf_perch2",10,1)
      write(1,*)
      write(1,*) "     integer icf_optch1(",noptch1,")"
      write(1,*) "     integer icf_optch2(",noptch2,")"
      call data_array(icf_optch1,noptch1,"icf_optch1",10,1)
      call data_array(icf_optch2,noptch2,"icf_optch2",10,1)
      close(1)

      return
      end


      subroutine make_icfch_3(nfin)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
      include 'opt_info.inc'
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nfin
C     LOCAL VARIABLES 
      integer i
      integer perpos1,perpos2,perpos3,nperch1,nperch2,nperch3
      integer icf_perch1(1000),icf_perch2(1000),icf_perch3(1000)
      integer optpos1,optpos2,optpos3,noptch1,noptch2,noptch3
      integer icf_optch1(1000),icf_optch2(1000),icf_optch3(1000)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      perpos1 = 1
      perpos2 = 1
      perpos3 = 1
      do i = 1,nch_per
         if (perch_type(i).eq.0) then
            icf_perch1(perpos1) = i
            perpos1 = perpos1 +1
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
         elseif (perch_type(i).eq.1) then
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
         elseif ((perch_type(i).eq.2).or.(perch_type(i).eq.20)) then
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
         endif
      enddo

      optpos1 = 1
      optpos2 = 1
      optpos3 = 1
      do i = 1,nch_opt
         if (optch_type(i).eq.0) then
            icf_optch1(optpos1) = i
            optpos1 = optpos1 +1
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
         elseif (optch_type(i).eq.1) then
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
         elseif ((optch_type(i).eq.2).or.(optch_type(i).eq.20)) then
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
         endif
      enddo

      nperch1 = ntch_per
      nperch2 = ntch_per +nsch_per
      nperch3 = ntch_per +nsch_per +n2sch_per

      noptch1 = ntch_opt
      noptch2 = ntch_opt +nsch_opt
      noptch3 = ntch_opt +nsch_opt +n2sch_opt

      open(1,file='icf_ch.inc',status='replace')
      write(1,*) "     integer icf_perch1(",nperch1,")"
      write(1,*) "     integer icf_perch2(",nperch2,")"
      write(1,*) "     integer icf_perch3(",nperch3,")"
      call data_array(icf_perch1,nperch1,"icf_perch1",10,1)
      call data_array(icf_perch2,nperch2,"icf_perch2",10,1)
      call data_array(icf_perch3,nperch3,"icf_perch3",10,1)
      write(1,*)
      write(1,*) "     integer icf_optch1(",noptch1,")"
      write(1,*) "     integer icf_optch2(",noptch2,")"
      write(1,*) "     integer icf_optch3(",noptch3,")"
      call data_array(icf_optch1,noptch1,"icf_optch1",10,1)
      call data_array(icf_optch2,noptch2,"icf_optch2",10,1)
      call data_array(icf_optch3,noptch3,"icf_optch3",10,1)
      close(1)

      return
      end


      subroutine make_icfch_4(nfin)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
      include 'opt_info.inc'
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nfin
C     LOCAL VARIABLES 
      integer i
      integer perpos1,perpos2,perpos3,perpos4
      integer nperch1,nperch2,nperch3,nperch4
      integer icf_perch1(1000),icf_perch2(1000),icf_perch3(1000)
      integer icf_perch4(1000)
      integer optpos1,optpos2,optpos3,optpos4
      integer noptch1,noptch2,noptch3,noptch4
      integer icf_optch1(1000),icf_optch2(1000),icf_optch3(1000)
      integer icf_optch4(1000)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      perpos1 = 1
      perpos2 = 1
      perpos3 = 1
      perpos4 = 1
      do i = 1,nch_per
         if (perch_type(i).eq.0) then
            icf_perch1(perpos1) = i
            perpos1 = perpos1 +1
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
         elseif (perch_type(i).eq.1) then
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
         elseif ((perch_type(i).eq.2).or.(perch_type(i).eq.20)) then
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
         elseif (perch_type(i).eq.3) then
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
         endif
      enddo

      optpos1 = 1
      optpos2 = 1
      optpos3 = 1
      optpos4 = 1
      do i = 1,nch_opt
         if (optch_type(i).eq.0) then
            icf_optch1(optpos1) = i
            optpos1 = optpos1 +1
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
         elseif (optch_type(i).eq.1) then
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
         elseif ((optch_type(i).eq.2).or.(optch_type(i).eq.20)) then
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
         elseif (optch_type(i).eq.3) then
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
         endif
      enddo

      nperch1 = ntch_per
      nperch2 = ntch_per +nsch_per
      nperch3 = ntch_per +nsch_per +n2sch_per
      nperch4 = ntch_per +nsch_per +n2sch_per +n3sch_per 

      noptch1 = ntch_opt
      noptch2 = ntch_opt +nsch_opt
      noptch3 = ntch_opt +nsch_opt +n2sch_opt
      noptch4 = ntch_opt +nsch_opt +n2sch_opt +n3sch_opt 

      open(1,file='icf_ch.inc',status='replace')
      write(1,*) "     integer icf_perch1(",nperch1,")"
      write(1,*) "     integer icf_perch2(",nperch2,")"
      write(1,*) "     integer icf_perch3(",nperch3,")"
      write(1,*) "     integer icf_perch4(",nperch4,")"
      call data_array(icf_perch1,nperch1,"icf_perch1",10,1)
      call data_array(icf_perch2,nperch2,"icf_perch2",10,1)
      call data_array(icf_perch3,nperch3,"icf_perch3",10,1)
      call data_array(icf_perch4,nperch4,"icf_perch4",10,1)
      write(1,*)
      write(1,*) "     integer icf_optch1(",noptch1,")"
      write(1,*) "     integer icf_optch2(",noptch2,")"
      write(1,*) "     integer icf_optch3(",noptch3,")"
      write(1,*) "     integer icf_optch4(",noptch4,")"
      call data_array(icf_optch1,noptch1,"icf_optch1",10,1)
      call data_array(icf_optch2,noptch2,"icf_optch2",10,1)
      call data_array(icf_optch3,noptch3,"icf_optch3",10,1)
      call data_array(icf_optch4,noptch4,"icf_optch4",10,1)
      close(1)

      return
      end


      subroutine make_icfch_5(nfin)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
      include 'opt_info.inc'
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nfin
C     LOCAL VARIABLES 
      integer i
      integer perpos1,perpos2,perpos3,perpos4,perpos5
      integer nperch1,nperch2,nperch3,nperch4,nperch5
      integer icf_perch1(1000),icf_perch2(1000),icf_perch3(1000)
      integer icf_perch4(1000),icf_perch5(1000)
      integer optpos1,optpos2,optpos3,optpos4,optpos5
      integer noptch1,noptch2,noptch3,noptch4,noptch5
      integer icf_optch1(1000),icf_optch2(1000),icf_optch3(1000)
      integer icf_optch4(1000),icf_optch5(1000)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      perpos1 = 1
      perpos2 = 1
      perpos3 = 1
      perpos4 = 1
      perpos5 = 1
      do i = 1,nch_per
         if (perch_type(i).eq.0) then
            icf_perch1(perpos1) = i
            perpos1 = perpos1 +1
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
            icf_perch5(perpos5) = i
            perpos5 = perpos5 +1            
         elseif (perch_type(i).eq.1) then
            icf_perch2(perpos2) = i
            perpos2 = perpos2 +1
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
            icf_perch5(perpos5) = i
            perpos5 = perpos5 +1            
         elseif ((perch_type(i).eq.2).or.(perch_type(i).eq.20)) then
            icf_perch3(perpos3) = i
            perpos3 = perpos3 +1
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
            icf_perch5(perpos5) = i
            perpos5 = perpos5 +1            
         elseif (perch_type(i).eq.3) then
            icf_perch4(perpos4) = i
            perpos4 = perpos4 +1            
            icf_perch5(perpos5) = i
            perpos5 = perpos5 +1            
         elseif (perch_type(i).eq.4) then
            icf_perch5(perpos5) = i
            perpos5 = perpos5 +1            
         endif
      enddo

      optpos1 = 1
      optpos2 = 1
      optpos3 = 1
      optpos4 = 1
      optpos5 = 1
      do i = 1,nch_opt
         if (optch_type(i).eq.0) then
            icf_optch1(optpos1) = i
            optpos1 = optpos1 +1
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
            icf_optch5(optpos5) = i
            optpos5 = optpos5 +1            
         elseif (optch_type(i).eq.1) then
            icf_optch2(optpos2) = i
            optpos2 = optpos2 +1
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
            icf_optch5(optpos5) = i
            optpos5 = optpos5 +1            
         elseif ((optch_type(i).eq.2).or.(optch_type(i).eq.20)) then
            icf_optch3(optpos3) = i
            optpos3 = optpos3 +1
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
            icf_optch5(optpos5) = i
            optpos5 = optpos5 +1            
         elseif (optch_type(i).eq.3) then
            icf_optch4(optpos4) = i
            optpos4 = optpos4 +1            
            icf_optch5(optpos5) = i
            optpos5 = optpos5 +1            
         elseif (optch_type(i).eq.4) then
            icf_optch5(optpos5) = i
            optpos5 = optpos5 +1            
         endif
      enddo

      nperch1 = ntch_per
      nperch2 = ntch_per +nsch_per
      nperch3 = ntch_per +nsch_per +n2sch_per
      nperch4 = ntch_per +nsch_per +n2sch_per +n3sch_per 
      nperch5 = ntch_per +nsch_per +n2sch_per +n3sch_per +n4sch_per 

      noptch1 = ntch_opt
      noptch2 = ntch_opt +nsch_opt
      noptch3 = ntch_opt +nsch_opt +n2sch_opt
      noptch4 = ntch_opt +nsch_opt +n2sch_opt +n3sch_opt 
      noptch5 = ntch_opt +nsch_opt +n2sch_opt +n3sch_opt +n4sch_opt 

      open(1,file='icf_ch.inc',status='replace')
      write(1,*) "     integer icf_perch1(",nperch1,")"
      write(1,*) "     integer icf_perch2(",nperch2,")"
      write(1,*) "     integer icf_perch3(",nperch3,")"
      write(1,*) "     integer icf_perch4(",nperch4,")"
      write(1,*) "     integer icf_perch5(",nperch5,")"
      call data_array(icf_perch1,nperch1,"icf_perch1",10,1)
      call data_array(icf_perch2,nperch2,"icf_perch2",10,1)
      call data_array(icf_perch3,nperch3,"icf_perch3",10,1)
      call data_array(icf_perch4,nperch4,"icf_perch4",10,1)
      call data_array(icf_perch5,nperch5,"icf_perch5",10,1)
      write(1,*)
      write(1,*) "     integer icf_optch1(",noptch1,")"
      write(1,*) "     integer icf_optch2(",noptch2,")"
      write(1,*) "     integer icf_optch3(",noptch3,")"
      write(1,*) "     integer icf_optch4(",noptch4,")"
      write(1,*) "     integer icf_optch5(",noptch5,")"
      call data_array(icf_optch1,noptch1,"icf_optch1",10,1)
      call data_array(icf_optch2,noptch2,"icf_optch2",10,1)
      call data_array(icf_optch3,noptch3,"icf_optch3",10,1)
      call data_array(icf_optch4,noptch4,"icf_optch4",10,1)
      call data_array(icf_optch5,noptch5,"icf_optch5",10,1)
      close(1)

      return
      end

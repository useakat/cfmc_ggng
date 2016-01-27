      program get_nevent
C     ****************************************************
C     
C     The driver for reading the number of channels from
C     cparam.inc.
C
C     By Yoshitaro Takaesu @KEK Nov.4 2012
C     ****************************************************
      implicitnone
C     
C     GLOBAL VARIABLES
C     
      include 'cparam.inc'
      integer i,j
      character*4 iii
      character*10 filename1
      character*11 filename2
      character*12 filename3
      real*8 xsec(nch),sigma(nch),total_xsec,total_error2,total_error
      real*8 total_sigma
      character*6 itotal_events,cchstart,cchp1
      integer total_events,ichstart,ichp1
C     ----------
C     BEGIN CODE
C     ----------
      call getarg(1,itotal_events)
      call getarg(2,cchstart)
      call getarg(3,cchp1)
      read(itotal_events,*) total_events
      read(cchstart,*) ichstart
      read(cchp1,*) ichp1

      open(2,file="log_xsec.dat",status="replace")
      open(3,file="log_nevent.dat",status="replace")

      total_xsec = 0d0
      total_error2 = 0d0
      write(2,*) "ch ","Xsec [pb]", " Error [pb]"
      do i = ichstart,ichp1-1
         write(iii,*) i
         if (i.lt.10) then
            filename1 = "xsec_"//iii(2:2)//".dat"
            open(10,file=filename1,status="old")
         elseif (i.lt.100) then
            filename2 = "xsec_"//iii(2:3)//".dat"
            open(10,file=filename2,status="old")
         elseif (i.lt.1000) then
            filename3 = "xsec_"//iii(2:4)//".dat"
            open(10,file=filename3,status="old")
         endif
         read(10,*) xsec(i),sigma(i)
         write(2,'(i3,1x,E10.4,1x,E10.4,a3,f5.2,a3)') i,xsec(i)
     &        ,sigma(i)," ( ",sigma(i)/xsec(i)*100,"% )"
         total_xsec = total_xsec + xsec(i)
         total_error2 = total_error2 +sigma(i)**2
         close(10)
      enddo
      total_sigma = dsqrt(total_error2)

      write(3,*) "ch ","nevents"
      do i = ichstart,ichp1-1
         write(iii,*) i
         if (i.lt.10) then
            filename1 = "nevt_"//iii(2:2)//".dat"
            open(20,file=filename1,status="replace")
         elseif (i.lt.100) then
            filename2 = "nevt_"//iii(2:3)//".dat"
            open(20,file=filename2,status="replace")
         elseif (i.lt.1000) then
            filename3 = "nevt_"//iii(2:4)//".dat"
            open(20,file=filename3,status="replace")
         endif
         write(20,*) int(xsec(i)/total_xsec*total_events) +1
         write(3,*) i,int(xsec(i)/total_xsec*total_events) +1
         close(20)
      enddo

      open(1,file="uwgt.dat",status="replace")
      write(1,*) total_xsec/dble(total_events)
      close(1)

      open(1,file="xsec_total.dat",status="replace")
      write(1,'(E10.4,1x,E10.4,1x,E10.4)') total_xsec,total_sigma
     &     ,total_sigma/total_xsec
      write(2,'(a14,E10.4,a4,E10.4,a2,E10.4,a6)') "total xsec = "
     &     ,total_xsec," +- ",total_sigma," (",total_sigma/total_xsec
     &     ,") [pb]"
      close(1)

      open(1,file="nevents_total.dat",status="replace")
      write(1,*) total_events
      write(3,*) "total events = ",total_events
      close(1)

      close(2)
      close(3)

      end

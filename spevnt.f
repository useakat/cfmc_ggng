      subroutine spevnt(iret)
      implicit none
      
      integer iret

c      include 'bsfile.inc'
      include 'cparam.inc'
      include 'lhe.inc'
      include 'bfunc.inc'

      integer ievent/ 0/
      save ievent

      iret = 0
      ievent = ievent + 1

c      IDUP(1) = 2
c      IDUP(2) = -2
c      IDUP(3) = 11
c      IDUP(4) = -11
c      IDUP(5) = 1
c      IDUP(6) = -1
      do i = 1,nini
         PUP(1,i) = Q(1,i)
         PUP(2,i) = Q(2,i)
         PUP(3,i) = Q(3,i)
         PUP(4,i) = Q(0,i)
         PUP(5,i) = 0d0
         VTIMUP(i) = 0
         SPINUP(i) = 1
         ISTUP(i) = -1
         MOTHUP(1,i) = 0
         MOTHUP(2,i) = 0
         ICOLUP(1,i) = 503
         ICOLUP(2,i) = 501
      enddo
      do i = nini+1,next
         PUP(1,i) = Q(1,i)
         PUP(2,i) = Q(2,i)
         PUP(3,i) = Q(3,i)
         PUP(4,i) = Q(0,i)
         PUP(5,i) = 0d0
         VTIMUP(i) = 0
         SPINUP(i) = 1
         ISTUP(i) = 1
         MOTHUP(1,i) = 1
         MOTHUP(2,i) = 2
         ICOLUP(1,i) = 503
         ICOLUP(2,i) = 501
      enddo

      open(1,file="uwgt.dat",status="old")
      read(1,*) XWGTUP
      close(1)
      
      write(luevtf,*) "<event>" 
      write(luevtf,50) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP     

      do i=1,nini+nfin
         write(luevtf,100) IDUP(i),ISTUP(i),MOTHUP(1,i)
     &        ,MOTHUP(2,i),ICOLUP(1,i),ICOLUP(2,i),PUP(1,i),PUP(2,i)
     &        ,PUP(3,i),PUP(4,i),PUP(5,i),VTIMUP(i),SPINUP(i)
      enddo
      write(luevtf,*) "</event>"      

 50   format(i2,i4,4(E16.7))
 100  format(i10,i5,i5,i5,i5,i5,e19.11,e19.11,
     &     e19.11,e19.11,e19.11,f3.0,f4.0)

      return
      end

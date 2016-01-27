      INTEGER I,J,K
      real*8 z(30)
      REAL*8 PMASS(NEXT)	
      integer ierr,ierr2
      real*8 pdf1(-7:7),pdf2(-7:7)
      real*8 pt(nfin),eta(nfin),dr,ptij
      integer fact
      external fact
      real*8 Q(0:3,next)
      common /Q/ Q
      integer ptevent,ptevent1
      common /ptevent/ ptevent,ptevent1
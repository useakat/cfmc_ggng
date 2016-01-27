      subroutine fill_histo(z,Q,amp2,bfunc)
      implicit none
C     GLOBAL VARIABLES
      include 'cparam.inc'
      include 'pshk.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 z(30),Q(0:3,next),amp2,bfunc
C     LOCAL VARIABLES
      integer i,j,k
      real*8 dist(12),OQ(0:3,next),m2summ2 
C     EXTERNAL FUNCTIONS
      real*8 pdeltar,ppt,py,peta,pphi,pmij,pktij
      external pdeltar,ppt,py,peta,pphi,pmij,pktij
C     ----------
C     BEGIN CODE
C     ---------
      m2summ2 = amp2/AAA
      call order_pt(nfin,Q(0,3),OQ(0,3))
      i = 1
      do k = 3,next
         dist(i) = ppt(OQ(0,k))
         i = i +1
      enddo
      do k = 3,next
         dist(i) = 5 +py(OQ(0,k))
         i = i +1
      enddo
      do k = 3,next-1
         do j = k+1,next
            dist(i) = pdeltar(OQ(0,k),OQ(0,j))
            i = i +1
         enddo
      enddo
      do k = 3,next-1
         do j = k+1,next
            dist(i) = pmij(OQ(0,k),OQ(0,j))
            i = i +1
         enddo
      enddo

      call xhfill(31,x1,bfunc)
      call xhfill(32,x2,bfunc)
      call xhfill(33,m2summ2,bfunc)
c      call dhfill(33,x1,x2,bfunc)
      do i=1,ndim
         call xhfill(33+i,z(i),bfunc)
      enddo
      do i = 1,ndist
         call xhfill(40+i,dist(i),bfunc)
      enddo

      return
      end

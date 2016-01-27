      subroutine userin

      implicit none

      include 'bases.inc'

      integer maxdim   
      parameter (maxdim=50)
      double precision xl(maxdim),xu(maxdim)
 
      integer ig(maxdim)

      integer ncall, itmx1, itmx2
      double precision acc1, acc2
      real*8 tau
      common /para/ tau

      integer idim

      ndim = 1
      nwild = 0
                
      xl(1) = 0.5*dlog(tau)
      xu(1) = -0.5*dlog(tau)
c      xl(1) = 0
c      xu(1) = 1
      ig(1) = 1

      call bssetd(ndim,nwild,xl,xu,ig)

      ncall = 10000
      itmx1 = 5
      itmx2 = 5
      acc1 = 0.01d0
      acc2 = 0.01d0      
      call bssetp(ncall, itmx1, itmx2, acc1, acc2)

*      call xhinit(1,-1.d0,1.d0,50,'d sigma / d x1')
*      call xhinit(2,-1.d0,1.d0,50,'d sigma / d x2')

*      call dhinit(100,-1.d0,1.d0,50,-1.d0,1.d0,50, 
*     &     'd sigma/ d x1/d x2')

      return
      end

      subroutine ps_rambo(rdn,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i
      real*8 rdn(30),Q(0:3,next),pmass(next),preQ(0:3,next),labQ1(0:3)
      real*8 labQ2(0:3),labQ12(0:3),prambo(4,next)
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat

      ierr = 0

      whatmin = nfin*ptcut
      taulmax = 0d0
      taulmin = 2*dlog(whatmin/sqrts)
      taul = taulmin +(taulmax -taulmin)*rdn(1)
      tau = dexp(taul)

      ymax = min(-0.5d0*taul,ayecut)
      ymin = max(0.5d0*taul,-ayecut)
      y = ymin +(ymax-ymin)*rdn(2)

      x1 = dsqrt(tau)*dexp(y)
      x2 = dsqrt(tau)*dexp(-y)
      shat = tau*ss

      labQ1(0) = x1*sqrts/2d0
      labQ1(1) = 0d0
      labQ1(2) = 0d0
      labQ1(3) = x1*sqrts/2d0
      labQ1(0) = x2*sqrts/2d0
      labQ1(1) = 0d0
      labQ1(2) = 0d0
      labQ1(3) = -x2*sqrts/2d0
      do i = 0,3
         labQ12(i) = labQ1(i) +labQ2(i)
      enddo

      do i = 1,next
         pmass(i) = 0d0
      enddo
      call rambo(next-2,shat,pmass(3),prambo,phffin)
      DO I=3, NEXT
         preQ(0,I)=PRAMBO(4,I-2)	
         preQ(1,I)=PRAMBO(1,I-2)
         preQ(2,I)=PRAMBO(2,I-2)
         preQ(3,I)=PRAMBO(3,I-2)	
      ENDDO
      do i = 1,next
         call boostx(preQ(0,i),labQ12,Q(0,i))
      enddo

      flux = 1/(2*shat)
      jacob = tau
      wpsn = (taulmax -taulmin)*(ymax -ymin)*phffin*flux*0.38937966d9

      return
      end

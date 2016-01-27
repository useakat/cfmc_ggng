      subroutine wrap_ps2bd(rnd,sgncos,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i,sgncos
      real*8 rnd(30),Q(0:3,next),labQ(0:3),jjacob,wgt1
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat
      real*8 jjacob1,jjacob2,jjacob3,wwgt1,wwgt2,wwgt3
      real*8 smin,smax

      ierr = 0

      if (ipdf.eq.0) then
         jjacob1 = 1d0
         wwgt1 = 1d0
         jjacob2 = 1d0
         wwgt2 = 1d0
         x1 = 1d0
         x2 = 1d0
      elseif (ipdf.eq.1) then
         whatmin = 2*ptcut
         taulmin = 2*dlog(whatmin/sqrts)
         taulmax = 2*dlog(sqrts/sqrts)
         taul = taulmin +(taulmax -taulmin)*rnd(1)
         tau = dexp(taul)
         jjacob1 = tau
         wwgt1 = taulmax -taulmin

c         ymax = min(-0.5d0*taul,ayecut)
c         ymin = max(0.5d0*taul,-ayecut)
         ymax = -0.5d0*taul
         ymin = 0.5d0*taul
         y = ymin +(ymax-ymin)*rnd(2)
         jjacob2 = 1d0
         wwgt2 = ymax -ymin
         
         x1 = dsqrt(tau)*dexp(y)
         x2 = dsqrt(tau)*dexp(-y)
         
      else
         ierr = 1
         return
      endif

      ebini1 = x1*sqrts/2d0
      ebini2 = x2*sqrts/2d0

      call momntx(ebini1,0d0, 1d0,0d0, Q(0,1))
      call momntx(ebini2,0d0,-1d0,pi, Q(0,2))
      do i = 0,3
         labQ(i) = Q(i,1) +Q(i,2)
      enddo

      if (ipdf.eq.0) then 
         call ps2bdt(rnd,labQ,mass(3),sgncos,Q(0,3),jjacob3
     &        ,wwgt3,ierr)
      elseif (ipdf.eq.1) then
         call ps2bdt(rnd(3),labQ,mass(3),sgncos,Q(0,3),jjacob3
     &        ,wwgt3,ierr)
      else
         write(99,*) "ERROR: ipdf should be 0 or 1. (wrap_ps2bd)"
         stop
      endif

      shat = x1*x2*sqrts**2
      flux = 1/(2*shat)
      jacob = jjacob1*jjacob2*jjacob3
      wpsn = wwgt1*wwgt2*wwgt3*flux*0.38937966d9

      return
      end


      subroutine wrap_ps2bdt(rnd,sgncos,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i,sgncos,ipw
      real*8 rnd(10),Q(0:3,next),labQ(0:3),jjacob,wgt1
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat
      real*8 jjacob1,jjacob2,jjacob3,jjacob4,jjacob5,jjacob6
      real*8 wwgt1,wwgt2,wwgt3,wwgt4,wwgt5,wwgt6
      real*8 smin1,smax1,smin2,smax2,diffcut
      real*8 s1,s1l,s1lmin,s1lmax
      real*8 s2,s2l,s2lmin,s2lmax
      real*8 smass(2),K(0:3,2)

      ierr = 0
      ipw = 2
      
      do i = 1,2
         if (mass(i+2).eq.0d0) then
            smass(i) = ptcut
         elseif (mass(i+2).gt.0d0) then
            smass(i) = mass(i+2)
         endif
      enddo

      if (ipdf.eq.0) then
         jjacob2 = 1d0
         wwgt2 = 1d0
         jjacob3 = 1d0
         wwgt3 = 1d0
         x1 = 1d0
         x2 = 1d0
      elseif (ipdf.eq.1) then
c         whatmin = 2*Qcut
c         taulmin = 2*dlog(whatmin/sqrts)
c         taulmax = 2*dlog(sqrts/sqrts)
c         taul = taulmin +(taulmax -taulmin)*rnd(1)
c         tau = dexp(taul)
c         jjacob1 = tau
c         wwgt1 = taulmax -taulmin

c         diffcut = 1d-8
         whatmin = smass(1) +smass(2)
         taulmin = 2*dlog(whatmin/sqrts)
         taulmax = 0d0
         taul = taulmin +(taulmax -taulmin)*rnd(1)
         tau = dexp(taul)
         jjacob1 = dexp(taul)
         wwgt1 = taulmax -taulmin

c         diffcut = 1d-8
c         whatmin = (s1 +s2)/sqrts**2
c         taulmin = 1d0/diffcut
c         taulmax = 1d0/(1d0 -whatmin +diffcut)
c         taul = taulmin +(taulmax -taulmin)*rnd(3)
c         tau = 1d0/taul +whatmin -diffcut
c         jjacob2 = -(1d0/taul)**2
c         wwgt2 = taulmax -taulmin

c         whatmin = 2*Qcut
c         taulmin = 2*dlog(whatmin/sqrts)
c         taulmax = 2*dlog(sqrts/sqrts)
c         taul = taulmin +(taulmax -taulmin)*rnd(1)
c         tau = dexp(taul)
c         jjacob1 = tau
c         wwgt1 = taulmax -taulmin


c         ymax = min(-0.5d0*taul,ayecut)
c         ymin = max(0.5d0*taul,-ayecut)
         ymax = -0.5d0*dlog(tau)
         ymin = 0.5d0*dlog(tau)
         y = ymin +( ymax -ymin )*rnd(2)
         jjacob2 = 1d0
         wwgt2 = ymax -ymin
         
         x1 = dsqrt(tau)*dexp(y)
         x2 = dsqrt(tau)*dexp(-y)
      else
         ierr = 1
         return
      endif

      ebini1 = x1*sqrts/2d0
      ebini2 = x2*sqrts/2d0

      call momntx(ebini1,0d0, 1d0,0d0, Q(0,1))
      call momntx(ebini2,0d0,-1d0,pi, Q(0,2))
      do i = 0,3
         labQ(i) = Q(i,1) +Q(i,2)
      enddo

      if (ipdf.eq.0) then 
         call ps2bdt(rnd(1),labQ,mass,sgncos,Q,jjacob3,wwgt3,ierr)
      elseif (ipdf.eq.1) then
         call ps2bdt(rnd(3),labQ,smass,sgncos,Q,jjacob3,wwgt3,ierr)
      else
         ierr = 1
         return
      endif

      shat = x1*x2*sqrts**2
      flux = 1/(2*shat)
      jacob = jjacob1*jjacob2*jjacob3
      wpsn = wwgt1*wwgt2*wwgt3*flux*0.38937966d9

      return
      end


      subroutine wrap_ps3bd(rnd,sgncos,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i,sgncos
      real*8 rnd(30),Q(0:3,next),labQ(0:3),jjacob,wgt1
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat
      real*8 jjacob1,jjacob2,jjacob3,wwgt1,wwgt2,wwgt3
      real*8 smin,smax

      ierr = 0

      if (ipdf.eq.0) then
         jjacob1 = 1d0
         wwgt1 = 1d0
         jjacob2 = 1d0
         wwgt2 = 1d0
         x1 = 1d0
         x2 = 1d0
      elseif (ipdf.eq.1) then
         whatmin = Qcut +ptcut
         taulmin = 2*dlog(whatmin/sqrts)
         taulmax = 2*dlog(sqrts/sqrts)
         taul = taulmin +(taulmax -taulmin)*rnd(1)
         tau = dexp(taul)
         jjacob1 = tau
         wwgt1 = taulmax -taulmin

c         ymax = min(-0.5d0*taul,ayecut)
c         ymin = max(0.5d0*taul,-ayecut)
         ymax = -0.5d0*taul
         ymin = 0.5d0*taul
         y = ymin +(ymax-ymin)*rnd(2)
         jjacob2 = 1d0
         wwgt2 = ymax -ymin
         
         x1 = dsqrt(tau)*dexp(y)
         x2 = dsqrt(tau)*dexp(-y)
         
      else
         ierr = 1
         return
      endif

      ebini1 = x1*sqrts/2d0
      ebini2 = x2*sqrts/2d0

      call momntx(ebini1,0d0, 1d0,0d0, Q(0,1))
      call momntx(ebini2,0d0,-1d0,pi, Q(0,2))
      do i = 0,3
         labQ(i) = Q(i,1) +Q(i,2)
      enddo

      smin = Qcut**2
      smax = Qmax**2

      if (ipdf.eq.0) then 
         call ps3bd2(rnd,labQ,mass(3),smin,smax,sgncos,Q(0,3),jjacob3
     &        ,wwgt3,ierr)
      elseif (ipdf.eq.1) then
         call ps3bd2(rnd(3),labQ,mass(3),smin,smax,sgncos,Q(0,3),jjacob3
     &        ,wwgt3,ierr)
      else
         ierr = 1
         return
      endif

      shat = x1*x2*sqrts**2
      flux = 1/(2*shat)
      jacob = jjacob1*jjacob2*jjacob3
      wpsn = wwgt1*wwgt2*wwgt3*flux*0.38937966d9

      return
      end


      subroutine wrap_ps4bd(rnd,sgncos,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i,sgncos
      real*8 rnd(10),Q(0:3,next),labQ(0:3),jjacob,wgt1
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat
      real*8 jjacob1,jjacob2,jjacob3,wwgt1,wwgt2,wwgt3
      real*8 smin1,smax1,smin2,smax2,taull

      ierr = 0

      if (ipdf.eq.0) then
         jjacob1 = 1d0
         wwgt1 = 1d0
         jjacob2 = 1d0
         wwgt2 = 1d0
         x1 = 1d0
         x2 = 1d0
      elseif (ipdf.eq.1) then
         whatmin = 2*Qcut
         taulmin = 2*dlog(whatmin/sqrts)
         taulmax = 2*dlog(sqrts/sqrts)
         taul = taulmin +(taulmax -taulmin)*rnd(1)
         tau = dexp(taul)
         jjacob1 = tau
         wwgt1 = taulmax -taulmin

c         whatmin = 2*Qcut
c         taulmin = 2*dlog(whatmin/sqrts)
c         taulmax = 2*dlog(sqrts/sqrts)
c         taul = taulmin +(taulmax -taulmin)*rnd(1)
c         tau = dexp(taul)
c         jjacob1 = tau
c         wwgt1 = taulmax -taulmin


c         ymax = min(-0.5d0*taul,ayecut)
c         ymin = max(0.5d0*taul,-ayecut)
         ymax = -0.5d0*taul
         ymin = 0.5d0*taul
         y = ymin +(ymax-ymin)*rnd(2)
         jjacob2 = 1d0
         wwgt2 = ymax -ymin
         
         x1 = dsqrt(tau)*dexp(y)
         x2 = dsqrt(tau)*dexp(-y)
      else
         ierr = 1
         return
      endif

      ebini1 = x1*sqrts/2d0
      ebini2 = x2*sqrts/2d0

      call momntx(ebini1,0d0, 1d0,0d0, Q(0,1))
      call momntx(ebini2,0d0,-1d0,pi, Q(0,2))
      do i = 0,3
         labQ(i) = Q(i,1) +Q(i,2)
      enddo

      smin1 = Qcut**2
      smax1 = Qmax**2
      smin2 = Qcut**2
      smax2 = Qmax**2

      if (ipdf.eq.0) then 
         call ps4bd22(rnd,labQ,mass(3),smin1,smax1,smin2,smax2,sgncos
     &        ,Q(0,3),jjacob3,wwgt3,ierr)
      elseif (ipdf.eq.1) then
         call ps4bd22(rnd(3),labQ,mass(3),smin1,smax1,smin2,smax2,sgncos
     &        ,Q(0,3),jjacob3,wwgt3,ierr)
      else
         ierr = 1
         return
      endif

      shat = x1*x2*sqrts**2
      flux = 1/(2*shat)
      jacob = jjacob1*jjacob2*jjacob3
      wpsn = wwgt1*wwgt2*wwgt3*flux*0.38937966d9

      return
      end


      subroutine wrap_ps4bd22(rnd,sgncos,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i,sgncos
      real*8 rnd(10),Q(0:3,next),labQ(0:3),jjacob,wgt1
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat
      real*8 jjacob1,jjacob2,jjacob3,wwgt1,wwgt2,wwgt3
      real*8 smin1,smax1,smin2,smax2

      ierr = 0

      if (ipdf.eq.0) then
         jjacob1 = 1d0
         wwgt1 = 1d0
         jjacob2 = 1d0
         wwgt2 = 1d0
         x1 = 1d0
         x2 = 1d0
      elseif (ipdf.eq.1) then
         whatmin = 2*Qcut
         taulmin = 2*dlog(whatmin/sqrts)
         taulmax = 2*dlog(sqrts/sqrts)
         taul = taulmin +(taulmax -taulmin)*rnd(1)
         tau = dexp(taul)
         jjacob1 = tau
         wwgt1 = taulmax -taulmin

c         ymax = min(-0.5d0*taul,ayecut)
c         ymin = max(0.5d0*taul,-ayecut)
         ymax = -0.5d0*taul
         ymin = 0.5d0*taul
         y = ymin +(ymax-ymin)*rnd(2)
         jjacob2 = 1d0
         wwgt2 = ymax -ymin
         
         x1 = dsqrt(tau)*dexp(y)
         x2 = dsqrt(tau)*dexp(-y)
      else
         ierr = 1
         return
      endif

      ebini1 = x1*sqrts/2d0
      ebini2 = x2*sqrts/2d0

      call momntx(ebini1,0d0, 1d0,0d0, Q(0,1))
      call momntx(ebini2,0d0,-1d0,pi, Q(0,2))
      do i = 0,3
         labQ(i) = Q(i,1) +Q(i,2)
      enddo

      smin1 = Qcut**2
      smax1 = Qmax**2
      smin2 = Qcut**2
      smax2 = Qmax**2

      if (ipdf.eq.0) then 
         call ps4bd22(rnd,labQ,mass(3),smin1,smax1,smin2,smax2,sgncos
     &        ,Q(0,3),jjacob3,wwgt3,ierr)
      elseif (ipdf.eq.1) then
         call ps4bd22(rnd(3),labQ,mass(3),smin1,smax1,smin2,smax2,sgncos
     &        ,Q(0,3),jjacob3,wwgt3,ierr)
      else
         ierr = 1
         return
      endif

      shat = x1*x2*sqrts**2
      flux = 1/(2*shat)
      jacob = jjacob1*jjacob2*jjacob3
      wpsn = wwgt1*wwgt2*wwgt3*flux*0.38937966d9

      return
      end


      subroutine wrap_ps4bd22_2(rnd,sgncos,Q,ierr)
      implicitnone

      include 'cparam.inc'
      include 'pshk.inc'

      integer ierr,i,sgncos,ipw
      real*8 rnd(10),Q(0:3,6),labQ(0:3),jjacob,wgt1
      real*8 taul,taulmin,tau,phfpdf,ymax,ymin,etacut,y,ebini1,ebini2
      real*8 taulmax,flux,phffin,whatmin,shat
      real*8 jjacob1,jjacob2,jjacob3,jjacob4,jjacob5,jjacob6
      real*8 wwgt1,wwgt2,wwgt3,wwgt4,wwgt5,wwgt6
      real*8 smin1,smax1,smin2,smax2,diffcut
      real*8 s1,s1l,s1lmin,s1lmax
      real*8 s2,s2l,s2lmin,s2lmax
      real*8 smass(2),K(0:3,2)

      ierr = 0
      ipw = 2

      smin1 = Qcut**2
      smax1 = (Qmax-Qcut)**2

      if (ipw.eq.0) then
         s1lmin = smin1
         s1lmax = smax1
         s2lmin = smin2
         s2lmax = smax2
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = s1l
         s2 = s2l
         wwgt1 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jjacob1 = 1d0
      elseif(ipw.eq.1) then
         s1lmin = dsqrt(smin1)
         s1lmax = dsqrt(smax1)
         s2lmin = dsqrt(smin2)
         s2lmax = dsqrt(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = s1l**2
         s2 = s2l**2
         wwgt1 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jjacob1 = 2*dsqrt(s1)*2*dsqrt(s2)
      elseif(ipw.eq.2) then
         s1lmin = dlog(smin1)
         s1lmax = dlog(smax1)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s1 = dexp(s1l)
         smin2 = Qcut**2
         smax2 = (Qmax -dsqrt(s1))**2
         s2lmin = dlog(smin2)
         s2lmax = dlog(smax2)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s2 = dexp(s2l)
         wwgt1 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jjacob1 = s1*s2
      elseif(ipw.eq.3) then
         s1lmin = 1d0/dsqrt(smin1)
         s1lmax = 1d0/dsqrt(smax1)
         s2lmin = 1d0/dsqrt(smin2)
         s2lmax = 1d0/dsqrt(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = 1d0/s1l**2
         s2 = 1d0/s2l**2
         wwgt1 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jjacob1 = (-2*s1**1.5d0)*(-2*s2**1.5d0)
      elseif(ipw.eq.4) then
         s1lmin = 1d0/smin1
         s1lmax = 1d0/smax1
         s2lmin = 1d0/smin2
         s2lmax = 1d0/smax2
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = 1d0/s1l
         s2 = 1d0/s2l
         wwgt1 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jjacob1 = (-s1**2)*(-s2**2)
      endif
      
      smass(1) = dsqrt(s1)
      smass(2) = dsqrt(s2)

      if (ipdf.eq.0) then
         jjacob2 = 1d0
         wwgt2 = 1d0
         jjacob3 = 1d0
         wwgt3 = 1d0
         x1 = 1d0
         x2 = 1d0
      elseif (ipdf.eq.1) then
c         whatmin = 2*Qcut
c         taulmin = 2*dlog(whatmin/sqrts)
c         taulmax = 2*dlog(sqrts/sqrts)
c         taul = taulmin +(taulmax -taulmin)*rnd(1)
c         tau = dexp(taul)
c         jjacob1 = tau
c         wwgt1 = taulmax -taulmin

         diffcut = 1d-8
         whatmin = (s1 +s2)/sqrts**2
         taulmin = dlog(diffcut)
         taulmax = dlog(1d0 -whatmin +diffcut)
         taul = taulmin +(taulmax -taulmin)*rnd(3)
         tau = dexp(taul) +whatmin -diffcut
         jjacob2 = dexp(taul)
         wwgt2 = taulmax -taulmin

c         diffcut = 1d-8
c         whatmin = (s1 +s2)/sqrts**2
c         taulmin = 1d0/diffcut
c         taulmax = 1d0/(1d0 -whatmin +diffcut)
c         taul = taulmin +(taulmax -taulmin)*rnd(3)
c         tau = 1d0/taul +whatmin -diffcut
c         jjacob2 = -(1d0/taul)**2
c         wwgt2 = taulmax -taulmin

c         whatmin = 2*Qcut
c         taulmin = 2*dlog(whatmin/sqrts)
c         taulmax = 2*dlog(sqrts/sqrts)
c         taul = taulmin +(taulmax -taulmin)*rnd(1)
c         tau = dexp(taul)
c         jjacob1 = tau
c         wwgt1 = taulmax -taulmin


c         ymax = min(-0.5d0*taul,ayecut)
c         ymin = max(0.5d0*taul,-ayecut)
         ymax = -0.5d0*dlog(tau)
         ymin = 0.5d0*dlog(tau)
         y = ymin +( ymax -ymin )*rnd(4)
         jjacob3 = 1d0
         wwgt3 = ymax -ymin
         
         x1 = dsqrt(tau)*dexp(y)
         x2 = dsqrt(tau)*dexp(-y)
      else
         ierr = 1
         return
      endif

      ebini1 = x1*sqrts/2d0
      ebini2 = x2*sqrts/2d0

      call momntx(ebini1,0d0, 1d0,0d0, Q(0,1))
      call momntx(ebini2,0d0,-1d0,pi, Q(0,2))
      do i = 0,3
         labQ(i) = Q(i,1) +Q(i,2)
      enddo

      if (ipdf.eq.0) then 
         call ps2bdt(rnd(3),labQ,smass,sgncos,K,jjacob4,wwgt4,ierr)
         call ps2bd(rnd(5),K(0,1),mass(3),Q(0,3),jjacob5,wwgt5,ierr)
         call ps2bd(rnd(7),K(0,2),mass(5),Q(0,5),jjacob6,wwgt6,ierr)
      elseif (ipdf.eq.1) then
         call ps2bdt(rnd(5),labQ,smass,sgncos,K,jjacob4,wwgt4,ierr)
         call ps2bd(rnd(7),K(0,1),mass(3),Q(0,3),jjacob5,wwgt5,ierr)
         call ps2bd(rnd(9),K(0,2),mass(5),Q(0,5),jjacob6,wwgt6,ierr)
      else
         ierr = 1
         return
      endif

      shat = x1*x2*sqrts**2
      flux = 1/(2*shat)
      jacob = jjacob1*jjacob2*jjacob3*jjacob4*jjacob5*jjacob6
      wpsn = wwgt1*wwgt2*wwgt3*wwgt4*wwgt5*wwgt6*flux*0.38937966d9

      return
      end


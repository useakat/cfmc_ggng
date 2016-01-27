      real*8 function tjacob(n1,n2,ipw1,ipw2,P)
      implicit none

      include 'cparam.inc'

      integer n1,n2
      real*8 P(0:3,next),ipw1,ipw2,tjacob1,tjacob2

      real*8 pt
      external pt

      if (n2.eq.0) then      
         if (ipw1.eq.1) then
            tjacob = 1d0/pt(P(0,n1))**2
         elseif (ipw1.gt.1) then
            tjacob = 1d0/((pt(P(0,n1))**2)**ipw1/(ipw1-1d0))
         else
            write(luerr,*) "ERROR: ipw1 should be >= 1."
         endif
      else
         if (ipw1.eq.1) then
            tjacob1 = 1d0/pt(P(0,n1))**2
         elseif (ipw1.gt.1) then
            tjacob1 = 1d0/((pt(P(0,n1))**2)**ipw1/(ipw1-1d0))
         else
            write(luerr,*) "ERROR: ipw1 should be >= 1."
         endif
         if (ipw2.eq.1) then
            tjacob2 = 1d0/pt(P(0,n2))**2
         elseif (ipw2.gt.1) then
            tjacob2 = 1d0/((pt(P(0,n2))**2)**ipw2/(ipw2-1d0))
         else
            write(luerr,*) "ERROR: ipw2 should be >= 1."
         endif
         tjacob = tjacob1*tjacob2
      endif

      return
      end


      real*8 function sjacob(n1,n2,ipwm,P)
      implicit none

      include 'cparam.inc'

      integer n1,n2,i
      real*8 P(0:3,next),P12(0:3),ipwm,mmx

      real*8 mij
      external mij

      do i = 0,3
         P12(i) = P(i,n1) +P(i,n2)
      enddo

      if (ipwm.eq.0) then
         mmx = 1d0
      elseif (ipwm.eq.1) then
         mmx = mij(P(0,n1),P(0,n2))**2
      elseif (ipwm.gt.1) then
         mmx = (mij(P(0,n1),P(0,n2))**2)**ipwm/(ipwm-1d0)
      endif

      sjacob = 1d0/( mmx
c     &     *P(0,n1)/P12(0)*P(0,n2)/P12(0)/6d0 )
     &     )

      return
      end


      real*8 function sptjacob(n1,n2,ipwpt,ipwm,Q)
      implicit none

      include 'cparam.inc'
      include 'pshk.inc'

      integer n1,n2,i
      real*8 Q(0:3,next),Q12(0:3),ipwpt,mmx,prop,ipwm
      real*8 meff2

      real*8 pmij,ppt
      external pmij,ppt

      do i = 0,3
         Q12(i) = Q(i,n1) +Q(i,n2)
      enddo

      if (ipwm.eq.0) then
         mmx = 1d0
      elseif (ipwm.eq.1) then
         mmx = pmij(Q(0,n1),Q(0,n2))**2
      elseif (ipwm.gt.1) then
         mmx = (pmij(Q(0,n1),Q(0,n2))**2)**ipwm/(ipwm-1d0)
      endif
      meff2 = pmij(Q(0,n1),Q(0,n2))**2
c      meff2 = Qcut**2
      if (ipwpt.eq.0) then
         prop = 1d0
      elseif (ipwpt.eq.1) then
         prop = ppt(Q12)**2 +meff2
      elseif (ipwpt.gt.1) then
         prop = (ppt(Q12)**2 +meff2)**ipwpt/(ipwpt-1d0)
      elseif (ipwpt.lt.1) then
         prop = (ppt(Q12)**2 +meff2)**ipwpt/(1d0-ipwpt)
      endif

      sptjacob = 1d0/( mmx*prop
c     &        *Q(0,n1)/Q12(0)*Q(0,n2)/Q12(0)/6d0 )
     &        )

      return
      end


      real*8 function ps4bd22jacob(n1,n2,n3,n4,Q)
      implicit none

      include 'cparam.inc'
      include 'pshk.inc'

      integer n1,n2,n3,n4,i
      real*8 Q(0:3,next),Q12(0:3),s1,s2,prop,P1(0:3)
      real*8 shatmin,diffcut,tau,taul,shat,Q122(0:3)
      real*8 bbeta,a,adenom,costh

      real*8 pmij,ppt,esbeta,pcosth
      external pmij,ppt,esbeta,pcosth

      P1(0) = Q(0,1) +Q(0,2)
      P1(1) = -( Q(1,1) +Q(1,2) )
      P1(2) = -( Q(2,1) +Q(2,2) )
      P1(3) = -( Q(3,1) +Q(3,2) )

      do i = 0,3
         Q12(i) = Q(i,n1) +Q(i,n2)
      enddo
      call boostx(Q12,P1,Q122)

      s1 = pmij(Q(0,n1),Q(0,n2))**2
      s2 = pmij(Q(0,n3),Q(0,n4))**2

      diffcut = 1d-8
      tau = x1*x2
      shat = tau*ss
      shatmin = (s1+s2)/ss
      taul = tau -shatmin +diffcut

      bbeta = esbeta(s1/shat,s2/shat)
      adenom = 1d0 -(s1+s2)/shat +bbeta
      a = 4*s1*s2/shat**2/adenom
      costh = pcosth(Q122)
      prop = a/bbeta +1d0 -costh

      ps4bd22jacob = 1d0/( s1*s2*taul*prop )

      return
      end


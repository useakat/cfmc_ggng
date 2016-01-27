      subroutine vvtxxx(wm,wp,ft,g,vmass , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex wm(6)         : vector               flow-in  V
c       complex wp(6)         : vector               flow-out V~
c       complex ft(18)        : tensor               KK mode T
c       real    g             : coupling constant    -kappa/2
c       real    vmass         : V boson mass          m_V
c
c output:
c       complex vertex        : amplitude            gamma(wm,wp,tc)
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex wm(6), wp(6), ft(18),tc(6,4), vertex
      double precision g, vmass

      double complex T12, T13, T14, T23, T24, T34
      double complex V1V2, K1V2, K2V1
c     new
     &,K1V1,K2V2
c     new
      double complex TKK, TVV, TK1V2, TK2V1, dum
      double precision pwm(4), pwp(4), F

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      tc(1,1) = ft(1)
      tc(1,2) = ft(2)
      tc(1,3) = ft(3)
      tc(1,4) = ft(4)
      tc(2,1) = ft(5)
      tc(2,2) = ft(6)
      tc(2,3) = ft(7)
      tc(2,4) = ft(8)
      tc(3,1) = ft(9)
      tc(3,2) = ft(10)
      tc(3,3) = ft(11)
      tc(3,4) = ft(12)
      tc(4,1) = ft(13)
      tc(4,2) = ft(14)
      tc(4,3) = ft(15)
      tc(4,4) = ft(16)
      tc(5,1) = ft(17)
      tc(6,1) = ft(18)

      pwm(1) = dreal(wm(5))
      pwm(2) = dreal(wm(6))
      pwm(3) = dimag(wm(6))
      pwm(4) = dimag(wm(5))
      pwp(1) = dreal(wp(5))
      pwp(2) = dreal(wp(6))
      pwp(3) = dimag(wp(6))
      pwp(4) = dimag(wp(5))

      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      V1V2 =  wm(1)*wp(1) -  wm(2)*wp(2) -  wm(3)*wp(3) -  wm(4)*wp(4)
      K1V2 = pwm(1)*wp(1) - pwm(2)*wp(2) - pwm(3)*wp(3) - pwm(4)*wp(4)
      K2V1 = pwp(1)*wm(1) - pwp(2)*wm(2) - pwp(3)*wm(3) - pwp(4)*wm(4)
c     new
      K1V1 = pwm(1)*wm(1) - pwm(2)*wm(2) - pwm(3)*wm(3) - pwm(4)*wm(4)
      K2V2 = pwp(1)*wp(1) - pwp(2)*wp(2) - pwp(3)*wp(3) - pwp(4)*wp(4)
c     new

      F = pwm(1)*pwp(1) - pwm(2)*pwp(2) - pwm(3)*pwp(3) - pwm(4)*pwp(4)
      if ( vmass.ne.rZero ) then
         F = F + vmass**2
      end if

      TKK   = cZero
      TVV   = cZero
      TK1V2 = cZero
      TK2V1 = cZero

      do i = 1,4
         dum   = tc(i,i)*pwm(i)
         TKK   = TKK   + dum*pwp(i)
         TK1V2 = TK1V2 + dum*wp(i)
         dum   = tc(i,i)*wm(i)
         TVV   = TVV   + dum*wp(i)
         TK2V1 = TK2V1 + dum*pwp(i)
      end do

      TKK   = rTwo*TKK
      TVV   = rTwo*TVV
      TK1V2 = rTwo*TK1V2
      TK2V1 = rTwo*TK2V1

      TKK = TKK - T12*(pwm(1)*pwp(2) + pwm(2)*pwp(1))
     &          - T13*(pwm(1)*pwp(3) + pwm(3)*pwp(1))
     &          - T14*(pwm(1)*pwp(4) + pwm(4)*pwp(1))
     &          + T23*(pwm(2)*pwp(3) + pwm(3)*pwp(2))
     &          + T24*(pwm(2)*pwp(4) + pwm(4)*pwp(2))
     &          + T34*(pwm(3)*pwp(4) + pwm(4)*pwp(3))

      TK1V2 = TK1V2 - T12*(pwm(1)*wp(2) + pwm(2)*wp(1))
     &              - T13*(pwm(1)*wp(3) + pwm(3)*wp(1))
     &              - T14*(pwm(1)*wp(4) + pwm(4)*wp(1))
     &              + T23*(pwm(2)*wp(3) + pwm(3)*wp(2))
     &              + T24*(pwm(2)*wp(4) + pwm(4)*wp(2))
     &              + T34*(pwm(3)*wp(4) + pwm(4)*wp(3))

      TVV = TVV - T12*(wm(1)*wp(2) + wm(2)*wp(1))
     &          - T13*(wm(1)*wp(3) + wm(3)*wp(1))
     &          - T14*(wm(1)*wp(4) + wm(4)*wp(1))
     &          + T23*(wm(2)*wp(3) + wm(3)*wp(2))
     &          + T24*(wm(2)*wp(4) + wm(4)*wp(2))
     &          + T34*(wm(3)*wp(4) + wm(4)*wp(3))

      TK2V1 = TK2V1 - T12*(wm(1)*pwp(2) + wm(2)*pwp(1))
     &              - T13*(wm(1)*pwp(3) + wm(3)*pwp(1))
     &              - T14*(wm(1)*pwp(4) + wm(4)*pwp(1))
     &              + T23*(wm(2)*pwp(3) + wm(3)*pwp(2))
     &              + T24*(wm(2)*pwp(4) + wm(4)*pwp(2))
     &              + T34*(wm(3)*pwp(4) + wm(4)*pwp(3))

      vertex =  (tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4))*( K1V2*K2V1 - V1V2*F )
     &        + F*TVV + V1V2*TKK - K2V1*TK1V2 - K1V2*TK2V1

C      vertex = F*TVV + V1V2*TKK - K2V1*TK1V2 - K1V2*TK2V1

c     new, additonal gauge fixing term in Feyman gauge
      if ( vmass.eq.rZero ) then
         vertex = vertex 
     &+(tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4))*(K1V1*K1V2+K2V1*K2V2+K1V1*K2V2)
     &-K1V1*TK1V2-K2V2*TK2V1	   
      endif	
c     new    

      vertex = vertex * g

      return
      end

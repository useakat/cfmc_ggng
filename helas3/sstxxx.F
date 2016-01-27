      subroutine sstxxx(s1,s2,ft,g,smass , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two scalar and a Kaluza-Klein tensor boson.
c
c input:
c       complex s1(3)          : first scalar               s1
c       complex s2(3)          : second scalar              s2
c       complex ft(18)         : tensor               KK mode T
c       real    g              : coupling constant    -kappa/2
c       real    smass          : scalar mass          m_s
c
c output:
c       complex vertex         : amplitude            gamma(s1,s2,tc)
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex s1(3), s2(3),ft(18),tc(6,4), vertex
      double precision g, smass

      double complex T12, T13, T14, T23, T24, T34
      double complex TKK
      double precision ps1(4), ps2(4)
	integer i
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

      ps1(1) = dreal(s1(2))
      ps1(2) = dreal(s1(3))
      ps1(3) = dimag(s1(3))
      ps1(4) = dimag(s1(2))

      ps2(1) = -dreal(s2(2))
      ps2(2) = -dreal(s2(3))
      ps2(3) = -dimag(s2(3))
      ps2(4) = -dimag(s2(2))

      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)


      TKK   = cZero
    
      do i = 1,4
         TKK=TKK+tc(i,i)*ps1(i)*ps2(i)
      end do

      TKK   = rTwo*TKK

      TKK = TKK - T12*(ps1(1)*ps2(2) + ps1(2)*ps2(1))
     &          - T13*(ps1(1)*ps2(3) + ps1(3)*ps2(1))
     &          - T14*(ps1(1)*ps2(4) + ps1(4)*ps2(1))
     &          + T23*(ps1(2)*ps2(3) + ps1(3)*ps2(2))
     &          + T24*(ps1(2)*ps2(4) + ps1(4)*ps2(2))
     &          + T34*(ps1(3)*ps2(4) + ps1(4)*ps2(3))


      vertex = TKK+(tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4))
     &	*(smass**2-ps1(1)*ps2(1)+ps1(2)*ps2(2)
     &      +ps1(3)*ps2(3)+ps1(4)*ps2(4))

      vertex = vertex * g*s1(1)*s2(1)

      return
      end

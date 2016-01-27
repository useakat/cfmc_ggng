      subroutine hstxxx(ft,sc,g,smass,swidth , hts)
c
c This subroutine computes an off-shell scalar wavfunction from
c the scalar-scalar-Kaluza-Klein tensor boson coupling.
c
c input:
c       complex sc(3)          : input scalar               s
c       complex ft(18)         : tensor               KK mode T
c       real    g              : coupling constant    -kappa/2
c       real    smass          : mass  of OUTPUT scalar s'
c       real    swidth         : width of OUTPUT scalar s'
c output:
c       complex hts(3)         : scalar current       j(s':T,s)     
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex sc(3),ft(18),tc(6,4),hts(3)
      double precision g, smass,swidth

      double complex T12, T13, T14, T23, T24, T34
      double complex TKK
      double precision ps1(4), ps2(4)
      integer i
      double complex cZero, d
      double precision rZero, rTwo, pf2
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

      hts(2) = sc(2)+tc(5,1)
      hts(3) = sc(3)+tc(6,1)

      ps1(1) = dreal(sc(2))
      ps1(2) = dreal(sc(3))
      ps1(3) = dimag(sc(3))
      ps1(4) = dimag(sc(2))

      ps2(1) = dreal(hts(2))
      ps2(2) = dreal(hts(3))
      ps2(3) = dimag(hts(3))
      ps2(4) = dimag(hts(2))

      pf2 = ps2(1)**2 - ps2(2)**2 - ps2(3)**2 - ps2(4)**2
      
      if ( smass.gt.rZero ) then
         d = - g/dcmplx( pf2-smass**2, smass*swidth )
      else
         d = - g/dcmplx( pf2, rZero )
      end if

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

      hts(1) = TKK+(tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4))
     &	*(smass**2-ps1(1)*ps2(1)+ps1(2)*ps2(2)
     &      +ps1(3)*ps2(3)+ps1(4)*ps2(4))

      hts(1) = hts(1) * d*sc(1)

      return
      end

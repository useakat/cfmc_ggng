      subroutine ftoxxx(fo,ft,g,fmass,fwidth , fto)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a KK tensor boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex ft(18)         : input    tensor                      t
c       real    g              : coupling constant                   gtf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,t,f'|
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fo(6),ft(18),tc(6,4), fto(6)
      double precision g, fmass, fwidth

      double complex Tii,d
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision p(4), k(4), p2
      integer i

      double precision rZero, r2,rtwo
      parameter( rZero = 0.0d0, r2 = 2.0d0,rtwo=2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )


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

      fto(5) = fo(5) + tc(5,1)
      fto(6) = fo(6) + tc(6,1)

      p(1) = dreal(fto(5))
      p(2) = dreal(fto(6))
      p(3) = dimag(fto(6))
      p(4) = dimag(fto(5))
      p2 = p(1)**2 - p(2)**2 - p(3)**2 - p(4)**2

      k(1) = dreal(fo(5)) + p(1)
      k(2) = dreal(fo(6)) + p(2)
      k(3) = dimag(fo(6)) + p(3)
      k(4) = dimag(fo(5)) + p(4)

      T11 = tc(1,1)
      T22 = tc(2,2)
      T33 = tc(3,3)
      T44 = tc(4,4)
      Tii = T11-T22-T33-T44
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)


      if ( fmass.gt.rZero ) then
         d = - g/dcmplx( p2-fmass**2, fmass*fwidth )
      else
         d = - g/dcmplx( p2, rZero )
      end if

      fto(1) =fmass*(fmass*fo(1)*r2**2*Tii 
     &	- fo(4)*r2*Tii*(k(2) + ci*k(3))
     & - fo(3)*r2*Tii*(k(1) + k(4))) + 
     &    fmass*(fo(4)*(T12*k(1) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) - k(3)) - T24*k(4) - ci*T34*k(4)) + 
     &       fo(3)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4)))) + 
     &    (fmass*fo(4)*r2**2*Tii - fo(1)*r2*Tii*(-k(2) 
     &+ ci*k(3)) - fo(2)*r2*Tii*(k(1) + k(4)))*(p(2) + ci*p(3)) + 
     &    (fo(1)*(-(T12*k(1)) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) + k(3)) + T24*k(4) - ci*T34*k(4)) + 
     &       fo(2)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4))))*(p(2) + ci*p(3)) + 
     &    (fmass*fo(3)*r2**2*Tii - fo(2)*r2*Tii*(-k(2) 
     &- ci*k(3)) - fo(1)*r2*Tii*(k(1) - k(4)))*(p(1) + p(4)) + 
     &    (fo(2)*(-(T12*k(1)) - ci*T13*k(1) 
     &+ T23*(ci*k(2) + k(3)) + T24*k(4) + ci*T34*k(4)) + 
     &       fo(1)*(-(T12*k(2)) + T24*k(2) 
     &- T13*k(3) + T34*k(3) - T14*(k(1) + k(4))))*(p(1) + p(4)) + 
     &    rtwo*(fmass*(fo(4)*(-(T22*k(2)) 
     &- ci*T33*k(3)) + fo(3)*(T11*k(1) - T44*k(4))) + 
     &       (fo(1)*(T22*k(2) - ci*T33*k(3)) 
     &+ fo(2)*(T11*k(1) - T44*k(4)))*(p(2) + ci*p(3)) + 
     &       (fo(2)*(T22*k(2) + ci*T33*k(3))
     & + fo(1)*(T11*k(1) + T44*k(4)))*(p(1) + p(4)))

      fto(2) = fmass*(fmass*fo(2)*r2**2*Tii 
     &	- fo(3)*r2*Tii*(k(2) - ci*k(3))
     & - fo(4)*r2*Tii*(k(1) - k(4))) + 
     &    fmass*(fo(3)*(T12*k(1) - ci*T13*k(1)
     & + T23*(ci*k(2) - k(3)) - T24*k(4) + ci*T34*k(4)) + 
     &       fo(4)*(-(T12*k(2)) + T24*k(2) 
     &- T13*k(3) + T34*k(3) - T14*(k(1) + k(4)))) + 
     &    (fmass*fo(3)*r2**2*Tii - fo(2)*r2*Tii*(-k(2) 
     &- ci*k(3)) - fo(1)*r2*Tii*(k(1) - k(4)))*(p(2) - ci*p(3)) + 
     &    (fo(2)*(-(T12*k(1)) - ci*T13*k(1) 
     &+ T23*(ci*k(2) + k(3)) + T24*k(4) + ci*T34*k(4)) + 
     &       fo(1)*(-(T12*k(2)) + T24*k(2) - T13*k(3)
     & + T34*k(3) - T14*(k(1) + k(4))))*(p(2) - ci*p(3)) + 
     &    rtwo*(fmass*(fo(3)*(-(T22*k(2)) + ci*T33*k(3)) 
     &+ fo(4)*(T11*k(1) + T44*k(4))) + 
     &       (fo(2)*(T22*k(2) + ci*T33*k(3)) 
     &+ fo(1)*(T11*k(1) + T44*k(4)))*(p(2) - ci*p(3)) + 
     &       (fo(1)*(T22*k(2) - ci*T33*k(3)) 
     &+ fo(2)*(T11*k(1) - T44*k(4)))*(p(1) - p(4))) + 
     &    (fmass*fo(4)*r2**2*Tii - fo(1)*r2*Tii*(-k(2) 
     &+ ci*k(3)) - fo(2)*r2*Tii*(k(1) + k(4)))*(p(1) - p(4)) + 
     &    (fo(1)*(-(T12*k(1)) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) + k(3)) + T24*k(4) - ci*T34*k(4)) + 
     &       fo(2)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4))))*(p(1) - p(4))

      fto(3) = fmass*(fmass*fo(3)*r2**2*Tii 
     &	- fo(2)*r2*Tii*(-k(2) - ci*k(3)) 
     &- fo(1)*r2*Tii*(k(1) - k(4))) + 
     &    fmass*(fo(2)*(-(T12*k(1)) - ci*T13*k(1) 
     &+ T23*(ci*k(2) + k(3)) + T24*k(4) + ci*T34*k(4)) + 
     &       fo(1)*(-(T12*k(2)) + T24*k(2) - T13*k(3)
     & + T34*k(3) - T14*(k(1) + k(4)))) + 
     &    (fmass*fo(2)*r2**2*Tii - fo(3)*r2*Tii*(k(2) 
     &- ci*k(3)) - fo(4)*r2*Tii*(k(1) - k(4)))*(-p(2) - ci*p(3)) + 
     &    (fo(3)*(T12*k(1) - ci*T13*k(1) + T23*(ci*k(2)
     & - k(3)) - T24*k(4) + ci*T34*k(4)) + 
     &       fo(4)*(-(T12*k(2)) + T24*k(2) - T13*k(3) 
     &+ T34*k(3) - T14*(k(1) + k(4))))*(-p(2) - ci*p(3)) + 
     &    rtwo*(fmass*(fo(2)*(T22*k(2) + ci*T33*k(3)) 
     &+ fo(1)*(T11*k(1) + T44*k(4))) + 
     &       (fo(3)*(-(T22*k(2)) + ci*T33*k(3)) 
     &+ fo(4)*(T11*k(1) + T44*k(4)))*(-p(2) - ci*p(3)) + 
     &       (fo(4)*(-(T22*k(2)) - ci*T33*k(3))
     & + fo(3)*(T11*k(1) - T44*k(4)))*(p(1) - p(4))) + 
     &    (fmass*fo(1)*r2**2*Tii - fo(4)*r2*Tii*(k(2) 
     &+ ci*k(3)) - fo(3)*r2*Tii*(k(1) + k(4)))*(p(1) - p(4)) + 
     &    (fo(4)*(T12*k(1) + ci*T13*k(1) +
     & T23*(-(ci*k(2)) - k(3)) - T24*k(4) - ci*T34*k(4)) + 
     &       fo(3)*(-(T12*k(2)) - T24*k(2) - T13*k(3) 
     &- T34*k(3) - T14*(-k(1) + k(4))))*(p(1) - p(4))

      fto(4) =fmass*(fmass*fo(4)*r2**2*Tii
     &	 - fo(1)*r2*Tii*(-k(2) + ci*k(3))
     & - fo(2)*r2*Tii*(k(1) + k(4))) + 
     &    fmass*(fo(1)*(-(T12*k(1)) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) + k(3)) + T24*k(4) - ci*T34*k(4)) + 
     &       fo(2)*(-(T12*k(2)) - T24*k(2) - T13*k(3) 
     &- T34*k(3) - T14*(-k(1) + k(4)))) + 
     &    (fmass*fo(1)*r2**2*Tii - fo(4)*r2*Tii*(k(2) 
     &+ ci*k(3)) - fo(3)*r2*Tii*(k(1) + k(4)))*(-p(2) + ci*p(3)) + 
     &    (fo(4)*(T12*k(1) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) - k(3)) - T24*k(4) - ci*T34*k(4)) + 
     &       fo(3)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4))))*(-p(2) + ci*p(3)) + 
     &    (fmass*fo(2)*r2**2*Tii - fo(3)*r2*Tii*(k(2) 
     &- ci*k(3)) - fo(4)*r2*Tii*(k(1) - k(4)))*(p(1) + p(4)) + 
     &    (fo(3)*(T12*k(1) - ci*T13*k(1) + T23*(ci*k(2)
     & - k(3)) - T24*k(4) + ci*T34*k(4)) + 
     &       fo(4)*(-(T12*k(2)) + T24*k(2) - T13*k(3) 
     &+ T34*k(3) - T14*(k(1) + k(4))))*(p(1) + p(4)) + 
     &    rtwo*(fmass*(fo(1)*(T22*k(2) - ci*T33*k(3)) 
     &+ fo(2)*(T11*k(1) - T44*k(4))) + 
     &       (fo(4)*(-(T22*k(2)) - ci*T33*k(3))
     & + fo(3)*(T11*k(1) - T44*k(4)))*(-p(2) + ci*p(3)) + 
     &       (fo(3)*(-(T22*k(2)) + ci*T33*k(3))
     & + fo(4)*(T11*k(1) + T44*k(4)))*(p(1) + p(4)))
     
      do i = 1,4
         fto(i) = fto(i)*d
      end do

      return
      end

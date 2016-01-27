      subroutine fvtoxx(fo,vc,ft,g,g2,fmass,fwidth , fvto)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-out external fermion, a gaguge boson and a KK tensor boson.
c
c input:
c       complex fo(6)          : flow-out  fermion                   <fo|
c       complex ft(18)         : input    tensor                      t
c       complex vc(6)          : vector                               v
c       complex g(2)*g2        : coupling constant            -g(L,R)*kappa/4
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvto(6)         : off-shell fermion            <f',t,v,fo|
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex g(2),g2,fo(6),vc(6),ft(18),tc(6,4),fvto(6)
      double precision  fmass, fwidth

      double complex d, T00, T12, T13, T14, T23, T24, T34
      double precision pi(4), pi2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ))


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

      fvto(5) = fo(5) + tc(5,1) +vc(5)
      fvto(6) = fo(6) + tc(6,1) +vc(6)

      pi(1) = dreal(fvto(5))
      pi(2) = dreal(fvto(6))
      pi(3) = dimag(fvto(6))
      pi(4) = dimag(fvto(5))
 
      pi2 = pi(1)**2-pi(2)**2-pi(3)**2-pi(4)**2
     

      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      if ( fmass.gt.rZero ) then
         d =  -1.0d0/dcmplx( pi2-fmass**2, fmass*fwidth )
      else
         d =  -1.0d0/dcmplx( pi2, rZero )
      end if

      fvto(1) = g(2)*(pi(2) + cone*pi(3))*(fo(1)*(T12*vc(1)
     &	 - cone*T13*vc(1) - T23*(-(cone*vc(2)) + vc(3)) + 
     &          2.d0*T00*(-vc(2) + cone*vc(3)) - T24*vc(4)
     & + cone*T34*vc(4)) + 
     &       fo(2)*(T12*vc(2) + T24*vc(2) + T13*vc(3) 
     &+ T34*vc(3) + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))))  
     &   + fmass*g(1)*(fo(4)*(-(T12*vc(1)) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) - vc(3)) + 2.d0*T00*(vc(2) + cone*vc(3)) + 
     &          T24*vc(4) + cone*T34*vc(4)) + fo(3)*
     &        (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    g(2)*(pi(1) + pi(4))*(fo(2)*(T12*vc(1) + cone*T13*vc(1)
     &- T23*(cone*vc(2) + vc(3)) + 2.d0*T00*(-vc(2) - cone*vc(3)) - 
     &          T24*vc(4) - cone*T34*vc(4)) + fo(1)*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3) 
     &+ 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))


      fvto(2) = g(2)*(pi(1) - pi(4))*(fo(1)*(T12*vc(1) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) + vc(3)) + 2.d0*T00*(-vc(2) + cone*vc(3)) 
     &        -   T24*vc(4) + cone*T34*vc(4)) + fo(2)*
     &        (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3) 
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    g(2)*(pi(2) - cone*pi(3))*(fo(2)*(T12*vc(1) 
     & + cone*T13*vc(1) - T23*(cone*vc(2) + vc(3)) + 
     &          2.d0*T00*(-vc(2) - cone*vc(3)) - T24*vc(4)
     & - cone*T34*vc(4)) + 
     &       fo(1)*(T12*vc(2) - T24*vc(2) + T13*vc(3) 
     & - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))  
     &    +fmass*g(1)*(fo(3)*(-(T12*vc(1)) + cone*T13*vc(1)
     & - T23*(cone*vc(2) - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &          T24*vc(4) - cone*T34*vc(4)) + fo(4)*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3)
     & + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))



      fvto(3) = g(1)*(pi(1) - pi(4))*(fo(4)*(-(T12*vc(1))
     &	 - cone*T13*vc(1) - T23*(-(cone*vc(2)) - vc(3)) + 
     &    2.d0*T00*(vc(2) + cone*vc(3)) + T24*vc(4) + cone*T34*vc(4))  
     &       +fo(3)*(T12*vc(2) + T24*vc(2) + T13*vc(3) 
     & + T34*vc(3) + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))))  
     &   + fmass*g(2)*(fo(2)*(T12*vc(1) + cone*T13*vc(1) 
     &- T23*(cone*vc(2) + vc(3)) + 2.d0*T00*(-vc(2) - cone*vc(3))
     & - T24*vc(4) - 
     &          cone*T34*vc(4)) + fo(1)*(T12*vc(2) 
     & - T24*vc(2) + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4))  
     &         + T14*(vc(1) + vc(4)))) + g(1)*(-pi(2) - cone*pi(3))*
     &     (fo(3)*(-(T12*vc(1)) + cone*T13*vc(1) 
     &- T23*(cone*vc(2) - vc(3))
     & + 2.d0*T00*(vc(2) - cone*vc(3)) + T24*vc(4) - 
     &          cone*T34*vc(4)) + fo(4)*(T12*vc(2) - T24*vc(2)
     & + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + 
     &          T14*(vc(1) + vc(4))))

      fvto(4) = fmass*g(2)*(fo(1)*
     &        (T12*vc(1) - cone*T13*vc(1) - T23*(-(cone*vc(2)) + vc(3))
     & + 2.d0*T00*(-vc(2) + cone*vc(3)) - T24*vc(4) + 
     &          cone*T34*vc(4)) + fo(2)*(T12*vc(2) + T24*vc(2) 
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) + 
     &   2.d0*T00*(vc(1) + vc(4)))) + g(1)*(-pi(2) + cone*pi(3))*
     &     (fo(4)*(-(T12*vc(1)) - cone*T13*vc(1) - T23*(-(cone*vc(2))
     & - vc(3)) + 2.d0*T00*(vc(2) + cone*vc(3)) + T24*vc(4) + 
     &          cone*T34*vc(4)) + fo(3)*(T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) + 
     &          2.d0*T00*(vc(1) + vc(4)))) + g(1)*(pi(1) + pi(4))*
     &     (fo(3)*(-(T12*vc(1)) + cone*T13*vc(1) - T23*(cone*vc(2) 
     & - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + T24*vc(4) - 
     &          cone*T34*vc(4)) + fo(4)*(T12*vc(2) - T24*vc(2)
     & + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + 
     &          T14*(vc(1) + vc(4))))

      fvto(1) = fvto(1)+
     &g(2)*(pi(2) + cone*pi(3))*(-2*fo(1)*(tc(2,2)*vc(2)
     & - cone*tc(3,3)*vc(3)) - 
     &       2*fo(2)*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    fmass*g(1)*(-2*fo(4)*(-(tc(2,2)*vc(2)) - cone*tc(3,3)*vc(3))
     & - 2*fo(3)*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    g(2)*(pi(1) + pi(4))*(-2*fo(2)*(tc(2,2)*vc(2) 
     &+ cone*tc(3,3)*vc(3)) - 2*fo(1)*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))

      fvto(2) = fvto(2)+
     &g(2)*(pi(1) - pi(4))*(-2*fo(1)*(tc(2,2)*vc(2)-cone*tc(3,3)*vc(3))
     & - 2*fo(2)*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    g(2)*(pi(2)-cone*pi(3))*(-2*fo(2)*(tc(2,2)*vc(2)
     & + cone*tc(3,3)*vc(3))-2*fo(1)*(tc(1,1)*vc(1) + tc(4,4)*vc(4))) + 
     &    fmass*g(1)*(-2*fo(3)*(-(tc(2,2)*vc(2)) + cone*tc(3,3)*vc(3))
     & - 2*fo(4)*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))
     	
      fvto(3) = fvto(3)+
     &g(1)*(pi(1) - pi(4))*(-2*fo(4)*(-(tc(2,2)*vc(2))
     & - cone*tc(3,3)*vc(3)) - 2*fo(3)*(tc(1,1)*vc(1)-tc(4,4)*vc(4)))+ 
     &    fmass*g(2)*(-2*fo(2)*(tc(2,2)*vc(2) + cone*tc(3,3)*vc(3))
     & - 2*fo(1)*(tc(1,1)*vc(1) + tc(4,4)*vc(4))) + 
     &    g(1)*(-pi(2) - cone*pi(3))*(-2*fo(3)*(-(tc(2,2)*vc(2))
     & + cone*tc(3,3)*vc(3)) - 2*fo(4)*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))
	
      fvto(4) = fvto(4)+
     &fmass*g(2)*(-2*fo(1)*(tc(2,2)*vc(2) - cone*tc(3,3)*vc(3))
     & - 2*fo(2)*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    g(1)*(-pi(2) + cone*pi(3))*(-2*fo(4)*(-(tc(2,2)*vc(2))
     & - cone*tc(3,3)*vc(3)) - 2*fo(3)*(tc(1,1)*vc(1)-tc(4,4)*vc(4)))+ 
     &    g(1)*(pi(1) + pi(4))*(-2*fo(3)*(-(tc(2,2)*vc(2)) 
     &+ cone*tc(3,3)*vc(3)) - 2*fo(4)*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))

      do i = 1,4
         fvto(i) = -fvto(i)*d*g2
      end do

      return
      end

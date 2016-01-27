      subroutine fvtixx(fi,vc,ft,g,g2,fmass,fwidth , fvti)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion, a gaguge boson and a KK tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ft(18)         : input    tensor                      t
c       complex vc(6)          : vector                               v
c       complex g(2)*g2        : coupling constant            -g(L,R)*kappa/4
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvti(6)         : off-shell fermion             |f',t,v,fi>
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex g(2),g2,fi(6),vc(6),ft(18),tc(6,4), fvti(6)
      double precision  fmass, fwidth

      double complex d, T00, T12, T13, T14, T23, T24, T34
      double precision po(4), po2
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

      fvti(5) = fi(5) - tc(5,1) -vc(5)
      fvti(6) = fi(6) - tc(6,1) -vc(6)

      po(1) = dreal(fvti(5))
      po(2) = dreal(fvti(6))
      po(3) = dimag(fvti(6))
      po(4) = dimag(fvti(5))
 
      po2=po(1)**2-po(2)**2-po(3)**2-po(4)**2
     
      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)
   
      if ( fmass.gt.rZero ) then
         d =  -1.0d0/dcmplx( po2-fmass**2, fmass*fwidth )
      else
         d =  -1.0d0/dcmplx( po2, rZero )
      end if

      fvti(1) =  fmass*fi(4)*g(2)*(T12*vc(1) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) + vc(3)) + 2.d0*T00*(-vc(2) + cone*vc(3))
     & - T24*vc(4) + cone*T34*vc(4)) + fmass*fi(3)*g(2)*
     &     (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3)
     &+ 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))) + 
     &    fi(1)*g(1)*((-po(2) + cone*po(3))*(-(T12*vc(1))
     & - cone*T13*vc(1) - T23*(-(cone*vc(2)) - vc(3)) + 
     &          2.d0*T00*(vc(2) + cone*vc(3)) + T24*vc(4) 
     &+ cone*T34*vc(4)) + 
     &       (po(1) - po(4))*(T12*vc(2) + T24*vc(2) 
     &+ T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) 
     &+ 2.d0*T00*(vc(1) + vc(4)))) + 
     &    fi(2)*g(1)*((po(1) - po(4))*(-(T12*vc(1))
     & + cone*T13*vc(1) - T23*(cone*vc(2) - vc(3)) 
     &+ 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &          T24*vc(4) - cone*T34*vc(4)) 
     &+ (-po(2) + cone*po(3))*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) 
     &- T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))

      fvti(2) = fmass*fi(3)*g(2)*(T12*vc(1) 
     &+ cone*T13*vc(1) - T23*(cone*vc(2) + vc(3))
     & + 2.d0*T00*(-vc(2) - cone*vc(3)) - 
     &       T24*vc(4) - cone*T34*vc(4)) + fmass*fi(4)*g(2)*
     &     (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))) + 
     &    fi(1)*g(1)*((po(1) + po(4))*(-(T12*vc(1)) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) - vc(3)) + 
     &          2.d0*T00*(vc(2) + cone*vc(3))
     & + T24*vc(4) + cone*T34*vc(4)) + 
     &       (-po(2) - cone*po(3))*(T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) + 
     &          2.d0*T00*(vc(1) + vc(4)))) + fi(2)*g(1)*
     &     ((-po(2) - cone*po(3))*(-(T12*vc(1)) + cone*T13*vc(1)
     & - T23*(cone*vc(2) - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &          T24*vc(4) - cone*T34*vc(4)) + (po(1) + po(4))*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3) 
     &+ 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))

      fvti(3) = fmass*fi(2)*g(1)*(-(T12*vc(1)) + cone*T13*vc(1)
     & - T23*(cone*vc(2) - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &       T24*vc(4) - cone*T34*vc(4)) + fmass*fi(1)*g(1)*
     &     (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))) + 
     &    fi(4)*g(2)*((po(1) + po(4))*(T12*vc(1) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) + vc(3)) + 2.d0*T00*(-vc(2) + cone*vc(3))  
     &       -T24*vc(4) + cone*T34*vc(4)) + (po(2) - cone*po(3))*
     &        (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    fi(3)*g(2)*((po(2) - cone*po(3))*(T12*vc(1)
     & + cone*T13*vc(1) - T23*(cone*vc(2) + vc(3)) + 
     &          2.d0*T00*(-vc(2) - cone*vc(3))
     & - T24*vc(4) - cone*T34*vc(4)) + 
     &       (po(1) + po(4))*(T12*vc(2) - T24*vc(2)
     & + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) 
     & + T14*(vc(1) + vc(4))))

      fvti(4) = fmass*fi(1)*g(1)*(-(T12*vc(1)) - cone*T13*vc(1)
     &	 - T23*(-(cone*vc(2)) - vc(3)) 
     & + 2.d0*T00*(vc(2) + cone*vc(3)) + 
     &       T24*vc(4) + cone*T34*vc(4)) + fmass*fi(2)*g(1)*
     &     (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3)
     & + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))) + 
     &    fi(4)*g(2)*((po(2) + cone*po(3))*(T12*vc(1) 
     & - cone*T13*vc(1) - T23*(-(cone*vc(2)) + vc(3)) + 
     &          2.d0*T00*(-vc(2) + cone*vc(3)) - T24*vc(4)
     & + cone*T34*vc(4)) + 
     &       (po(1) - po(4))*(T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4))
     & + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    fi(3)*g(2)*((po(1) - po(4))*(T12*vc(1) 
     &+ cone*T13*vc(1) - T23*(cone*vc(2) + vc(3))
     & + 2.d0*T00*(-vc(2) - cone*vc(3)) - 
     &          T24*vc(4) - cone*T34*vc(4))
     & + (po(2) + cone*po(3))*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) 
     &- T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))  

      fvti(1) = fvti(1)
     &	-2*fmass*fi(4)*g(2)*(tc(2,2)*vc(2) - cone*tc(3,3)*vc(3))
     & - 2*fmass*fi(3)*g(2)*(tc(1,1)*vc(1) + tc(4,4)*vc(4)) + 
     &    fi(1)*g(1)*(-2*(-po(2) + cone*po(3))*(-(tc(2,2)*vc(2))
     & - cone*tc(3,3)*vc(3)) - 
     &       2*(po(1) - po(4))*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    fi(2)*g(1)*(-2*(po(1) - po(4))*(-(tc(2,2)*vc(2)) 
     &+ cone*tc(3,3)*vc(3)) - 
     &       2*(-po(2) + cone*po(3))*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))

      fvti(2) = fvti(2)
     &-2*fmass*fi(3)*g(2)*(tc(2,2)*vc(2) + cone*tc(3,3)*vc(3))
     & - 2*fmass*fi(4)*g(2)*(tc(1,1)*vc(1) - tc(4,4)*vc(4)) + 
     &   fi(1)*g(1)*(-2*(po(1) + po(4))*(-(tc(2,2)*vc(2))
     & - cone*tc(3,3)*vc(3)) - 
     &       2*(-po(2) - cone*po(3))*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    fi(2)*g(1)*(-2*(-po(2) - cone*po(3))*(-(tc(2,2)*vc(2))
     & + cone*tc(3,3)*vc(3)) - 
     &       2*(po(1) + po(4))*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))
     	
      fvti(3) = fvti(3)
     &-2*fmass*fi(2)*g(1)*(-(tc(2,2)*vc(2)) + cone*tc(3,3)*vc(3))
     & - 2*fmass*fi(1)*g(1)*(tc(1,1)*vc(1) - tc(4,4)*vc(4)) + 
     &    fi(4)*g(2)*(-2*(po(1) + po(4))*(tc(2,2)*vc(2) 
     &- cone*tc(3,3)*vc(3)) - 
     &       2*(po(2) - cone*po(3))*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    fi(3)*g(2)*(-2*(po(2) - cone*po(3))*(tc(2,2)*vc(2) 
     &+ cone*tc(3,3)*vc(3)) - 
     &       2*(po(1) + po(4))*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))	
	
      fvti(4) = fvti(4)
     &-2*fmass*fi(1)*g(1)*(-(tc(2,2)*vc(2)) - cone*tc(3,3)*vc(3))
     & - 2*fmass*fi(2)*g(1)*(tc(1,1)*vc(1) + tc(4,4)*vc(4)) + 
     &    fi(4)*g(2)*(-2*(po(2) + cone*po(3))*(tc(2,2)*vc(2)
     & - cone*tc(3,3)*vc(3)) - 
     &       2*(po(1) - po(4))*(tc(1,1)*vc(1) - tc(4,4)*vc(4))) + 
     &    fi(3)*g(2)*(-2*(po(1) - po(4))*(tc(2,2)*vc(2)
     & + cone*tc(3,3)*vc(3)) - 
     &       2*(po(2) + cone*po(3))*(tc(1,1)*vc(1) + tc(4,4)*vc(4)))	

      do i = 1,4
         fvti(i) = -fvti(i)*d*g2
      end do

      return
      end

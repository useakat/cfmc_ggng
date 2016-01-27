      subroutine iovtxx(fi,fo,vc,ft,g,g1 , vertex)
c
c This subroutine computes an amplitude of the four-point coupling of
c a vector boson, two fermions and a Kaluza-Klein tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion     SM |fi>
c       complex fo(6)          : flow-out fermion     SM <fo|
c       complex vc(6)          : vector               SM   v
c       complex ft(18)         : tensor               KK   t
c       complex g(2)           : coupling constant    g(L,R)
c       compex  g1                                    -kappa/4
c
c output:
c       complex vertex         : amplitude            gamma(fi,fo,vc,tc)
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6),fo(6),vc(6),ft(18),tc(6,4), vertex
      double complex g(2)
      double complex g1

      double complex T00,T12, T13, T14, T23, T24, T34

      double precision rZero, r2
      parameter( rZero = 0.0d0, r2 = 2.0d0 )
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ) )

      
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

      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)


      vertex =-g1*(fi(4)*g(2)*(fo(1)*(T12*vc(1) - cone*T13*vc(1)
     &	 - r2*tc(2,2)*vc(2) + cone*r2*tc(3,3)*vc(3) - 
     &          T23*(-(cone*vc(2)) + vc(3)) + r2*T00*(-vc(2) 
     &+ cone*vc(3)) - T24*vc(4) + cone*T34*vc(4)) + 
     &       fo(2)*(-(r2*tc(1,1)*vc(1)) + T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + r2*tc(4,4)*vc(4) + 
     &          T14*(-vc(1) + vc(4)) + r2*T00*(vc(1) + vc(4)))) + 
     &    fi(1)*g(1)*(fo(4)*(-(T12*vc(1)) - cone*T13*vc(1) 
     &+ r2*tc(2,2)*vc(2) - T23*(-(cone*vc(2)) - vc(3)) + 
     &          cone*r2*tc(3,3)*vc(3) + r2*T00*(vc(2) 
     &+ cone*vc(3)) + T24*vc(4) + cone*T34*vc(4)) + 
     &       fo(3)*(-(r2*tc(1,1)*vc(1)) + T12*vc(2) 
     &+ T24*vc(2) + T13*vc(3) + T34*vc(3) + r2*tc(4,4)*vc(4) + 
     &          T14*(-vc(1) + vc(4)) + r2*T00*(vc(1) + vc(4)))) + 
     &    fi(3)*g(2)*(fo(2)*(T12*vc(1) + cone*T13*vc(1)
     & - r2*tc(2,2)*vc(2) - cone*r2*tc(3,3)*vc(3)
     & - T23*(cone*vc(2) + vc(3)) + 
     &          r2*T00*(-vc(2) - cone*vc(3)) - T24*vc(4) 
     &- cone*T34*vc(4)) + 
     &       fo(1)*(-(r2*tc(1,1)*vc(1)) + T12*vc(2) 
     &- T24*vc(2) + T13*vc(3) - T34*vc(3) + r2*T00*(vc(1) - vc(4)) - 
     &          r2*tc(4,4)*vc(4) + T14*(vc(1) + vc(4)))) + 
     &    fi(2)*g(1)*(fo(3)*(-(T12*vc(1)) + cone*T13*vc(1)
     & + r2*tc(2,2)*vc(2) - T23*(cone*vc(2) - vc(3)) - 
     &          cone*r2*tc(3,3)*vc(3) + r2*T00*(vc(2) 
     &- cone*vc(3)) + T24*vc(4) - cone*T34*vc(4)) + 
     &       fo(4)*(-(r2*tc(1,1)*vc(1)) + T12*vc(2) 
     &- T24*vc(2) + T13*vc(3) - T34*vc(3) + r2*T00*(vc(1) - vc(4)) - 
     &          r2*tc(4,4)*vc(4) + T14*(vc(1) + vc(4)))))

      return
      end

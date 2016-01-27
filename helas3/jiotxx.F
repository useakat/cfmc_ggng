      subroutine jiotxx(fi,fo,ft,g,g2,vmass,vwidth , jiot)
c
c This subroutine computes an off-shell vector boson wavefunction from a
c flowing-out fermion, a flowing-in fermion and a KK tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                    |fi>
c       complex fo(6)          : flow-out  fermion                   <fo|
c       complex ft(18)         : input    tensor                      t
c       complex g(2)           : coupling constant                  g(L,R)
c       complex g2             : coupling constant                 -kappa/4
c       real    vmass          : mass  of output vector boson
c       real    vwidth         : width of output vector boson
c
c output:
c       complex jiot(6)        : off-shell vector boson           <v,t,f',fo|
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6),fo(6),tc(6,4),ft(18),jiot(6),g(2),g2
      double precision vmass, vwidth

      double complex d, T00, T12, T13, T14, T23, T24, T34
      double precision pv(4), pv2
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
      jiot(5) = fo(5) + tc(5,1) -fi(5)
      jiot(6) = fo(6) + tc(6,1) -fi(6)

      pv(1) = dreal(jiot(5))
      pv(2) = dreal(jiot(6))
      pv(3) = dimag(jiot(6))
      pv(4) = dimag(jiot(5))
 
      pv2=pv(1)**2-pv(2)**2-pv(3)**2-pv(4)**2
     

      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)
      
      if ( vmass.gt.rZero ) then
         d =  1.0d0/dcmplx( pv2-vmass**2, vmass*vwidth )
      else
         d =  1.0d0/dcmplx( pv2, rZero )
      end if

      if (vmass.gt.rzero) then
         jiot(1) = fi(3)*g(2)*(-((rtwo*T00*fo(2)*pv(1)*(-pv(2)
     &	 - cone*pv(3)))/vmass**2) + 
     &       fo(1)*((rtwo*T00) 
     &- (rtwo*T00*pv(1)*(pv(1) - pv(4)))/vmass**2)) + 
     &    fi(2)*g(1)*(-((rtwo*T00*fo(3)*pv(1)
     &*(pv(2) - cone*pv(3)))/vmass**2) + 
     &       fo(4)*((rtwo*T00) - (rtwo*T00*pv(1)
     &*(pv(1) - pv(4)))/vmass**2)) + 
     &    fi(4)*g(2)*(-((rtwo*T00*fo(1)*pv(1)
     &*(-pv(2) + cone*pv(3)))/vmass**2) + 
     &       fo(2)*((rtwo*T00) 
     &- (rtwo*T00*pv(1)*(pv(1) + pv(4)))/vmass**2)) + 
     &    fi(1)*g(1)*(-((rtwo*T00*fo(4)*pv(1)*(pv(2)
     & + cone*pv(3)))/vmass**2) + 
     &       fo(3)*((rtwo*T00) 
     &- (rtwo*T00*pv(1)*(pv(1) + pv(4)))/vmass**2)) + 
     &    fi(4)*g(2)*(fo(1)*(tc(1,2) - cone*tc(1,3) 
     &+ tc(2,1) - cone*tc(3,1))
     & + fo(2)*(-2*tc(1,1) - tc(1,4) - tc(4,1))) + 
     &    fi(1)*g(1)*(fo(4)*(-tc(1,2) - cone*tc(1,3) - tc(2,1) 
     &- cone*tc(3,1)) + fo(3)*(-2*tc(1,1) - tc(1,4) - tc(4,1))) + 
     &    fi(3)*g(2)*(fo(2)*(tc(1,2) + cone*tc(1,3) + tc(2,1) 
     &+ cone*tc(3,1)) + fo(1)*(-2*tc(1,1) + tc(1,4) + tc(4,1))) + 
     &    fi(2)*g(1)*(fo(3)*(-tc(1,2) + cone*tc(1,3) - tc(2,1)
     & + cone*tc(3,1)) + fo(4)*(-2*tc(1,1) + tc(1,4) + tc(4,1))) + 
     &    (pv(1)*(fi(4)*g(2)*(fo(1)*(-(T12*pv(1)) + cone*T13*pv(1)
     & + T23*(-(cone*pv(2)) + pv(3)) + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3) - 
     &T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(1)*g(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3) 
     &- T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(3)*g(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4))) + fi(2)*g(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) - pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3) 
     &+ T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4)))))/vmass**2
      
         jiot(2) = fi(3)*g(2)*(fo(2)*((rtwo*T00) 
     &- (rtwo*T00*pv(2)*(-pv(2) - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(1)*pv(2)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(2)*g(1)*(fo(3)*(-rtwo*T00 - (rtwo*T00*pv(2)*(pv(2)
     & - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(4)*pv(2)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(4)*g(2)*(fo(1)*((rtwo*T00)
     & - (rtwo*T00*pv(2)*(-pv(2) + cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(2)*pv(2)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(1)*g(1)*(fo(4)*(-rtwo*T00 - (rtwo*T00*pv(2)*(pv(2) 
     &+ cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(3)*pv(2)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(4)*g(2)*(fo(1)*(2*tc(2,2) - cone*tc(2,3) - cone*tc(3,2))
     & + fo(2)*(-tc(1,2) - tc(2,1) - tc(2,4) - tc(4,2))) + 
     &    fi(1)*g(1)*(fo(4)*(-2*tc(2,2) - cone*tc(2,3) - cone*tc(3,2)) 
     &+ fo(3)*(-tc(1,2) - tc(2,1) - tc(2,4) - tc(4,2))) + 
     &    fi(3)*g(2)*(fo(2)*(2*tc(2,2) + cone*tc(2,3) + cone*tc(3,2)) 
     &+ fo(1)*(-tc(1,2) - tc(2,1) + tc(2,4) + tc(4,2))) + 
     &    fi(2)*g(1)*(fo(3)*(-2*tc(2,2) + cone*tc(2,3) + cone*tc(3,2))
     & + fo(4)*(-tc(1,2) - tc(2,1) + tc(2,4) + tc(4,2))) + 
     &    (pv(2)*(fi(4)*g(2)*(fo(1)*(-(T12*pv(1)) + cone*T13*pv(1)
     & + T23*(-(cone*pv(2)) + pv(3)) + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(1)*g(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) + T23*(-(cone*pv(2))
     & - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(3)*g(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4))) + fi(2)*g(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) - pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3) 
     &+ T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4)))))/vmass**2

         jiot(3) =fi(3)*g(2)*(fo(2)*((cone*rtwo*T00)
     &	 - (rtwo*T00*pv(3)*(-pv(2) - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(1)*pv(3)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(2)*g(1)*(fo(3)*((cone*rtwo*T00) 
     &- (rtwo*T00*pv(3)*(pv(2) - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(4)*pv(3)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(4)*g(2)*(fo(1)*(-cone*rtwo*T00 
     &- (rtwo*T00*pv(3)*(-pv(2) + cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(2)*pv(3)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(1)*g(1)*(fo(4)*(-cone*rtwo*T00
     & - (rtwo*T00*pv(3)*(pv(2) + cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(3)*pv(3)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(4)*g(2)*(fo(1)*(tc(2,3) + tc(3,2) 
     &- 2*cone*tc(3,3)) + fo(2)*(-tc(1,3) - tc(3,1) 
     &- tc(3,4) - tc(4,3))) + 
     &    fi(1)*g(1)*(fo(4)*(-tc(2,3) - tc(3,2) 
     &- 2*cone*tc(3,3)) + fo(3)*(-tc(1,3) - tc(3,1) 
     &- tc(3,4) - tc(4,3))) + 
     &    fi(3)*g(2)*(fo(2)*(tc(2,3) + tc(3,2) + 2*cone*tc(3,3))
     & + fo(1)*(-tc(1,3) - tc(3,1) + tc(3,4) + tc(4,3))) + 
     &    fi(2)*g(1)*(fo(3)*(-tc(2,3) - tc(3,2) + 2*cone*tc(3,3)) 
     &+ fo(4)*(-tc(1,3) - tc(3,1) + tc(3,4) + tc(4,3))) + 
     &    (pv(3)*(fi(4)*g(2)*(fo(1)*(-(T12*pv(1)) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) + pv(3)) + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     &- T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(1)*g(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(3)*g(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4))) + fi(2)*g(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) - pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4)))))/vmass**2

         jiot(4) =fi(3)*g(2)*(-((rtwo*T00*fo(2)*(-pv(2) 
     &	- cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(1)*((rtwo*T00) - (rtwo*T00*(pv(1)
     & - pv(4))*pv(4))/vmass**2)) + 
     &    fi(2)*g(1)*(-((rtwo*T00*fo(3)*(pv(2)
     & - cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(4)*((rtwo*T00) - (rtwo*T00*(pv(1)
     & - pv(4))*pv(4))/vmass**2)) + 
     &    fi(4)*g(2)*(-((rtwo*T00*fo(1)*(-pv(2) 
     &+ cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(2)*(-rtwo*T00
     & - (rtwo*T00*pv(4)*(pv(1) + pv(4)))/vmass**2)) + 
     &    fi(1)*g(1)*(-((rtwo*T00*fo(4)*(pv(2)
     & + cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(3)*(-rtwo*T00 - (rtwo*T00*pv(4)*(pv(1)
     & + pv(4)))/vmass**2)) + 
     &    fi(4)*g(2)*(fo(1)*(tc(2,4) - cone*tc(3,4)
     & + tc(4,2) - cone*tc(4,3)) + fo(2)*(-tc(1,4) 
     &- tc(4,1) - 2*tc(4,4))) + 
     &    fi(1)*g(1)*(fo(4)*(-tc(2,4) - cone*tc(3,4)
     & - tc(4,2) - cone*tc(4,3)) + fo(3)*(-tc(1,4)
     & - tc(4,1) - 2*tc(4,4))) + 
     &    fi(3)*g(2)*(fo(2)*(tc(2,4) + cone*tc(3,4) 
     &+ tc(4,2) + cone*tc(4,3)) + fo(1)*(-tc(1,4) 
     &- tc(4,1) + 2*tc(4,4))) + 
     &    fi(2)*g(1)*(fo(3)*(-tc(2,4) + cone*tc(3,4)
     & - tc(4,2) + cone*tc(4,3)) + fo(4)*(-tc(1,4)
     & - tc(4,1) + 2*tc(4,4))) + 
     &    (pv(4)*(fi(4)*g(2)*(fo(1)*(-(T12*pv(1)) 
     &+ cone*T13*pv(1) + T23*(-(cone*pv(2)) + pv(3))
     & + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(1)*g(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) - cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) - 
     &               rtwo*pv(4)*tc(4,4))) + fi(3)*g(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3) 
     &+ T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4))) + fi(2)*g(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) + T23*(cone*pv(2) 
     &- pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*tc(2,2) + cone*rtwo*pv(3)*tc(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*tc(1,1) + 
     &               rtwo*pv(4)*tc(4,4)))))/vmass**2

      else

         jiot(1) = rtwo*T00*fi(1)*fo(3)*g(1)
     & + rtwo*T00*fi(2)*fo(4)*g(1) + rtwo*T00*fi(3)*fo(1)*g(2) + 
     &    rtwo*T00*fi(4)*fo(2)*g(2) + fi(4)*g(2)*(fo(1)*(tc(1,2) 
     &- cone*tc(1,3) + tc(2,1) - cone*tc(3,1)) + 
     &       fo(2)*(-2*tc(1,1) - tc(1,4) - tc(4,1))) + 
     &    fi(1)*g(1)*(fo(4)*(-tc(1,2) - cone*tc(1,3) - tc(2,1) 
     &- cone*tc(3,1)) + fo(3)*(-2*tc(1,1) - tc(1,4) - tc(4,1))) + 
     &    fi(3)*g(2)*(fo(2)*(tc(1,2) + cone*tc(1,3) + tc(2,1) 
     &+ cone*tc(3,1)) + fo(1)*(-2*tc(1,1) + tc(1,4) + tc(4,1))) + 
     &    fi(2)*g(1)*(fo(3)*(-tc(1,2) + cone*tc(1,3) - tc(2,1)
     & + cone*tc(3,1)) + fo(4)*(-2*tc(1,1) + tc(1,4) + tc(4,1)))

         jiot(2) = -(rtwo*T00*fi(2)*fo(3)*g(1)) 
     &  - rtwo*T00*fi(1)*fo(4)*g(1)
     &	 + rtwo*T00*fi(4)*fo(1)*g(2) + 
     &    rtwo*T00*fi(3)*fo(2)*g(2) + fi(4)*g(2)*(fo(1)*(2*tc(2,2) 
     &- cone*tc(2,3) - cone*tc(3,2)) + 
     &       fo(2)*(-tc(1,2) - tc(2,1) - tc(2,4) - tc(4,2))) + 
     &    fi(1)*g(1)*(fo(4)*(-2*tc(2,2) - cone*tc(2,3) - cone*tc(3,2))
     & + fo(3)*(-tc(1,2) - tc(2,1) - tc(2,4) - tc(4,2))) + 
     &    fi(3)*g(2)*(fo(2)*(2*tc(2,2) + cone*tc(2,3) + cone*tc(3,2))
     & + fo(1)*(-tc(1,2) - tc(2,1) + tc(2,4) + tc(4,2))) + 
     &    fi(2)*g(1)*(fo(3)*(-2*tc(2,2) + cone*tc(2,3) + cone*tc(3,2)) 
     &+ fo(4)*(-tc(1,2) - tc(2,1) + tc(2,4) + tc(4,2)))

         jiot(3) = cone*rtwo*T00*fi(2)*fo(3)*g(1) 
     &	- cone*rtwo*T00*fi(1)*fo(4)*g(1)
     & - cone*rtwo*T00*fi(4)*fo(1)*g(2) + 
     &    cone*rtwo*T00*fi(3)*fo(2)*g(2) + fi(4)*g(2)*
     &     (fo(1)*(tc(2,3) + tc(3,2) - 2*cone*tc(3,3))
     & + fo(2)*(-tc(1,3) - tc(3,1) - tc(3,4) - tc(4,3))) + 
     &    fi(1)*g(1)*(fo(4)*(-tc(2,3) - tc(3,2) - 2*cone*tc(3,3))
     & + fo(3)*(-tc(1,3) - tc(3,1) - tc(3,4) - tc(4,3))) + 
     &    fi(3)*g(2)*(fo(2)*(tc(2,3) + tc(3,2) + 2*cone*tc(3,3)) 
     &+ fo(1)*(-tc(1,3) - tc(3,1) + tc(3,4) + tc(4,3))) + 
     &    fi(2)*g(1)*(fo(3)*(-tc(2,3) - tc(3,2) + 2*cone*tc(3,3)) 
     &+ fo(4)*(-tc(1,3) - tc(3,1) + tc(3,4) + tc(4,3)))
       
         jiot(4) = -(rtwo*T00*fi(1)*fo(3)*g(1))
     &	 + rtwo*T00*fi(2)*fo(4)*g(1) + rtwo*T00*fi(3)*fo(1)*g(2) - 
     &    rtwo*T00*fi(4)*fo(2)*g(2) + fi(4)*g(2)*(fo(1)*(tc(2,4) 
     &- cone*tc(3,4) + tc(4,2) - cone*tc(4,3)) + 
     &       fo(2)*(-tc(1,4) - tc(4,1) - 2*tc(4,4))) + 
     &    fi(1)*g(1)*(fo(4)*(-tc(2,4) - cone*tc(3,4) - tc(4,2)
     & - cone*tc(4,3)) + fo(3)*(-tc(1,4) - tc(4,1) - 2*tc(4,4))) + 
     &    fi(3)*g(2)*(fo(2)*(tc(2,4) + cone*tc(3,4) + tc(4,2)
     & + cone*tc(4,3)) + fo(1)*(-tc(1,4) - tc(4,1) + 2*tc(4,4))) + 
     &    fi(2)*g(1)*(fo(3)*(-tc(2,4) + cone*tc(3,4) - tc(4,2) 
     &+ cone*tc(4,3)) + fo(4)*(-tc(1,4) - tc(4,1) + 2*tc(4,4)))  

      endif

      do i = 1,4
         jiot(i) = -jiot(i)*d*g2
      end do

      return
      end

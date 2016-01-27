      subroutine uiovxx(fi,fo,vc,g,g2,tmass,twidth , ft)
c
c This subroutine computes an off-shell KK tensor boson wavefunction from a
c flowing-out fermion, a flowing-in fermion and a gauge vector boson.
c
c input:
c       complex fi(6)      : flow-in  fermion                    |fi>
c       complex fo(6)      : flow-out  fermion                   <fo|
c       complex vc(6)      : input  gauge vector boson            v
c       complex g2*g(2)    : coupling constant            -g(L,R)*kappa/4
c       real    tmass      : mass  of output KK tensor boson
c       real    twidth     : width of output KK tensor boson
c
c output:
c       complex ft(18)    : off-shell KK tensor boson          <t, v,fi,fo|
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6), fo(6),vc(6),ft(18),yiov(6,4),g(2),g2
      double precision  tmass, twidth

      double complex d
      double precision pt(4), pt2, KTVC
      double precision MET(4,4)
      integer i,j
      
      double precision rZero, rTwo,r4,r8,r3
      parameter(rZero=0.0d0, rTwo=2.0d0, r4=4.0d0, r8=8.0d0, r3=3.0d0)
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ))


      yiov(5,1) = fo(5) + vc(5) -fi(5)
      yiov(6,1) = fo(6) + vc(6) -fi(6)

      pt(1) = dreal(yiov(5,1))
      pt(2) = dreal(yiov(6,1))
      pt(3) = dimag(yiov(6,1))
      pt(4) = dimag(yiov(5,1))
      
      pt2 = pt(1)**2-pt(2)**2-pt(3)**2-pt(4)**2
      KTVC = pt(1)*vc(1)-pt(2)*vc(2)-pt(3)*vc(3)-pt(4)*vc(4)
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      if ( tmass.gt.rZero ) then
         d =  -1.0d0/dcmplx( pt2-tmass**2, tmass*twidth )
      else
         d =  -1.0d0/dcmplx( pt2, rZero )
      end if

      do i = 1,4
	 do j=1,4
            yiov(i,j) = fi(3)*g(2)*(fo(2)*(-((KTVC*r4*MET(i,j)*(-pt(2)
     &	    - cone*pt(3)))/(r3*Tmass**2)) - 
     &      (KTVC*r8*(-pt(2) - cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) +
     &          (pt2*r4*MET(i,j)*(-vc(2) - cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(-vc(2) - cone*vc(3)))/(r3*Tmass**4)
     & - (r4*pt(i)*pt(j)*(-vc(2) - cone*vc(3)))/Tmass**2 + 
     &    (rtwo*(-pt(2) - cone*pt(3))*(pt(j)*vc(i) + pt(i)*vc(j)))
     &/Tmass**2) + 
     &       fo(1)*(-((KTVC*r4*MET(i,j)*(pt(1) - pt(4)))/(r3*Tmass**2))
     & - (KTVC*r8*(pt(1) - pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) - vc(4)))/(r3*Tmass**2)
     & + (pt2*r8*pt(i)*pt(j)*(vc(1) - vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) - vc(4)))/Tmass**2
     & + (rtwo*(pt(1) - pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2))  
     &   + fi(2)*g(1)*(fo(3)*(-((KTVC*r4*MET(i,j)*(pt(2) 
     &- cone*pt(3)))/(r3*Tmass**2)) - 
     &     (KTVC*r8*(pt(2) - cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(2) - cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(vc(2) - cone*vc(3)))/(r3*Tmass**4)
     & - (r4*pt(i)*pt(j)*(vc(2) - cone*vc(3)))/Tmass**2 + 
     &    (rtwo*(pt(2) - cone*pt(3))*(pt(j)*vc(i) 
     &+ pt(i)*vc(j)))/Tmass**2) + 
     &       fo(4)*(-((KTVC*r4*MET(i,j)*(pt(1) - pt(4)))/(r3*Tmass**2)) 
     &- (KTVC*r8*(pt(1) - pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) - vc(4)))/(r3*Tmass**2) 
     &+ (pt2*r8*pt(i)*pt(j)*(vc(1) - vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) - vc(4)))/Tmass**2 
     &+ (rtwo*(pt(1) - pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2)) + 
     &    fi(4)*g(2)*(fo(1)*(-((KTVC*r4*MET(i,j)*(-pt(2) + cone*pt(3)))
     &/(r3*Tmass**2)) - 
     &      (KTVC*r8*(-pt(2) + cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(-vc(2) + cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(-vc(2) + cone*vc(3)))/(r3*Tmass**4) 
     &- (r4*pt(i)*pt(j)*(-vc(2) + cone*vc(3)))/Tmass**2 + 
     &          (rtwo*(-pt(2) + cone*pt(3))*(pt(j)*vc(i) 
     &+ pt(i)*vc(j)))/Tmass**2) + 
     &       fo(2)*(-((KTVC*r4*MET(i,j)*(pt(1) + pt(4)))/(r3*Tmass**2))
     & - (KTVC*r8*(pt(1) + pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) + vc(4)))/(r3*Tmass**2)
     & + (pt2*r8*pt(i)*pt(j)*(vc(1) + vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) + vc(4)))/Tmass**2
     & + (rtwo*(pt(1) + pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2)) + 
     &    fi(1)*g(1)*(fo(4)*(-((KTVC*r4*MET(i,j)*(pt(2) + cone*pt(3)))
     &/(r3*Tmass**2)) - 
     &     (KTVC*r8*(pt(2) + cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(2) + cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(vc(2) + cone*vc(3)))/(r3*Tmass**4)
     & - (r4*pt(i)*pt(j)*(vc(2) + cone*vc(3)))/Tmass**2 + 
     &          (rtwo*(pt(2) + cone*pt(3))*(pt(j)*vc(i)
     & + pt(i)*vc(j)))/Tmass**2) + 
     &       fo(3)*(-((KTVC*r4*MET(i,j)*(pt(1) + pt(4)))/(r3*Tmass**2))
     & - (KTVC*r8*(pt(1) + pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) + vc(4)))/(r3*Tmass**2) 
     &+ (pt2*r8*pt(i)*pt(j)*(vc(1) + vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) + vc(4)))/Tmass**2 
     &+ (rtwo*(pt(1) + pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2))

         end do
      end do 

      yiov(1,1) = yiov(1,1)
     & +fi(1)*fo(3)*g(1)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1)) + 
     &  fi(2)*fo(4)*g(1)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1)) + 
     &  fi(3)*fo(1)*g(2)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1)) 
     &+ fi(4)*fo(2)*g(2)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1))

      yiov(2,2) = yiov(2,2)
     &	+fi(4)*fo(1)*g(2)*((4*KTVC*pt(2))/Tmass**2 - 2*rtwo*vc(2)) + 
     &    fi(3)*fo(2)*g(2)*((4*KTVC*pt(2))/Tmass**2 - 2*rtwo*vc(2)) + 
     &    fi(2)*fo(3)*g(1)*((-4*KTVC*pt(2))/Tmass**2 + 2*rtwo*vc(2)) + 
     &    fi(1)*fo(4)*g(1)*((-4*KTVC*pt(2))/Tmass**2 + 2*rtwo*vc(2))

      yiov(3,3) = yiov(3,3)
     &+fi(2)*fo(3)*g(1)*((4*cone*KTVC*pt(3))/Tmass**2
     & - 2*cone*rtwo*vc(3)) + 
     &    fi(3)*fo(2)*g(2)*((4*cone*KTVC*pt(3))/Tmass**2 
     &- 2*cone*rtwo*vc(3)) + 
     &    fi(1)*fo(4)*g(1)*((-4*cone*KTVC*pt(3))/Tmass**2
     & + 2*cone*rtwo*vc(3)) + 
     &    fi(4)*fo(1)*g(2)*((-4*cone*KTVC*pt(3))/Tmass**2
     & + 2*cone*rtwo*vc(3))

      yiov(4,4) = yiov(4,4)
     &+fi(2)*fo(4)*g(1)*((4*KTVC*pt(4))/Tmass**2 - 2*rtwo*vc(4)) + 
     &    fi(3)*fo(1)*g(2)*((4*KTVC*pt(4))/Tmass**2 - 2*rtwo*vc(4)) + 
     &    fi(1)*fo(3)*g(1)*((-4*KTVC*pt(4))/Tmass**2 + 2*rtwo*vc(4)) + 
     &    fi(4)*fo(2)*g(2)*((-4*KTVC*pt(4))/Tmass**2 + 2*rtwo*vc(4))

      yiov(1,2) = yiov(1,2)
     &	+fi(3)*g(2)*(fo(2)*((2*KTVC*pt(1))/Tmass**2 - rtwo*vc(1))
     & + fo(1)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))) + 
     &    fi(4)*g(2)*(fo(1)*((2*KTVC*pt(1))/Tmass**2 - rtwo*vc(1)) 
     &+ fo(2)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))) + 
     &    fi(1)*g(1)*(fo(4)*((-2*KTVC*pt(1))/Tmass**2 + rtwo*vc(1))
     & + fo(3)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))) + 
     &    fi(2)*g(1)*(fo(3)*((-2*KTVC*pt(1))/Tmass**2 + rtwo*vc(1))
     & + fo(4)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2)))

      yiov(1,3) = yiov(1,3)
     &	+fi(3)*g(2)*(fo(2)*((2*cone*KTVC*pt(1))/Tmass**2 
     &- cone*rtwo*vc(1)) + 
     &       fo(1)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3))) + 
     &    fi(4)*g(2)*(fo(1)*((-2*cone*KTVC*pt(1))/Tmass**2 
     &+ cone*rtwo*vc(1)) + fo(2)*((2*KTVC*pt(3))/Tmass**2 
     &- rtwo*vc(3))) + 
     &    fi(1)*g(1)*(fo(4)*((-2*cone*KTVC*pt(1))/Tmass**2
     & + cone*rtwo*vc(1)) + fo(3)*((2*KTVC*pt(3))/Tmass**2 
     &- rtwo*vc(3))) + 
     &    fi(2)*g(1)*(fo(3)*((2*cone*KTVC*pt(1))/Tmass**2
     & -cone*rtwo*vc(1))+ fo(4)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3)))

      yiov(1,4) = yiov(1,4)
     &+fi(1)*fo(3)*g(1)*((2*KTVC*(-pt(1) + pt(4)))/Tmass**2
     & - rtwo*(-vc(1) + vc(4))) + 
     &     fi(4)*fo(2)*g(2)*((2*KTVC*(-pt(1) + pt(4)))/Tmass**2
     & - rtwo*(-vc(1) + vc(4))) + 
     &    fi(2)*fo(4)*g(1)*((2*KTVC*(pt(1) + pt(4)))/Tmass**2 
     &- rtwo*(vc(1) + vc(4))) + 
     &    fi(3)*fo(1)*g(2)*((2*KTVC*(pt(1) + pt(4)))/Tmass**2 
     &- rtwo*(vc(1) + vc(4)))

      yiov(2,3) = yiov(2,3)
     &+fi(1)*fo(4)*g(1)*((2*KTVC*(-(cone*pt(2)) - pt(3)))/Tmass**2
     & - rtwo*(-(cone*vc(2)) - vc(3))) + 
     &    fi(2)*fo(3)*g(1)*((2*KTVC*(cone*pt(2) - pt(3)))/Tmass**2
     & - rtwo*(cone*vc(2) - vc(3))) + 
     &    fi(4)*fo(1)*g(2)*((2*KTVC*(-(cone*pt(2)) + pt(3)))/Tmass**2
     & - rtwo*(-(cone*vc(2)) + vc(3))) + 
     &    fi(3)*fo(2)*g(2)*((2*KTVC*(cone*pt(2) + pt(3)))/Tmass**2 
     &- rtwo*(cone*vc(2) + vc(3)))

      yiov(2,4) = yiov(2,4)
     &+fi(4)*g(2)*(fo(2)*((-2*KTVC*pt(2))/Tmass**2 + rtwo*vc(2))
     & + fo(1)*((2*KTVC*pt(4))/Tmass**2 - rtwo*vc(4))) + 
     &    fi(3)*g(2)*(fo(1)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))
     & + fo(2)*((2*KTVC*pt(4))/Tmass**2 - rtwo*vc(4))) + 
     &    fi(2)*g(1)*(fo(4)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2)) 
     &+ fo(3)*((-2*KTVC*pt(4))/Tmass**2 + rtwo*vc(4))) + 
     &    fi(1)*g(1)*(fo(3)*((-2*KTVC*pt(2))/Tmass**2 + rtwo*vc(2))
     & + fo(4)*((-2*KTVC*pt(4))/Tmass**2 + rtwo*vc(4)))

      yiov(3,4) = yiov(3,4)
     &	+fi(3)*g(2)*(fo(1)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3)) + 
     &       fo(2)*((2*cone*KTVC*pt(4))/Tmass**2 - cone*rtwo*vc(4))) + 
     &    fi(2)*g(1)*(fo(4)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3)) 
     &+ fo(3)*((2*cone*KTVC*pt(4))/Tmass**2 - cone*rtwo*vc(4))) + 
     &    fi(4)*g(2)*(fo(2)*((-2*KTVC*pt(3))/Tmass**2 + rtwo*vc(3)) + 
     &       fo(1)*((-2*cone*KTVC*pt(4))/Tmass**2 + cone*rtwo*vc(4))) + 
     &    fi(1)*g(1)*(fo(3)*((-2*KTVC*pt(3))/Tmass**2 + rtwo*vc(3)) 
     &+ fo(4)*((-2*cone*KTVC*pt(4))/Tmass**2 + cone*rtwo*vc(4)))

      yiov(2,1) = yiov(1,2)
      yiov(3,1) = yiov(1,3)
      yiov(4,1) = yiov(1,4)
      yiov(3,2) = yiov(2,3)
      yiov(4,2) = yiov(2,4)
      yiov(4,3) = yiov(3,4)
			
      do i = 1,4
	 do j=1,4
            yiov(i,j) = -yiov(i,j)*d/2.0d0*g2
         end do
      end do
  
      ft(1) = yiov(1,1)
      ft(2) = yiov(1,2)
      ft(3) = yiov(1,3)
      ft(4) = yiov(1,4)
      ft(5) = yiov(2,1)
      ft(6) = yiov(2,2)
      ft(7) = yiov(2,3)
      ft(8) = yiov(2,4)
      ft(9) = yiov(3,1)
      ft(10) = yiov(3,2)
      ft(11) = yiov(3,3)
      ft(12) = yiov(3,4)
      ft(13) = yiov(4,1)
      ft(14) = yiov(4,2)
      ft(15) = yiov(4,3)
      ft(16) = yiov(4,4)
      ft(17) = yiov(5,1)
      ft(18) = yiov(6,1)

      return
      end

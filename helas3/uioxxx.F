      subroutine uioxxx(fi,fo,g,fmass,Tmass,Twidth , ft)
c
c This subroutine computes an off-shell KK Tensor wavfunction from
c the fermion-anti-fermion-Kaluza-Klein tensor boson coupling.
c
c input:
c       complex fi(6)         : flow-in fermion               |fi>
c       complex fo(6)         : flow-out fermion              <fo|
c       real    g             : coupling constant            -kappa/8
c       real    fmass         : mass  of input fermion
c       real    Tmass         : mass  of OUTPUT KK Tensor T'
c       real    Twidth        : width of OUTPUT KK Tensor T'
c output:
c       complex ft(18)        : KK tensor current       j(T:fi,fo)     
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex fi(6),fo(6),yio(6,4),ft(18)
      double precision g,fmass,Tmass,Twidth
      integer i,j
      double precision pi(4),po(4),km(4),kp(4)
      double precision MET(4,4)
      double complex cone,cZero, d, tt1,tt2,tt3
      double precision rZero, rTwo
      double precision KT2,K1KT,K2KT
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ), cone=(0.0d0,1.0d0) )

      
      yio(5,1) = -fi(5)+fo(5)
      yio(6,1) = -fi(6)+fo(6)
      
      pi(1) = dreal(fi(5))
      pi(2) = dreal(fi(6))
      pi(3) = dimag(fi(6))
      pi(4) = dimag(fi(5))
      
      po(1) = dreal(fo(5))
      po(2) = dreal(fo(6))
      po(3) = dimag(fo(6))
      po(4) = dimag(fo(5))
      
      km(1) = dreal(yio(5,1))
      km(2) = dreal(yio(6,1))
      km(3) = dimag(yio(6,1))
      km(4) = dimag(yio(5,1))

      kp(1) = po(1)+pi(1)
      kp(2) = po(2)+pi(2)
      kp(3) = po(3)+pi(3)
      kp(4) = po(4)+pi(4)
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      K1KT = pi(1)*km(1)-pi(2)*km(2)-pi(3)*km(3)-pi(4)*km(4)
      K2KT = po(1)*km(1)-po(2)*km(2)-po(3)*km(3)-po(4)*km(4)
      KT2 = km(1)**2-km(2)**2-km(3)**2-km(4)**2
      
      if ( Tmass.gt.rZero ) then
         d = - g/dcmplx( KT2-Tmass**2, Tmass*Twidth )
      else
         d = - g/dcmplx( KT2, rZero )
      end if

      tt1 = fi(3)*(((-8.0d0*fmass)/3.d0 
     &	+ (8.0d0*fmass*KT2)/(3.d0*Tmass**2))*fo(3)
     & +  fo(2)*((4.0d0*(K1KT + K2KT)
     &            *(-km(2) - cone*km(3)))/(3.d0*Tmass**2)
     & - (4.0d0*KT2*(-kp(2) - cone*kp(3)))/(3.d0*Tmass**2)) + 
     &    fo(1)*((4.0d0*(K1KT + K2KT)*(km(1) - km(4)))/(3.d0*Tmass**2) 
     & - (4.0d0*KT2*(kp(1) - kp(4)))/(3.d0*Tmass**2))) + 
     &    fi(2)*(((-8.0d0*fmass)/3.d0 
     &        + (8.0d0*fmass*KT2)/(3.d0*Tmass**2))*fo(2)  
     & +  fo(3)*((4.0d0*(K1KT + K2KT)*(km(2) 
     &          - cone*km(3)))/(3.d0*Tmass**2) 
     &- (4.0d0*KT2*(kp(2) - cone*kp(3)))/(3.d0*Tmass**2)) + 
     &    fo(4)*((4.0d0*(K1KT + K2KT)*(km(1) - km(4)))/(3.d0*Tmass**2)
     & - (4.0d0*KT2*(kp(1) - kp(4)))/(3.d0*Tmass**2))) + 
     &    fi(4)*(((-8.0d0*fmass)/3.d0 
     &+ (8.0d0*fmass*KT2)/(3.d0*Tmass**2))*fo(4)  
     & +  fo(1)*((4.0d0*(K1KT + K2KT)
     &       *(-km(2) + cone*km(3)))/(3.d0*Tmass**2)
     & - (4.0d0*KT2*(-kp(2) + cone*kp(3)))/(3.d0*Tmass**2)) + 
     &    fo(2)*((4.0d0*(K1KT + K2KT)*(km(1) + km(4)))/(3.d0*Tmass**2)
     & - (4.0d0*KT2*(kp(1) + kp(4)))/(3.d0*Tmass**2))) + 
     &    fi(1)*(((-8.0d0*fmass)/3.d0 
     &+ (8.0d0*fmass*KT2)/(3.d0*Tmass**2))*fo(1) 
     & +  fo(4)*((4.0d0*(K1KT + K2KT)*(km(2)
     &      + cone*km(3)))/(3.d0*Tmass**2) 
     &- (4.0d0*KT2*(kp(2) + cone*kp(3)))/(3.d0*Tmass**2)) + 
     &    fo(3)*((4.0d0*(K1KT + K2KT)
     &*(km(1) + km(4)))/(3.d0*Tmass**2)
     & - (4.0d0*KT2*(kp(1) + kp(4)))/(3.d0*Tmass**2)))
 
      
      tt2 = (fi(3)*(((-16.0d0*KT2*fmass)/(3.0D0*Tmass**4)
     &	 - (16.d0*fmass)/(3.d0*Tmass**2))*fo(3) + 
     &       fo(2)*((8.d0*(K1KT + K2KT)*(-km(2)
     & - cone*km(3)))/(3.d0*Tmass**4) 
     &+ (8.d0*KT2*(-kp(2) - cone*kp(3)))/(3.0D0*Tmass**4) + 
     &          (4.d0*(-kp(2) - cone*kp(3)))/Tmass**2) + 
     &       fo(1)*((8.d0*(K1KT + K2KT)*(km(1)
     & - km(4)))/(3.d0*Tmass**4) 
     &+ (8.d0*KT2*(kp(1) - kp(4)))/(3.d0*Tmass**4) + 
     &          (4.d0*(kp(1) - kp(4)))/Tmass**2)) + fi(2)*
     &     (((-16.d0*KT2*fmass)/(3.d0*Tmass**4) 
     &- (16.d0*fmass)/(3.d0*Tmass**2))*fo(2) + 
     &       fo(3)*((8.d0*(K1KT + K2KT)*(km(2)
     & - cone*km(3)))/(3.d0*Tmass**4) 
     &+ (8.d0*KT2*(kp(2) - cone*kp(3)))/(3.d0*Tmass**4) + 
     &          (4.d0*(kp(2) - cone*kp(3)))/Tmass**2) + 
     &       fo(4)*((8.d0*(K1KT + K2KT)*(km(1) - km(4)))/(3.d0*Tmass**4)
     & + (8.d0*KT2*(kp(1) - kp(4)))/(3.d0*Tmass**4) + 
     &          (4.d0*(kp(1) - kp(4)))/Tmass**2)) + fi(4)*
     &     (((-16.d0*KT2*fmass)/(3.d0*Tmass**4) 
     &- (16.d0*fmass)/(3.d0*Tmass**2))*fo(4) + 
     &       fo(1)*((8.d0*(K1KT + K2KT)*(-km(2) 
     &+ cone*km(3)))/(3.d0*Tmass**4)
     & + (8.d0*KT2*(-kp(2) + cone*kp(3)))/(3.d0*Tmass**4) + 
     &          (4.d0*(-kp(2) + cone*kp(3)))/Tmass**2) + 
     &       fo(2)*((8.d0*(K1KT + K2KT)*(km(1) 
     &+ km(4)))/(3.d0*Tmass**4)
     & + (8.d0*KT2*(kp(1) + kp(4)))/(3.d0*Tmass**4) + 
     &          (4.d0*(kp(1) + kp(4)))/Tmass**2)) + fi(1)*
     &     (((-16.d0*KT2*fmass)/(3.d0*Tmass**4) 
     &- (16.d0*fmass)/(3.d0*Tmass**2))*fo(1) + 
     &       fo(4)*((8.d0*(K1KT + K2KT)*(km(2)
     & + cone*km(3)))/(3.d0*Tmass**4)
     & + (8.d0*KT2*(kp(2) + cone*kp(3)))/(3.d0*Tmass**4) + 
     &          (4.d0*(kp(2) + cone*kp(3)))/Tmass**2) + 
     &       fo(3)*((8.d0*(K1KT + K2KT)
     &*(km(1) + km(4)))/(3.d0*Tmass**4)
     & + (8.d0*KT2*(kp(1) + kp(4)))/(3.d0*Tmass**4) + 
     &          (4.d0*(kp(1) + kp(4)))/Tmass**2)))

      tt3 = fi(3)*((-2.0d0*fo(2)*(-km(2) - cone*km(3)))/Tmass**2
     &	 - (2.0d0*fo(1)*(km(1) - km(4)))/Tmass**2) + 
     &    fi(2)*((-2.0d0*fo(3)*(km(2) - cone*km(3)))/Tmass**2 
     &- (2.0d0*fo(4)*(km(1) - km(4)))/Tmass**2) + 
     &    fi(4)*((-2.0d0*fo(1)*(-km(2) + cone*km(3)))/Tmass**2
     & - (2.0d0*fo(2)*(km(1) + km(4)))/Tmass**2) + 
     &    fi(1)*((-2.0d0*fo(4)*(km(2) + cone*km(3)))/Tmass**2 
     &- (2.0d0*fo(3)*(km(1) + km(4)))/Tmass**2) 
      
      do i = 1,4
         do j=1,4
            yio(i,j) = tt1*MET(i,j)+tt2*km(i)*km(j)
     &+tt3*(km(i)*kp(j)+km(j)*kp(i))
         end do
      enddo
	
      yio(1,1) = yio(1,1)
     &+fi(3)*fo(1)*((-4.0d0*(K1KT + K2KT)*km(1))/Tmass**2 + 4.0d0*kp(1))  
     & +fi(4)*fo(2)*((-4.d0*(K1KT + K2KT)*km(1))/Tmass**2 + 4.0d0*kp(1)) 
     & +fi(1)*fo(3)*((-4.0d0*(K1KT + K2KT)*km(1))/Tmass**2+ 4.0d0*kp(1)) 
     &+ fi(2)*fo(4)*((-4.0d0*(K1KT + K2KT)*km(1))/Tmass**2+ 4.0d0*kp(1))
	
      yio(2,2) = yio(2,2)
     &+fi(2)*fo(3)*((4.0d0*(K1KT + K2KT)*km(2))/Tmass**2 - 4.d0*kp(2))
     &+fi(1)*fo(4)*((4.0d0*(K1KT + K2KT)*km(2))/Tmass**2 - 4.d0*kp(2))  
     &+fi(4)*fo(1)*((-4.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 4.d0*kp(2))
     &+fi(3)*fo(2)*((-4.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 4.d0*kp(2))

      yio(3,3) = yio(3,3)
     &+fi(4)*fo(1)*((4.0d0*cone*(K1KT + K2KT)*km(3))/Tmass**2
     & - 4.0d0*cone*kp(3)) + 
     &    fi(1)*fo(4)*((4.0d0*cone*(K1KT + K2KT)*km(3))/Tmass**2
     & - 4.0d0*cone*kp(3)) + 
     &    fi(3)*fo(2)*((-4.0d0*cone*(K1KT + K2KT)*km(3))/Tmass**2
     & + 4.d0*cone*kp(3)) + 
     &    fi(2)*fo(3)*((-4.d0*cone*(K1KT + K2KT)*km(3))/Tmass**2
     & + 4.d0*cone*kp(3))	
	
      yio(4,4) = yio(4,4)
     &+fi(4)*fo(2)*((4.d0*(K1KT + K2KT)*km(4))/Tmass**2 - 4.d0*kp(4)) 
     &+fi(1)*fo(3)*((4.d0*(K1KT + K2KT)*km(4))/Tmass**2 - 4.d0*kp(4)) 
     &+fi(3)*fo(1)*((-4.d0*(K1KT + K2KT)*km(4))/Tmass**2+ 4.d0*kp(4))
     &+fi(2)*fo(4)*((-4.d0*(K1KT + K2KT)*km(4))/Tmass**2+ 4.d0*kp(4))	

      yio(1,2) = yio(1,2)
     &+fi(3)*(fo(2)*((-2.d0*(K1KT + K2KT)*km(1))/Tmass**2 + 2.d0*kp(1)) 
     &+ fo(1)*((-2.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 2.d0*kp(2))) + 
     &  fi(4)*(fo(1)*((-2.d0*(K1KT + K2KT)*km(1))/Tmass**2 + 2.d0*kp(1))
     & + fo(2)*((-2.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 2.d0*kp(2))) + 
     &  fi(1)*(fo(4)*((2.d0*(K1KT + K2KT)*km(1))/Tmass**2 - 2.d0*kp(1))
     & + fo(3)*((-2.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 2.d0*kp(2))) + 
     & fi(2)*(fo(3)*((2.d0*(K1KT + K2KT)*km(1))/Tmass**2 - 2.d0*kp(1))
     & + fo(4)*((-2.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 2.d0*kp(2)))

      yio(1,3) = yio(1,3)
     &+fi(3)*(fo(2)*((-2.0d0*cone*(K1KT + K2KT)*km(1))/Tmass**2
     & + 2.d0*cone*kp(1)) + 
     &    fo(1)*((-2.d0*(K1KT + K2KT)*km(3))/Tmass**2 +2.d0*kp(3))) + 
     &    fi(4)*(fo(1)*((2.d0*cone*(K1KT + K2KT)*km(1))/Tmass**2
     & - 2.d0*cone*kp(1)) + 
     &    fo(2)*((-2.d0*(K1KT + K2KT)*km(3))/Tmass**2 +2.d0*kp(3))) + 
     &    fi(1)*(fo(4)*((2.d0*cone*(K1KT + K2KT)*km(1))/Tmass**2 
     &- 2.d0*cone*kp(1)) + 
     &    fo(3)*((-2.d0*(K1KT + K2KT)*km(3))/Tmass**2 +2.d0*kp(3))) + 
     &    fi(2)*(fo(3)*((-2.d0*cone*(K1KT + K2KT)*km(1))/Tmass**2 
     &+ 2.d0*cone*kp(1)) + 
     &    fo(4)*((-2.d0*(K1KT + K2KT)*km(3))/Tmass**2 + 2.d0*kp(3)))
  
      yio(1,4) = yio(1,4)
     &+fi(4)*fo(2)*((-2.d0*(K1KT + K2KT)*(-km(1) + km(4)))/Tmass**2
     & - 2.d0*kp(1) + 2.d0*kp(4)) + 
     & fi(1)*fo(3)*((-2.d0*(K1KT + K2KT)*(-km(1) + km(4)))/Tmass**2
     & - 2.d0*kp(1) + 2.d0*kp(4)) + 
     & fi(3)*fo(1)*((-2.d0*(K1KT + K2KT)*(km(1) + km(4)))/Tmass**2 
     &+ 2.d0*kp(1) + 2.d0*kp(4)) + 
     & fi(2)*fo(4)*((-2.d0*(K1KT + K2KT)*(km(1) + km(4)))/Tmass**2
     & + 2.d0*kp(1) + 2.d0*kp(4))

      yio(2,3) = yio(2,3)
     &+fi(1)*fo(4)*((-2.0d0*(K1KT + K2KT)*(-(cone*km(2))
     & - km(3)))/Tmass**2 - 2.0d0*cone*kp(2) - 2.d0*kp(3)) + 
     &    fi(2)*fo(3)*((-2.d0*(K1KT + K2KT)*(cone*km(2)
     & - km(3)))/Tmass**2 + 2.d0*cone*kp(2) - 2.d0*kp(3)) + 
     &    fi(4)*fo(1)*((-2.d0*(K1KT + K2KT)*(-(cone*km(2))
     & + km(3)))/Tmass**2 - 2.d0*cone*kp(2) + 2.d0*kp(3)) + 
     &    fi(3)*fo(2)*((-2.d0*(K1KT + K2KT)*(cone*km(2)
     & + km(3)))/Tmass**2 + 2.d0*cone*kp(2) + 2.d0*kp(3))


      yio(2,4) = yio(2,4)
     &+fi(2)*(fo(4)*((-2.d0*(K1KT + K2KT)*km(2))/Tmass**2+2.d0*kp(2)) + 
     &       fo(3)*((2.d0*(K1KT + K2KT)*km(4))/Tmass**2 - 2.d0*kp(4))) + 
     &  fi(1)*(fo(3)*((2.d0*(K1KT + K2KT)*km(2))/Tmass**2 - 2.d0*kp(2))
     & + fo(4)*((2.d0*(K1KT + K2KT)*km(4))/Tmass**2 - 2.d0*kp(4))) + 
     &  fi(4)*(fo(2)*((2.d0*(K1KT + K2KT)*km(2))/Tmass**2 - 2.d0*kp(2))
     & + fo(1)*((-2.d0*(K1KT + K2KT)*km(4))/Tmass**2 + 2.d0*kp(4))) + 
     & fi(3)*(fo(1)*((-2.d0*(K1KT + K2KT)*km(2))/Tmass**2 + 2.d0*kp(2))
     & + fo(2)*((-2.d0*(K1KT + K2KT)*km(4))/Tmass**2 + 2.d0*kp(4)))

      yio(3,4) = yio(3,4)
     &+fi(4)*(fo(2)*((2.d0*(K1KT + K2KT)*km(3))/Tmass**2- 2.d0*kp(3)) + 
     &       fo(1)*((2.d0*cone*(K1KT + K2KT)*km(4))/Tmass**2
     & - 2.d0*cone*kp(4))) + 
     &    fi(1)*(fo(3)*((2.d0*(K1KT + K2KT)*km(3))/Tmass**2
     & - 2.d0*kp(3)) + 
     &       fo(4)*((2.d0*cone*(K1KT + K2KT)*km(4))/Tmass**2
     & - 2.d0*cone*kp(4))) + 
     &    fi(3)*(fo(1)*((-2.d0*(K1KT + K2KT)*km(3))/Tmass**2
     & + 2.d0*kp(3)) + 
     &       fo(2)*((-2.d0*cone*(K1KT + K2KT)*km(4))/Tmass**2
     & + 2.d0*cone*kp(4))) + 
     &    fi(2)*(fo(4)*((-2.d0*(K1KT + K2KT)*km(3))/Tmass**2
     & + 2.d0*kp(3)) + 
     &       fo(3)*((-2.d0*cone*(K1KT + K2KT)*km(4))/Tmass**2
     &+ 2.d0*cone*kp(4)))

      yio(2,1) = yio(1,2)
      yio(3,1) = yio(1,3)
      yio(4,1) = yio(1,4)
      
      yio(3,2) = yio(2,3)
      yio(4,2) = yio(2,4)
      yio(4,3) = yio(3,4)
     
      do i = 1,4
         do j=1,4
            yio(i,j) = yio(i,j)*d/2.0d0
         end do
      enddo

      ft(1) = yio(1,1)
      ft(2) = yio(1,2)
      ft(3) = yio(1,3)
      ft(4) = yio(1,4)
      ft(5) = yio(2,1)
      ft(6) = yio(2,2)
      ft(7) = yio(2,3)
      ft(8) = yio(2,4)
      ft(9) = yio(3,1)
      ft(10) = yio(3,2)
      ft(11) = yio(3,3)
      ft(12) = yio(3,4)
      ft(13) = yio(4,1)
      ft(14) = yio(4,2)
      ft(15) = yio(4,3)
      ft(16) = yio(4,4)
      ft(17) = yio(5,1)
      ft(18) = yio(6,1)

      return
      end

      subroutine uvvxxx(vc1,vc2,g,vmass,Tmass,Twidth , ft)
c
c This subroutine computes an off-shell Tensor boson wavfunction from 
c the two gauge bosons-Kaluza-Klein tensor boson coupling.
c
c input:
c       complex vc1(3)          : first input vector boson          v1
c       complex vc2(3)          : second input vector boson         v2
c       real    g              : coupling constant    -kappa/2
c       real    vmass          : mass  of input vector bosons
c       real    Tmass          : mass  of OUTPUT KK Tensor T'
c       real    Twidth         : width of OUTPUT KK Tensor boson T'
c output:
c       complex ft(18)       : KK tensor current       j(T:v1,v2)     
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex vc1(6),vc2(6),yvv(6,4),ft(18)
      double complex KTE1,KTE2,K2E1,K1E2,K1E1,K2E2,E1E2
      
      double precision g, vmass,Tmass,Twidth
      integer i,j
      double precision pv1(4), pv2(4), pT(4)
      double precision MET(4,4)
      double complex cZero, d
      double precision rZero, rTwo
      double precision K1K2,KT2,K1KT,K2KT
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      yvv(5,1) = vc1(5)+vc2(5)
      yvv(6,1) = vc1(6)+vc2(6)
      
      pv1(1) = dreal(vc1(5))
      pv1(2) = dreal(vc1(6))
      pv1(3) = dimag(vc1(6))
      pv1(4) = dimag(vc1(5))
      
      pv2(1) = dreal(vc2(5))
      pv2(2) = dreal(vc2(6))
      pv2(3) = dimag(vc2(6))
      pv2(4) = dimag(vc2(5))
      
      pT(1) = dreal(yvv(5,1))
      pT(2) = dreal(yvv(6,1))
      pT(3) = dimag(yvv(6,1))
      pT(4) = dimag(yvv(5,1))
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      K1K2 = pv1(1)*pv2(1)-pv1(2)*pv2(2)-pv1(3)*pv2(3)-pv1(4)*pv2(4)
      K1KT = pv1(1)*pT(1)-pv1(2)*pT(2)-pv1(3)*pT(3)-pv1(4)*pT(4)
      K2KT = pT(1)*pv2(1)-pT(2)*pv2(2)-pT(3)*pv2(3)-pT(4)*pv2(4)
      KT2 = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
      
      KTE1 = pT(1)*vc1(1)-pT(2)*vc1(2)-pT(3)*vc1(3)-pT(4)*vc1(4)
      KTE2 = pT(1)*vc2(1)-pT(2)*vc2(2)-pT(3)*vc2(3)-pT(4)*vc2(4)
      K1E1 = pv1(1)*vc1(1)-pv1(2)*vc1(2)-pv1(3)*vc1(3)-pv1(4)*vc1(4)
      K1E2 = pv1(1)*vc2(1)-pv1(2)*vc2(2)-pv1(3)*vc2(3)-pv1(4)*vc2(4)
      K2E1 = pv2(1)*vc1(1)-pv2(2)*vc1(2)-pv2(3)*vc1(3)-pv2(4)*vc1(4)
      K2E2 = pv2(1)*vc2(1)-pv2(2)*vc2(2)-pv2(3)*vc2(3)-pv2(4)*vc2(4)
      E1E2 = vc2(1)*vc1(1)-vc2(2)*vc1(2)-vc2(3)*vc1(3)-vc2(4)*vc1(4)
      
      if ( Tmass.gt.rZero ) then
         d = - g/dcmplx( KT2-Tmass**2, Tmass*Twidth )
      else
         d = - g/dcmplx( KT2, rZero )
      end if

    
      do i = 1,4
         do j=1,4
            yvv(i,j) = 2.0d0*K1K2*(vc1(i)*vc2(j)+vc1(j)*vc2(i))
     &-2.0d0*K1K2*KTE2/Tmass**2*(PT(i)*vc1(j)+PT(j)*vc1(i))
     &-2.0d0*K1E2*(pv2(i)*vc1(j)+pv2(j)*vc1(i))
     &+2.0d0*K1E2*K2KT/Tmass**2*(PT(i)*vc1(j)+PT(j)*vc1(i))
     &-2.0d0/3.0d0*E1E2*K1K2*MET(i,j)
     &+8.0d0/3.0d0*K1K2*E1E2/Tmass**2*PT(i)*PT(j)
     &+2.0d0*E1E2*(pv1(i)*pv2(j)+pv1(j)*pv2(i)) 
     &-2.0d0*K1K2*KTE1/Tmass**2*(PT(i)*vc2(j)+PT(j)*vc2(i))
     &-2.0d0*K2E1*(pv1(i)*vc2(j)+pv1(j)*vc2(i))
     &+4.0d0*K1K2*KTE1*KTE2/3.0d0/Tmass**2*MET(i,j)
     &+8.0d0*K1K2*KTE1*KTE2/3.0d0/Tmass**4*PT(i)*PT(j)
     &+2.0d0*K2E1*KTE2/Tmass**2*(pv1(i)*PT(j)+pv1(j)*PT(i))
     &+2.0d0*KTE1*K1E2/Tmass**2*(pv2(i)*PT(j)+pv2(j)*PT(i))
     &+2.0d0*K2E1*K1E2*MET(i,j)
     &-4.0d0*K2E1*K1E2/Tmass**2*PT(i)*PT(j)
     &-2.0d0/3.0d0/Tmass**2*KT2*E1E2*K1K2*MET(i,j)
     &-4.0d0/3.0d0/Tmass**4*KT2*E1E2*K1K2*PT(i)*PT(j)
     &+2.0d0/3.0d0/Tmass**2*K2E1*K1E2*KT2*MET(i,j)
     &+4.0d0/3.0d0/Tmass**4*K2E1*K1E2*KT2*PT(i)*PT(j)
     &-2.0d0/Tmass**2*E1E2*K1KT*(pv2(i)*PT(j)+pv2(j)*PT(i))
     &+2.0d0/Tmass**2*K2E1*K1KT*(PT(i)*vc2(j)+PT(j)*vc2(i))
     &-4.0d0/3.0d0/Tmass**2*K2E1*KTE2*K1KT*MET(i,j)
     &-8.0d0/3.0d0/Tmass**4*K2E1*KTE2*K1KT*PT(i)*PT(j)
     &-2.0d0/Tmass**2*E1E2*K2KT*(PT(i)*pv1(j)+PT(j)*pv1(i))
     &-4.0d0/3.0d0/Tmass**2*KTE1*K1E2*K2KT*MET(i,j)
     &-8.0d0/3.0d0/Tmass**4*KTE1*K1E2*K2KT*PT(i)*PT(j)
     &+4.0d0/3.0d0/Tmass**2*E1E2*K2KT*K1KT*MET(i,j)
     &+8.0d0/3.0d0/Tmass**4*E1E2*K2KT*K1KT*PT(i)*PT(j)
     &-4.0d0/3.0d0*E1E2*K1K2*MET(i,j)
     &+4.0d0/3.0d0/Tmass**2*E1E2*K1K2
     &*PT(i)*PT(j)
            
            if ( vmass.ne.rZero ) then
               yvv(i,j) = 
     &              yvv(i,j)+2*vmass**2*(vc1(i)*vc2(j)+vc1(j)*vc2(i))
     &-2.0d0/3.0d0*vmass**2*E1E2*MET(i,j)
     &+8.0d0/3.0d0/Tmass**2*vmass**2*E1E2*PT(i)*PT(j)
     &-2.0d0/Tmass**2*vmass**2*KTE1*(PT(i)*vc2(j)+PT(j)*vc2(i)) 
     &-2.0d0/Tmass**2*vmass**2*KTE2*(PT(i)*vc1(j)+PT(j)*vc1(i)) 
     &+4.0d0/3.0d0/Tmass**2*vmass**2*KTE1*KTE2*MET(i,j)
     &+8.0d0/3.0d0/Tmass**4*vmass**2*KTE1*KTE2*PT(i)*PT(j)
     &-2.0d0/3.0d0/Tmass**2*vmass**2*KT2*E1E2*MET(i,j)
     &-4.0d0/3.0d0/Tmass**4*vmass**2*KT2*E1E2*PT(i)*PT(j)
            else
c     gauge fixing term for zero mass photon/gluon
               yvv(i,j) = 
     &              yvv(i,j)-2.0d0*K1E1*(pv1(i)*vc2(j)+pv1(j)*vc2(i))
     &+2.0d0/Tmass**2*KTE2*K1E1*(PT(i)*pv1(j)+PT(j)*pv1(i))
     &+2.0d0/3.0d0*K1E2*K1E1*MET(i,j)
     &-8.0d0/3.0d0/Tmass**2*K1E2*K1E1*PT(i)*PT(j)
     &-2.0d0/3.0d0*K2E2*K1E1*MET(i,j)
     &-4.0d0/3.0d0/Tmass**2*K2E2*K1E1*PT(i)*PT(j)
     &+2.0d0/3.0d0/Tmass**2*K1E2*KT2*K1E1*MET(i,j)
     &+4.0d0/3.0d0/Tmass**4*K1E2*KT2*K1E1*PT(i)*PT(j)
     &+2.0d0/3.0d0/Tmass**2*K2E2*KT2*K1E1*MET(i,j)
     &+4.0d0/3.0d0/Tmass**4*K2E2*KT2*K1E1*PT(i)*PT(j)
     &+2.0d0/Tmass**2*K1KT*K1E1*(PT(i)*vc2(j)+PT(j)*vc2(i))
     &-4.0d0/3.0d0/Tmass**2*K1KT*KTE2*K1E1*MET(i,j)
     &-8.0d0/3.0d0/Tmass**4*K1KT*KTE2*K1E1*PT(i)*PT(j)
     &-2.0d0*K2E2*(pv2(i)*vc1(j)+pv2(j)*vc1(i))
     &+2.0d0/Tmass**2*K2E2*KTE1*(pv2(i)*PT(j)+pv2(j)*PT(i))
     &+2.0d0/3.0d0*K2E1*K2E2*MET(i,j)
     &-8.0d0/3.0d0/Tmass**2*K2E1*K2E2*PT(i)*PT(j)
     &+2.0d0/3.0d0/Tmass**2*KT2*K2E2*K2E1*MET(i,j)
     &+4.0d0/3.0d0/Tmass**4*KT2*K2E2*K2E1*PT(i)*PT(j)
     &+2.0d0/Tmass**2*K2E2*K2KT*(PT(i)*vc1(j)+PT(j)*vc1(i))
     &-4.0d0/3.0d0/Tmass**2*K2E2*K2KT*KTE1*MET(i,j)
     &-8.0d0/3.0d0/Tmass**4*K2E2*K2KT*KTE1
     &*PT(i)*PT(j)

            endif

            yvv(i,j) = yvv(i,j)*d/2.0d0

         end do
      enddo
      
      ft(1) = yvv(1,1)
      ft(2) = yvv(1,2)
      ft(3) = yvv(1,3)
      ft(4) = yvv(1,4)
      ft(5) = yvv(2,1)
      ft(6) = yvv(2,2)
      ft(7) = yvv(2,3)
      ft(8) = yvv(2,4)
      ft(9) = yvv(3,1)
      ft(10) = yvv(3,2)
      ft(11) = yvv(3,3)
      ft(12) = yvv(3,4)
      ft(13) = yvv(4,1)
      ft(14) = yvv(4,2)
      ft(15) = yvv(4,3)
      ft(16) = yvv(4,4)
      ft(17) = yvv(5,1)
      ft(18) = yvv(6,1)

      return
      end

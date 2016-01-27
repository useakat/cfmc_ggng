      subroutine jvtxxx(vc,ft,g,vmass,vwidth , jtv)
c
c This subroutine computes an off-shell vector wavfunction from the coupling of
c two gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex vc(6)          : input vector               V
c       complex ft(18)         : tensor               KK mode T
c       real    g              : coupling constant    -kappa/2
c       real    vmass          : mass  of output vector v'
c       real    vwidth         : width of output vector v'
c output:
c       complex jtv(6)         : vector current         j^mu(v':v,tkk)
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex vc(6),ft(18),tc(6,4),jtv(6),TVM(4),TKM(4)
      double precision g,vmass,vwidth
      double precision MET(4,4)
      double complex T12, T13, T14, T23, T24, T34, T00
      double complex K2V1,K1V1
      double complex TKK,TK2V1, dum
      double precision pv1(4), pv2(4), F, pp2

      integer i, j

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
 
      jtv(5) = vc(5)+tc(5,1)
      jtv(6) = vc(6)+tc(6,1)

      pv1(1) = -dreal(vc(5))
      pv1(2) = -dreal(vc(6))
      pv1(3) = -dimag(vc(6))
      pv1(4) = -dimag(vc(5))

      pv2(1) = dreal(jtv(5))
      pv2(2) = dreal(jtv(6))
      pv2(3) = dimag(jtv(6))
      pv2(4) = dimag(jtv(5))

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      K2V1 = pv2(1)*vc(1) - pv2(2)*vc(2) - pv2(3)*vc(3) - pv2(4)*vc(4)
      K1V1 = pv1(1)*vc(1) - pv1(2)*vc(2) - pv1(3)*vc(3) - pv1(4)*vc(4)
      F = pv1(1)*pv2(1) - pv1(2)*pv2(2) - pv1(3)*pv2(3) - pv1(4)*pv2(4)
      pp2 = pv2(1)**2 - pv2(2)**2 - pv2(3)**2 - pv2(4)**2

      TKK   = cZero
      TK2V1 = cZero

      do i = 1,4
         dum   = tc(i,i)*pv2(i)
         TKK   = TKK   + dum*pv1(i)
         dum   = tc(i,i)*vc(i)
         TK2V1 = TK2V1 + dum*pv2(i)
      end do
      
      TKK   = rTwo*TKK
      TK2V1 = rTwo*TK2V1
      
      TKK = TKK - T12*(pv1(1)*pv2(2) + pv1(2)*pv2(1))
     &          - T13*(pv1(1)*pv2(3) + pv1(3)*pv2(1))
     &          - T14*(pv1(1)*pv2(4) + pv1(4)*pv2(1))
     &          + T23*(pv1(2)*pv2(3) + pv1(3)*pv2(2))
     &          + T24*(pv1(2)*pv2(4) + pv1(4)*pv2(2))
     &          + T34*(pv1(3)*pv2(4) + pv1(4)*pv2(3))

      TK2V1 = TK2V1 - T12*(vc(1)*pv2(2) + vc(2)*pv2(1))
     &              - T13*(vc(1)*pv2(3) + vc(3)*pv2(1))
     &              - T14*(vc(1)*pv2(4) + vc(4)*pv2(1))
     &              + T23*(vc(2)*pv2(3) + vc(3)*pv2(2))
     &              + T24*(vc(2)*pv2(4) + vc(4)*pv2(2))
     &              + T34*(vc(3)*pv2(4) + vc(4)*pv2(3))


      do j=1,4

         TVM(j) =
     &MET(j,1)*(tc(1,1)*vc(1)-tc(2,1)*vc(2)
     &-tc(3,1)*vc(3)-tc(4,1)*vc(4))
     &-MET(j,2)*(tc(1,2)*vc(1)-tc(2,2)*vc(2)
     &-tc(3,2)*vc(3)-tc(4,2)*vc(4))
     &-MET(j,3)*(tc(1,3)*vc(1)-tc(2,3)*vc(2)
     &-tc(3,3)*vc(3)-tc(4,3)*vc(4))
     &-MET(j,4)*(tc(1,4)*vc(1)-tc(2,4)*vc(2)
     &-tc(3,4)*vc(3)-tc(4,4)*vc(4))
     &+
     &MET(j,1)*(tc(1,1)*vc(1)-tc(1,2)*vc(2)
     &-tc(1,3)*vc(3)-tc(1,4)*vc(4))
     &-MET(j,2)*(tc(2,1)*vc(1)-tc(2,2)*vc(2)
     &-tc(2,3)*vc(3)-tc(2,4)*vc(4))
     &-MET(j,3)*(tc(3,1)*vc(1)-tc(3,2)*vc(2)
     &-tc(3,3)*vc(3)-tc(3,4)*vc(4))
     &-MET(j,4)*(tc(4,1)*vc(1)-tc(4,2)*vc(2)
     &-tc(4,3)*vc(3)-tc(4,4)*vc(4))

         TKM(j) =
     &MET(j,1)*(tc(1,1)*pv1(1)-tc(2,1)*pv1(2)
     &-tc(3,1)*pv1(3)-tc(4,1)*pv1(4))
     &-MET(j,2)*(tc(1,2)*pv1(1)-tc(2,2)*pv1(2)
     &-tc(3,2)*pv1(3)-tc(4,2)*pv1(4))
     &-MET(j,3)*(tc(1,3)*pv1(1)-tc(2,3)*pv1(2)
     &-tc(3,3)*pv1(3)-tc(4,3)*pv1(4))
     &-MET(j,4)*(tc(1,4)*pv1(1)-tc(2,4)*pv1(2)
     &-tc(3,4)*pv1(3)-tc(4,4)*pv1(4))
     &+
     &MET(j,1)*(tc(1,1)*pv1(1)-tc(1,2)*pv1(2)
     &-tc(1,3)*pv1(3)-tc(1,4)*pv1(4))
     &-MET(j,2)*(tc(2,1)*pv1(1)-tc(2,2)*pv1(2)
     &-tc(2,3)*pv1(3)-tc(2,4)*pv1(4))
     &-MET(j,3)*(tc(3,1)*pv1(1)-tc(3,2)*pv1(2)
     &-tc(3,3)*pv1(3)-tc(3,4)*pv1(4))
     &-MET(j,4)*(tc(4,1)*pv1(1)-tc(4,2)*pv1(2)
     &-tc(4,3)*pv1(3)-tc(4,4)*pv1(4))
     
      enddo

      if ( vmass.ne.rZero ) then

         do i=1,4

            jtv(i) = -(vmass**2+F)*T00*vc(i)
     &+T00*K2V1*(pv1(i)+(1.0d0+F/vmass**2)*pv2(i)
     &-F/vmass**2*pv2(i))
     &+(vmass**2+F)*TVM(i)
     &-TK2V1*pv1(i)
     &-TK2V1*pv2(i)*(1.0d0+F/vmass**2)
     &+TKK*vc(i)
     &-K2V1*TKM(i)
     &+F/vmass**2*TK2V1*pv2(i)
	 
            jtv(i)=jtv(i)*g/dcmplx(pp2-vmass**2, vmass*vwidth )

         enddo

      else
	
         do i=1,4

            jtv(i) = -F*T00*vc(i)
     &+K1V1*T00*(pv1(i)+pv2(i))
     &+T00*K2V1*(pv1(i)+pv2(i))
     &+F*TVM(i)
     &-TK2V1*(pv1(i)+pv2(i))
     &+TKK*vc(i)
     &-(K2V1+K1V1)*TKM(i)
  
            jtv(i)=jtv(i)*g/dcmplx(pp2, 0.0d0)
            
         enddo

         
      endif

      return
      end

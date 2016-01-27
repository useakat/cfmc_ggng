      subroutine jvvtxx(v1,v2,ft,g,g2,vmass,vwidth , jvvt)
c
c This subroutine computes an off-shell vector boson wavefunction from
c the four-point coupling of
c three gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex v1(6)          : first  vector           flow-out V
c       complex v2(6)          : second vector           flow-out V
c       complex ft(18)         : tensor                  KK mode T
c       complex g2             : coupling constant       -kappa/2
c       real g                 : coupling constant       (gs)
c       real    vmass          : mass  of output vector boson
c       real    vwidth         : width of output vector boson

c
c output:
c       complex jvvt(6)         : off-shell vector boson            
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex v1(6), v2(6),ft(18),tc(6,4), jvvt(6),
     &	TV1M(4),TV2M(4),TK12M(4)
	double complex g2 
      double precision g,vmass,vwidth
	double precision MET(4,4)
      double complex d,T00, T12, T13, T14, T23, T24, T34
      double complex V1V2,K1V2, K2V1
     &,K3V1,K3V2
      double precision K1K2

      double complex TV12,TKV1,TKV2,
     &TK3V1,TK3V2,TK312,dum
      double precision pv1(4), pv2(4), pv3(4),
     &pv32,p31(4),p23(4),p12(4),K1K3,K2K3

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


      jvvt(5) = v1(5)+v2(5)+tc(5,1)
      jvvt(6) = v1(6)+v2(6)+tc(6,1)

      pv1(1) = dreal(v1(5))
      pv1(2) = dreal(v1(6))
      pv1(3) = dimag(v1(6))
      pv1(4) = dimag(v1(5))

      pv2(1) = dreal(v2(5))
      pv2(2) = dreal(v2(6))
      pv2(3) = dimag(v2(6))
      pv2(4) = dimag(v2(5))

      pv3(1) = -dreal(jvvt(5))
      pv3(2) = -dreal(jvvt(6))
      pv3(3) = -dimag(jvvt(6))
      pv3(4) = -dimag(jvvt(5))
	
      pv32=pv3(1)**2-pv3(2)**2-pv3(3)**2-pv3(4)**2
      if ( vmass.gt.rZero ) then
         d =  g/dcmplx( pv32-vmass**2, vmass*vwidth )
      else
         d =  g/dcmplx( pv32, rZero )
      end if

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) = 1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      p31(1) = pv3(1)-pv1(1)
      p31(2) = pv3(2)-pv1(2)
      p31(3) = pv3(3)-pv1(3)
      p31(4) = pv3(4)-pv1(4)
      
      p12(1) = pv1(1)-pv2(1)
      p12(2) = pv1(2)-pv2(2)
      p12(3) = pv1(3)-pv2(3)
      p12(4) = pv1(4)-pv2(4)
      
      p23(1) = pv2(1)-pv3(1)
      p23(2) = pv2(2)-pv3(2)
      p23(3) = pv2(3)-pv3(3)
      p23(4) = pv2(4)-pv3(4)
      
      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)
      
      V1V2 =  v1(1)*v2(1) -  v1(2)*v2(2) -  v1(3)*v2(3) -  v1(4)*v2(4)
      K1V2 = pv1(1)*v2(1) - pv1(2)*v2(2) - pv1(3)*v2(3) - pv1(4)*v2(4)
      K2V1 = pv2(1)*v1(1) - pv2(2)*v1(2) - pv2(3)*v1(3) - pv2(4)*v1(4)
      K3V1 = pv3(1)*v1(1) - pv3(2)*v1(2) - pv3(3)*v1(3) - pv3(4)*v1(4)
      K3V2 = pv3(1)*v2(1) - pv3(2)*v2(2) - pv3(3)*v2(3) - pv3(4)*v2(4)
      
      K1K3 = pv1(1)*pv3(1)-pv1(2)*pv3(2)-pv1(3)*pv3(3)-pv1(4)*pv3(4)
      K2K3 = pv2(1)*pv3(1)-pv2(2)*pv3(2)-pv2(3)*pv3(3)-pv2(4)*pv3(4)
      
      TV12 = rtwo*(tc(1,1)*v1(1)*v2(1)+tc(2,2)*v1(2)*v2(2)
     &+tc(3,3)*v1(3)*v2(3)+tc(4,4)*v1(4)*v2(4))


      TKV1 = rtwo*(tc(1,1)*p23(1)*v1(1)+tc(2,2)*p23(2)*v1(2)
     &+tc(3,3)*p23(3)*v1(3)+tc(4,4)*p23(4)*v1(4))

      TKV2 = rtwo*(tc(1,1)*p31(1)*v2(1)+tc(2,2)*p31(2)*v2(2)
     &+tc(3,3)*p31(3)*v2(3)+tc(4,4)*p31(4)*v2(4))

     
      TK3V1 = rtwo*(tc(1,1)*pv3(1)*v1(1)+tc(2,2)*pv3(2)*v1(2)
     &+tc(3,3)*pv3(3)*v1(3)+tc(4,4)*pv3(4)*v1(4))	

      TK3V2 = rtwo*(tc(1,1)*pv3(1)*v2(1)+tc(2,2)*pv3(2)*v2(2)
     &+tc(3,3)*pv3(3)*v2(3)+tc(4,4)*pv3(4)*v2(4))	


      TK312 = rtwo*(tc(1,1)*pv3(1)*p12(1)+tc(2,2)*pv3(2)*p12(2)
     &+tc(3,3)*pv3(3)*p12(3)+tc(4,4)*pv3(4)*p12(4))	


      TV12 = TV12 - T12*(v1(1)*v2(2) + v1(2)*v2(1))
     &          - T13*(v1(1)*v2(3) + v1(3)*v2(1))
     &          - T14*(v1(1)*v2(4) + v1(4)*v2(1))
     &          + T23*(v1(2)*v2(3) + v1(3)*v2(2))
     &          + T24*(v1(2)*v2(4) + v1(4)*v2(2))
     &          + T34*(v1(3)*v2(4) + v1(4)*v2(3))


      TKV1 = TKV1 - T12*(p23(1)*v1(2) + p23(2)*v1(1))
     &              - T13*(p23(1)*v1(3) + p23(3)*v1(1))
     &              - T14*(p23(1)*v1(4) + p23(4)*v1(1))
     &              + T23*(p23(2)*v1(3) + p23(3)*v1(2))
     &              + T24*(p23(2)*v1(4) + p23(4)*v1(2))
     &              + T34*(p23(3)*v1(4) + p23(4)*v1(3))

      TKV2 = TKV2 - T12*(p31(1)*v2(2) + p31(2)*v2(1))
     &              - T13*(p31(1)*v2(3) + p31(3)*v2(1))
     &              - T14*(p31(1)*v2(4) + p31(4)*v2(1))
     &              + T23*(p31(2)*v2(3) + p31(3)*v2(2))
     &              + T24*(p31(2)*v2(4) + p31(4)*v2(2))
     &              + T34*(p31(3)*v2(4) + p31(4)*v2(3))

      TK3V1 = TK3V1 - T12*(pv3(1)*v1(2) + pv3(2)*v1(1))
     &              - T13*(pv3(1)*v1(3) + pv3(3)*v1(1))
     &              - T14*(pv3(1)*v1(4) + pv3(4)*v1(1))
     &              + T23*(pv3(2)*v1(3) + pv3(3)*v1(2))
     &              + T24*(pv3(2)*v1(4) + pv3(4)*v1(2))
     &              + T34*(pv3(3)*v1(4) + pv3(4)*v1(3))

      TK3V2 = TK3V2 - T12*(pv3(1)*v2(2) + pv3(2)*v2(1))
     &              - T13*(pv3(1)*v2(3) + pv3(3)*v2(1))
     &              - T14*(pv3(1)*v2(4) + pv3(4)*v2(1))
     &              + T23*(pv3(2)*v2(3) + pv3(3)*v2(2))
     &              + T24*(pv3(2)*v2(4) + pv3(4)*v2(2))
     &              + T34*(pv3(3)*v2(4) + pv3(4)*v2(3))

      TK312 = TK312 - T12*(pv3(1)*p12(2) + pv3(2)*p12(1))
     &              - T13*(pv3(1)*p12(3) + pv3(3)*p12(1))
     &              - T14*(pv3(1)*p12(4) + pv3(4)*p12(1))
     &              + T23*(pv3(2)*p12(3) + pv3(3)*p12(2))
     &              + T24*(pv3(2)*p12(4) + pv3(4)*p12(2))
     &              + T34*(pv3(3)*p12(4) + pv3(4)*p12(3))

      do j=1,4

         TV1M(j) =
     &MET(j,1)*(tc(1,1)*v1(1)-tc(2,1)*v1(2)
     &-tc(3,1)*v1(3)-tc(4,1)*v1(4))
     &-MET(j,2)*(tc(1,2)*v1(1)-tc(2,2)*v1(2)
     &-tc(3,2)*v1(3)-tc(4,2)*v1(4))
     &-MET(j,3)*(tc(1,3)*v1(1)-tc(2,3)*v1(2)
     &-tc(3,3)*v1(3)-tc(4,3)*v1(4))
     &-MET(j,4)*(tc(1,4)*v1(1)-tc(2,4)*v1(2)
     &-tc(3,4)*v1(3)-tc(4,4)*v1(4))
     &+
     &MET(j,1)*(tc(1,1)*v1(1)-tc(1,2)*v1(2)
     &-tc(1,3)*v1(3)-tc(1,4)*v1(4))
     &-MET(j,2)*(tc(2,1)*v1(1)-tc(2,2)*v1(2)
     &-tc(2,3)*v1(3)-tc(2,4)*v1(4))
     &-MET(j,3)*(tc(3,1)*v1(1)-tc(3,2)*v1(2)
     &-tc(3,3)*v1(3)-tc(3,4)*v1(4))
     &-MET(j,4)*(tc(4,1)*v1(1)-tc(4,2)*v1(2)
     &-tc(4,3)*v1(3)-tc(4,4)*v1(4))

         TV2M(j) =
     &MET(j,1)*(tc(1,1)*v2(1)-tc(2,1)*v2(2)
     &-tc(3,1)*v2(3)-tc(4,1)*v2(4))
     &-MET(j,2)*(tc(1,2)*v2(1)-tc(2,2)*v2(2)
     &-tc(3,2)*v2(3)-tc(4,2)*v2(4))
     &-MET(j,3)*(tc(1,3)*v2(1)-tc(2,3)*v2(2)
     &-tc(3,3)*v2(3)-tc(4,3)*v2(4))
     &-MET(j,4)*(tc(1,4)*v2(1)-tc(2,4)*v2(2)
     &-tc(3,4)*v2(3)-tc(4,4)*v2(4))
     &+
     &MET(j,1)*(tc(1,1)*v2(1)-tc(1,2)*v2(2)
     &-tc(1,3)*v2(3)-tc(1,4)*v2(4))
     &-MET(j,2)*(tc(2,1)*v2(1)-tc(2,2)*v2(2)
     &-tc(2,3)*v2(3)-tc(2,4)*v2(4))
     &-MET(j,3)*(tc(3,1)*v2(1)-tc(3,2)*v2(2)
     &-tc(3,3)*v2(3)-tc(3,4)*v2(4))
     &-MET(j,4)*(tc(4,1)*v2(1)-tc(4,2)*v2(2)
     &-tc(4,3)*v2(3)-tc(4,4)*v2(4))


         TK12M(j) =
     &MET(j,1)*(tc(1,1)*p12(1)-tc(2,1)*p12(2)
     &-tc(3,1)*p12(3)-tc(4,1)*p12(4))
     &-MET(j,2)*(tc(1,2)*p12(1)-tc(2,2)*p12(2)
     &-tc(3,2)*p12(3)-tc(4,2)*p12(4))
     &-MET(j,3)*(tc(1,3)*p12(1)-tc(2,3)*p12(2)
     &-tc(3,3)*p12(3)-tc(4,3)*p12(4))
     &-MET(j,4)*(tc(1,4)*p12(1)-tc(2,4)*p12(2)
     &-tc(3,4)*p12(3)-tc(4,4)*p12(4))
     &+
     &MET(j,1)*(tc(1,1)*p12(1)-tc(1,2)*p12(2)
     &-tc(1,3)*p12(3)-tc(1,4)*p12(4))
     &-MET(j,2)*(tc(2,1)*p12(1)-tc(2,2)*p12(2)
     &-tc(2,3)*p12(3)-tc(2,4)*p12(4))
     &-MET(j,3)*(tc(3,1)*p12(1)-tc(3,2)*p12(2)
     &-tc(3,3)*p12(3)-tc(3,4)*p12(4))
     &-MET(j,4)*(tc(4,1)*p12(1)-tc(4,2)*p12(2)
     &-tc(4,3)*p12(3)-tc(4,4)*p12(4))
     
      enddo
        
      do i=1,4
         jvvt(i) = TV12*p12(i)+TKV1*v2(i)+TKV2*v1(i)
     &+(-p12(i)*V1V2-v2(i)*K2V1+v2(i)*K3V1
     &  +v1(i)*K1V2-v1(i)*K3V2)*T00
     &-K1V2*TV1M(i)+K3V2*TV1M(i)
     &+V1V2*TK12M(i)
     &+K2V1*TV2M(i)-K3V1*TV2M(i)
      enddo 

      if ( vmass.gt.rZero ) then
         do i=1,4
            jvvt(i) = jvvt(i)
     &+(K1V2*pv3(i)
     &-K3V2*TKV1*pv3(i)-K3V2*TK3V1*pv3(i)
     &-K1K3*TV12*pv3(i)+K2K3*TV12*pv3(i)
     &-V1V2*TK312*pv3(i)
     &-K2V1*TK3V2*pv3(i)+K3V1*TK3V2*pv3(i)
     &-K3V1*TKV2*pv3(i))/vmass**2
     &+(-pv3(i)*K1V2*K3V1+pv3(i)*K2V1*K3V2
     &+pv3(i)*V1V2*K1K3-pv3(i)*V1V2*K2K3
     &)/vmass**2*T00
         enddo
           
      endif	
      
      do i=1,4
         jvvt(i) = -jvvt(i) * d*g2
      enddo
      
      return
      end

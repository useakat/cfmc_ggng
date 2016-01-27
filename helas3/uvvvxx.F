      subroutine uvvvxx(v1,v2,v3,g,g2,tmass,twidth , ft)
c
c This subroutine computes an off-shell KK tensor boson wavefunction from
c the four-point coupling of
c three gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex v1(6)          : first  vector           flow-out V
c       complex v2(6)          : second vector           flow-out V
c       complex v3(6)          : third vector           flow-out V
c       complex g2             : coupling constant       -kappa/2
c       real g                 : coupling constant        (gs)
c       real tmass             : mass  of output KK tensor boson
c       real twidth            : width of output KK tensor boson

c
c output:
c       complex ft(18)         : off-shell tensor            
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex v1(6), v2(6), v3(6), yvvv(6,4)
      double complex g2,ft(18)
      double precision g 
      double precision tmass,twidth
      double precision MET(4,4)
      double complex d
      double complex V1V2,V1V3,V2V3,KTV1,KTV2,KTV3,
     &K12V3,K23V1,K31V2
      double precision KTK12,KTK31,KTK23 

      double precision pv1(4), pv2(4), pv3(4),pt(4),
     &pt2,p31(4),p23(4),p12(4)

      integer i, j

      double complex cZero
      double precision rZero, r2, r3,r4
      parameter( rZero = 0.0d0, r2 = 2.0d0, r3=3.d0,r4=4.d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      yvvv(5,1) = v1(5)+v2(5)+v3(5)
      yvvv(6,1) = v1(6)+v2(6)+v3(6)

      pv1(1) = dreal(v1(5))
      pv1(2) = dreal(v1(6))
      pv1(3) = dimag(v1(6))
      pv1(4) = dimag(v1(5))

      pv2(1) = dreal(v2(5))
      pv2(2) = dreal(v2(6))
      pv2(3) = dimag(v2(6))
      pv2(4) = dimag(v2(5))

      pv3(1) = dreal(v3(5))
      pv3(2) = dreal(v3(6))
      pv3(3) = dimag(v3(6))
      pv3(4) = dimag(v3(5))

      pt(1) = dreal(yvvv(5,1))
      pt(2) = dreal(yvvv(6,1))
      pt(3) = dimag(yvvv(6,1))
      pt(4) = dimag(yvvv(5,1))
	
      pt2=pt(1)**2-pt(2)**2-pt(3)**2-pt(4)**2
      if ( tmass.gt.rZero ) then
         d =  -g/dcmplx( pt2-tmass**2, tmass*twidth )
      else
         d =  -g/dcmplx( pt2, rZero )
      end if

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
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
      
	
      V1V2 =  v1(1)*v2(1) -  v1(2)*v2(2) -  v1(3)*v2(3) -  v1(4)*v2(4)
      V1V3 =  v1(1)*v3(1) -  v1(2)*v3(2) -  v1(3)*v3(3) -  v1(4)*v3(4)
      V2V3 =  v3(1)*v2(1) -  v3(2)*v2(2) -  v3(3)*v2(3) -  v3(4)*v2(4)
      K31V2 = p31(1)*v2(1) - p31(2)*v2(2) - p31(3)*v2(3) - p31(4)*v2(4)
      K12V3 = p12(1)*v3(1) - p12(2)*v3(2) - p12(3)*v3(3) - p12(4)*v3(4)
      K23V1 = p23(1)*v1(1) - p23(2)*v1(2) - p23(3)*v1(3) - p23(4)*v1(4)
      
      KTV1 = pt(1)*v1(1) - pt(2)*v1(2) - pt(3)*v1(3) - pt(4)*v1(4)
      KTV2 = pt(1)*v2(1) - pt(2)*v2(2) - pt(3)*v2(3) - pt(4)*v2(4)
      KTV3 = pt(1)*v3(1) - pt(2)*v3(2) - pt(3)*v3(3) - pt(4)*v3(4)
      
      KTK12 =pt(1)*p12(1)-pt(2)*p12(2)-pt(3)*p12(3)-pt(4)*p12(4)
      KTK31 =pt(1)*p31(1)-pt(2)*p31(2)-pt(3)*p31(3)-pt(4)*p31(4)
      KTK23 =pt(1)*p23(1)-pt(2)*p23(2)-pt(3)*p23(3)-pt(4)*p23(4)
      
      do i=1,4
         do j=1,4

            yvvv(i,j) = (r2*K12V3*KTV1*KTV2*r2*MET(i,j))/(r3*tmass**2)
     &	 + (r2*K31V2*KTV1*KTV3*r2*MET(i,j))/(r3*tmass**2) + 
     &  (r2*K23V1*KTV2*KTV3*r2*MET(i,j))/(r3*tmass**2) 
     &- r2*K12V3*V1V2*MET(i,j) + 
     &  (r2*KTK12*KTV3*r2*V1V2*MET(i,j))/(r3*tmass**2) 
     &- (K12V3*pt2*r2*V1V2*MET(i,j))/(r3*tmass**2) - 
     &  r2*K31V2*V1V3*MET(i,j)
     & + (r2*KTK31*KTV2*r2*V1V3*MET(i,j))/(r3*tmass**2) - 
     &  (K31V2*pt2*r2*V1V3*MET(i,j))/(r3*tmass**2) 
     &- r2*K23V1*V2V3*MET(i,j) + 
     &  (r2*KTK23*KTV1*r2*V2V3*MET(i,j))/(r3*tmass**2)
     & - (K23V1*pt2*r2*V2V3*MET(i,j))/(r3*tmass**2) - 
     &  (r2*KTV3*V1V2*p12(j)*pt(i))/tmass**2 
     &- (r2*KTV1*V2V3*p23(j)*pt(i))/tmass**2 
     &- (r2*KTV2*V1V3*p31(j)*pt(i))/tmass**2 - 
     &  (r2*KTV3*V1V2*p12(i)*pt(j))/tmass**2
     & - (r2*KTV1*V2V3*p23(i)*pt(j))/tmass**2 
     &- (r2*KTV2*V1V3*p31(i)*pt(j))/tmass**2 + 
     &  (r4*K12V3*KTV1*KTV2*pt(i)*pt(j))/tmass**4
     & + (r4*K31V2*KTV1*KTV3*pt(i)*pt(j))/tmass**4 + 
     &  (r4*K23V1*KTV2*KTV3*pt(i)*pt(j))/tmass**4 
     &- (r2*K12V3*KTV1*KTV2*r2*pt(i)*pt(j))/(r3*tmass**4) - 
     &  (r2*K31V2*KTV1*KTV3*r2*pt(i)*pt(j))/(r3*tmass**4)
     & - (r2*K23V1*KTV2*KTV3*r2*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (r4*KTK12*KTV3*V1V2*pt(i)*pt(j))/tmass**4 
     &- (r2*K12V3*pt2*V1V2*pt(i)*pt(j))/tmass**4 - 
     &  (r2*KTK12*KTV3*r2*V1V2*pt(i)*pt(j))/(r3*tmass**4)
     & + (K12V3*pt2*r2*V1V2*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (r4*K12V3*V1V2*pt(i)*pt(j))/tmass**2 
     &+ (r4*KTK31*KTV2*V1V3*pt(i)*pt(j))/tmass**4 - 
     &  (r2*K31V2*pt2*V1V3*pt(i)*pt(j))/tmass**4
     & - (r2*KTK31*KTV2*r2*V1V3*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (K31V2*pt2*r2*V1V3*pt(i)*pt(j))/(r3*tmass**4)
     & + (r4*K31V2*V1V3*pt(i)*pt(j))/tmass**2 + 
     &  (r4*KTK23*KTV1*V2V3*pt(i)*pt(j))/tmass**4
     & - (r2*K23V1*pt2*V2V3*pt(i)*pt(j))/tmass**4 - 
     &  (r2*KTK23*KTV1*r2*V2V3*pt(i)*pt(j))/(r3*tmass**4)
     & + (K23V1*pt2*r2*V2V3*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (r4*K23V1*V2V3*pt(i)*pt(j))/tmass**2 
     &+ r2*V2V3*p23(j)*v1(i) - (r2*K12V3*KTV2*pt(j)*v1(i))/tmass**2 - 
     &  (r2*K31V2*KTV3*pt(j)*v1(i))/tmass**2 
     &- (r2*KTK23*V2V3*pt(j)*v1(i))/tmass**2 + 2*V2V3*p23(i)*v1(j) - 
     &  (r2*K12V3*KTV2*pt(i)*v1(j))/tmass**2 
     &- (r2*K31V2*KTV3*pt(i)*v1(j))/tmass**2 
     &- (r2*KTK23*V2V3*pt(i)*v1(j))/tmass**2 + 
     &  r2*V1V3*p31(j)*v2(i) - (r2*K12V3*KTV1*pt(j)*v2(i))/tmass**2 
     &- (r2*K23V1*KTV3*pt(j)*v2(i))/tmass**2 - 
     &  (r2*KTK31*V1V3*pt(j)*v2(i))/tmass**2 
     &+ r2*K12V3*v1(j)*v2(i) + r2*V1V3*p31(i)*v2(j)
     & - (r2*K12V3*KTV1*pt(i)*v2(j))/tmass**2 - 
     &  (r2*K23V1*KTV3*pt(i)*v2(j))/tmass**2 
     &- (r2*KTK31*V1V3*pt(i)*v2(j))/tmass**2 
     &+ r2*K12V3*v1(i)*v2(j) + r2*V1V2*p12(j)*v3(i) - 
     &  (r2*K31V2*KTV1*pt(j)*v3(i))/tmass**2 
     &- (r2*K23V1*KTV2*pt(j)*v3(i))/tmass**2 
     &- (r2*KTK12*V1V2*pt(j)*v3(i))/tmass**2 + 
     &  r2*K31V2*v1(j)*v3(i) + r2*K23V1*v2(j)*v3(i) 
     &+ r2*V1V2*p12(i)*v3(j) - (r2*K31V2*KTV1*pt(i)*v3(j))/tmass**2 - 
     &  (r2*K23V1*KTV2*pt(i)*v3(j))/tmass**2
     & - (r2*KTK12*V1V2*pt(i)*v3(j))/tmass**2
     & + r2*K31V2*v1(i)*v3(j) + r2*K23V1*v2(i)*v3(j)
	
            yvvv(i,j) = -yvvv(i,j)*d/2.0d0*g2
	
         enddo
      enddo

      ft(1) = yvvv(1,1)
      ft(2) = yvvv(1,2)
      ft(3) = yvvv(1,3)
      ft(4) = yvvv(1,4)
      ft(5) = yvvv(2,1)
      ft(6) = yvvv(2,2)
      ft(7) = yvvv(2,3)
      ft(8) = yvvv(2,4)
      ft(9) = yvvv(3,1)
      ft(10) = yvvv(3,2)
      ft(11) = yvvv(3,3)
      ft(12) = yvvv(3,4)
      ft(13) = yvvv(4,1)
      ft(14) = yvvv(4,2)
      ft(15) = yvvv(4,3)
      ft(16) = yvvv(4,4)
      ft(17) = yvvv(5,1)
      ft(18) = yvvv(6,1)

      return
      end

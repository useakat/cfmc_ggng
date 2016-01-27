      subroutine ussxxx(sc1,sc2,g,smass,Tmass,Twidth , ft)
c
c This subroutine computes an off-shell KK Tensor wavfunction from 
c the scalar-scalar-Kaluza-Klein tensor boson coupling.
c
c input:
c       complex sc1(3)        : first input scalar               s1
c       complex sc2(3)        : second input scalar              s2
c       real    g             : coupling constant            -kappa/2
c       real    smass         : mass  of input scalar 
c       real    Tmass         : mass  of OUTPUT KK Tensor T'
c       real    Twidth        : width of OUTPUT KK Tensor T'
c output:
c       complex ft(18)        : KK tensor current       j(T:s1,s2)     
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex sc1(3),sc2(3),yss(6,4),ft(18)
      double precision g, smass,Tmass,Twidth
      integer i,j
      double precision ps1(4), ps2(4), pT(4)
      double precision MET(4,4)
      double complex cZero, d
      double precision rZero, rTwo
      double precision p1p2,pT2,p1pT,p2pT
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      yss(5,1) = sc1(2)+sc2(2)
      yss(6,1) = sc1(3)+sc2(3)
      
      ps1(1) = dreal(sc1(2))
      ps1(2) = dreal(sc1(3))
      ps1(3) = dimag(sc1(3))
      ps1(4) = dimag(sc1(2))
      
      ps2(1) = -dreal(sc2(2))
      ps2(2) = -dreal(sc2(3))
      ps2(3) = -dimag(sc2(3))
      ps2(4) = -dimag(sc2(2))
      
      pT(1) = dreal(yss(5,1))
      pT(2) = dreal(yss(6,1))
      pT(3) = dimag(yss(6,1))
      pT(4) = dimag(yss(5,1))
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      p1p2 = ps1(1)*ps2(1)-ps1(2)*ps2(2)-ps1(3)*ps2(3)-ps1(4)*ps2(4)
      p1pT = ps1(1)*pT(1)-ps1(2)*pT(2)-ps1(3)*pT(3)-ps1(4)*pT(4)
      p2pT = pT(1)*ps2(1)-pT(2)*ps2(2)-pT(3)*ps2(3)-pT(4)*ps2(4)
      pT2  = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
      
      if ( Tmass.gt.rZero ) then
         d = - g/dcmplx( pT2-Tmass**2, Tmass*Twidth )
      else
         d = - g/dcmplx( pT2, rZero )
      end if
    
      do i = 1,4
         do j=1,4
            yss(i,j) = -2.0d0/3.0d0*MET(i,j)*smass**2
     &	-4.0d0/3.0d0*smass**2/Tmass**2*pT(i)*pT(j)
     &    +2.0d0/3.0d0*smass**2/Tmass**2*PT2*MET(i,j)
     &    +4.0d0/3.0d0*smass**2/Tmass**4*PT2*pT(i)*pT(j)
     &    +2.0d0*(ps1(i)*ps2(j)+ps1(j)*ps2(i))
     &    -2.0d0/Tmass**2*p1pT*(ps2(i)*pT(j)+ps2(j)*pT(i)) 
     &    -2.0d0/Tmass**2*p2pT*(ps1(i)*pT(j)+ps1(j)*pT(i)) 
     &    +4.0d0/3.0d0/Tmass**2*MET(i,j)*p1pT*p2pT
     &    +8.0d0/3.0d0/Tmass**4*p1pT*p2pT*pT(i)*pT(j)
     &    -2.0d0/3.0d0*p1p2*MET(i,j)
     &    +8.0d0/3.0d0/Tmass**2*p1p2*pT(i)*pT(j)
     &    -2.0d0/3.0d0/Tmass**2*p1p2*pT2*MET(i,j)
     &    -4.0d0/3.0d0/Tmass**4*PT2*p1p2*pT(i)*pT(j)

            yss(i,j) = yss(i,j)*d*sc1(1)*sc2(1)/2.0d0

         end do
      enddo

      ft(1) = yss(1,1)
      ft(2) = yss(1,2)
      ft(3) = yss(1,3)
      ft(4) = yss(1,4)
      ft(5) = yss(2,1)
      ft(6) = yss(2,2)
      ft(7) = yss(2,3)
      ft(8) = yss(2,4)
      ft(9) = yss(3,1)
      ft(10) = yss(3,2)
      ft(11) = yss(3,3)
      ft(12) = yss(3,4)
      ft(13) = yss(4,1)
      ft(14) = yss(4,2)
      ft(15) = yss(4,3)
      ft(16) = yss(4,4)
      ft(17) = yss(5,1)
      ft(18) = yss(6,1)

      return
      end

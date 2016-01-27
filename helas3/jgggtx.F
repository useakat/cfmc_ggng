      subroutine jgggtx(gb,gc,gd,ft,g,g2 , gx)
c
c This subroutine computes an off-shell Tensor boson wavfunction from 
c the 4 gluons-Kaluza-Klein tensor boson coupling, corresponding 
c to the color structure f^{x,b,e} f{c,d,e}. 
c
c To optain the complete amplitude, this subroutine must be called three
c times (once for each color structure) with the following permutations:
c     call jgggtx(gb,gc,gd,ft,g,g2,gx)
c     call jgggtx(gc,gd,gb,ft,g,g2,gx)
c     call jgggtx(gd,gb,gc,ft,g,g2,gx)
c     corresponding to 
c	f^{x,b,e} f{c,d,e}
c	f^{x,c,e} f{d,b,e}
c	f^{x,d,e} f{b,c,e}
c
c input:
c       complex gb(6)          : first input gluon          gb
c       complex gc(6)          : second input gluon         gc
c       complex gd(6)          : third input gluon          gd
c       complex ft(18)         : graviton KK tensor         ft
c
c       real    g              : coupling constant          gs
c       complex g2             : coupling constant        -k/2
c output:
c       complex gx(6)          : off-shell gluon current   j(T:gb,gc,gd,ft) 
c
c- by Q.Li - Jan. 2008
c     
      implicit none
      double complex gx(6),gb(6),gc(6),gd(6),tc(6,4),ft(18)
      double complex TVBVC, TVBVD, TVDM(4), TVCM(4)
      double complex T00,T12,T13,T14,T23,T24,T34,EBEC,EBED
      double precision g

      integer a,b,i,j

      double precision pb(4),pc(4),pd(4),pT(4),px(4)

      double precision MET(4,4)
      double complex cZero, d, g2
      double precision rZero, rTwo,px2,r3
      parameter( rZero = 0.0d0, rTwo = 2.0d0, r3=3.0d0)
      parameter(cZero=(0.0d0,0.0d0))


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

      gx(5)=gb(5)+gc(5)+gd(5)+tc(5,1)
      gx(6)=tc(6,1)+gb(6)+gc(6)+gd(6)

      px(1) = dreal(gx(5))
      px(2) = dreal(gx(6))
      px(3) = dimag(gx(6))
      px(4) = dimag(gx(5))

      pb(1) = dreal(gb(5))
      pb(2) = dreal(gb(6))
      pb(3) = dimag(gb(6))
      pb(4) = dimag(gb(5))

      pc(1) = dreal(gc(5))
      pc(2) = dreal(gc(6))
      pc(3) = dimag(gc(6))
      pc(4) = dimag(gc(5))

      pd(1) = dreal(gd(5))
      pd(2) = dreal(gd(6))
      pd(3) = dimag(gd(6))
      pd(4) = dimag(gd(5))

     
      do i=1,4
         do j=1,4
            MET(i,j)=0.0d0
         enddo 
      enddo

      MET(1,1)=1.0d0
      MET(2,2)=-1.0d0
      MET(3,3)=-1.0d0
      MET(4,4)=-1.0d0

      px2 = px(1)**2-px(2)**2-px(3)**2-px(4)**2
      EBEC = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      EBED = gb(1)*gd(1)-gb(2)*gd(2)-gb(3)*gd(3)-gb(4)*gd(4)
      
      d = 1.0d0/dcmplx( px2, 0.0d0 )
      
      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      TVBVC =2.0d0*(tc(1,1)*gb(1)*gc(1)
     &+tc(2,2)*gb(2)*gc(2)+tc(3,3)*gb(3)*gc(3)
     &+tc(4,4)*gb(4)*gc(4)
     &)
     & - T12*(gb(1)*gc(2) + gb(2)*gc(1))
     &          - T13*(gb(1)*gc(3) + gb(3)*gc(1))
     &          - T14*(gb(1)*gc(4) + gb(4)*gc(1))
     &          + T23*(gb(2)*gc(3) + gb(3)*gc(2))
     &          + T24*(gb(2)*gc(4) + gb(4)*gc(2))
     &          + T34*(gb(3)*gc(4) + gb(4)*gc(3))

      TVBVD =2.0d0*(tc(1,1)*gb(1)*gd(1)
     &+tc(2,2)*gb(2)*gd(2)+tc(3,3)*gb(3)*gd(3)
     &+tc(4,4)*gb(4)*gd(4)
     &)
     & - T12*(gb(1)*gd(2) + gb(2)*gd(1))
     &          - T13*(gb(1)*gd(3) + gb(3)*gd(1))
     &          - T14*(gb(1)*gd(4) + gb(4)*gd(1))
     &          + T23*(gb(2)*gd(3) + gb(3)*gd(2))
     &          + T24*(gb(2)*gd(4) + gb(4)*gd(2))
     &          + T34*(gb(3)*gd(4) + gb(4)*gd(3))

      do j=1,4
         TVDM(j)=
     &MET(j,1)*(tc(1,1)*gd(1)-tc(2,1)*gd(2)
     &-tc(3,1)*gd(3)-tc(4,1)*gd(4))
     &-MET(j,2)*(tc(1,2)*gd(1)-tc(2,2)*gd(2)
     &-tc(3,2)*gd(3)-tc(4,2)*gd(4))
     &-MET(j,3)*(tc(1,3)*gd(1)-tc(2,3)*gd(2)
     &-tc(3,3)*gd(3)-tc(4,3)*gd(4))
     &-MET(j,4)*(tc(1,4)*gd(1)-tc(2,4)*gd(2)
     &-tc(3,4)*gd(3)-tc(4,4)*gd(4))
     &+
     &MET(j,1)*(tc(1,1)*gd(1)-tc(1,2)*gd(2)
     &-tc(1,3)*gd(3)-tc(1,4)*gd(4))
     &-MET(j,2)*(tc(2,1)*gd(1)-tc(2,2)*gd(2)
     &-tc(2,3)*gd(3)-tc(2,4)*gd(4))
     &-MET(j,3)*(tc(3,1)*gd(1)-tc(3,2)*gd(2)
     &-tc(3,3)*gd(3)-tc(3,4)*gd(4))
     &-MET(j,4)*(tc(4,1)*gd(1)-tc(4,2)*gd(2)
     &-tc(4,3)*gd(3)-tc(4,4)*gd(4))
	
	TVCM(j)=
     &MET(j,1)*(tc(1,1)*gc(1)-tc(2,1)*gc(2)
     &-tc(3,1)*gc(3)-tc(4,1)*gc(4))
     &-MET(j,2)*(tc(1,2)*gc(1)-tc(2,2)*gc(2)
     &-tc(3,2)*gc(3)-tc(4,2)*gc(4))
     &-MET(j,3)*(tc(1,3)*gc(1)-tc(2,3)*gc(2)
     &-tc(3,3)*gc(3)-tc(4,3)*gc(4))
     &-MET(j,4)*(tc(1,4)*gc(1)-tc(2,4)*gc(2)
     &-tc(3,4)*gc(3)-tc(4,4)*gc(4))
     &+
     &MET(j,1)*(tc(1,1)*gc(1)-tc(1,2)*gc(2)
     &-tc(1,3)*gc(3)-tc(1,4)*gc(4))
     &-MET(j,2)*(tc(2,1)*gc(1)-tc(2,2)*gc(2)
     &-tc(2,3)*gc(3)-tc(2,4)*gc(4))
     &-MET(j,3)*(tc(3,1)*gc(1)-tc(3,2)*gc(2)
     &-tc(3,3)*gc(3)-tc(3,4)*gc(4))
     &-MET(j,4)*(tc(4,1)*gc(1)-tc(4,2)*gc(2)
     &-tc(4,3)*gc(3)-tc(4,4)*gc(4))

      enddo
	    
      do a=1,4
	
         gx(a)=gd(a)*TVBVC-gc(a)*TVBVD
     &        -T00*gd(a)*EBEC+T00*gc(a)*EBED
     &        +EBEC*TVDM(a)-EBED*TVCM(a)
	 
         gx(a)=-gx(a)*d*g*g*g2
         
      enddo 
 
      return
      end

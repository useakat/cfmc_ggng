      subroutine ggggtx(ga,gb,gc,gd,ft,g,g2 , vertex)
c      
c This subroutine computes the portion of the amplitude of the five-point 
c coupling of a KK tensor boson with 4 massless color octet gauge bosons
c (gluons) corresponding 
c to the color structure f^{a,b,e} f{c,d,e}. 
c To optain the complete amplitude, this coupling must be called three
c times (once for each color structure) with the following permutations:
c     call ggggtx(ga,gb,gc,gd,g,v1)
c     call ggggtx(ga,gc,gd,gb,g,v2)
c     call ggggtx(ga,gd,gb,gc,g,v3)
c
c	f^{a,b,e} f{c,d,e}
c	f^{a,c,e} f{d,b,e}
c	f^{a,d,e} f{b,c,e}
c                                                                       
c input:                                                                
c       complex ga(6)        : Boson with adjoint color index a 
c       complex gb(6)        : Boson with adjoint color index b
c       complex gc(6)        : Boson with adjoint color index c 
c       complex gd(6)        : Boson with adjoint color index d
c	complex ft(18)       : KK tensor boson
c       real    g            : gs            
c       complex g2           : -k/2            
c
c output:
c       complex vertex       : amplitude
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex ga(6),gb(6),gc(6),gd(6),tc(6,4),vertex,dvertx
      double complex vab,vac,vad,vbc,vbd,vcd,ft(18)
  
      double complex T00, T12, T13, T14, T23, T24, T34	
      double complex TV24,TV23,TV14,TV13
      double complex g2
      double complex cZero
      double precision rZero, rTwo,g
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

      T00 = tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      TV14 = rtwo*(tc(1,1)*ga(1)*gd(1)+tc(2,2)*ga(2)*gd(2)
     &+tc(3,3)*ga(3)*gd(3)+tc(4,4)*ga(4)*gd(4))

      TV13 = rtwo*(tc(1,1)*ga(1)*gc(1)+tc(2,2)*ga(2)*gc(2)
     &+tc(3,3)*ga(3)*gc(3)+tc(4,4)*ga(4)*gc(4))

      TV24 = rtwo*(tc(1,1)*gb(1)*gd(1)+tc(2,2)*gb(2)*gd(2)
     &+tc(3,3)*gb(3)*gd(3)+tc(4,4)*gb(4)*gd(4))

      TV23 = rtwo*(tc(1,1)*gb(1)*gc(1)+tc(2,2)*gb(2)*gc(2)
     &+tc(3,3)*gb(3)*gc(3)+tc(4,4)*gb(4)*gc(4))

	
      TV14 = TV14- T12*(ga(1)*gd(2) + ga(2)*gd(1))
     &          - T13*(ga(1)*gd(3) + ga(3)*gd(1))
     &          - T14*(ga(1)*gd(4) + ga(4)*gd(1))
     &          + T23*(ga(2)*gd(3) + ga(3)*gd(2))
     &          + T24*(ga(2)*gd(4) + ga(4)*gd(2))
     &          + T34*(ga(3)*gd(4) + ga(4)*gd(3))
      
      TV13 = TV13 - T12*(ga(1)*gc(2) + ga(2)*gc(1))
     &          - T13*(ga(1)*gc(3) + ga(3)*gc(1))
     &          - T14*(ga(1)*gc(4) + ga(4)*gc(1))
     &          + T23*(ga(2)*gc(3) + ga(3)*gc(2))
     &          + T24*(ga(2)*gc(4) + ga(4)*gc(2))
     &          + T34*(ga(3)*gc(4) + ga(4)*gc(3))

      TV24 = TV24 - T12*(gb(1)*gd(2) + gb(2)*gd(1))
     &          - T13*(gb(1)*gd(3) + gb(3)*gd(1))
     &          - T14*(gb(1)*gd(4) + gb(4)*gd(1))
     &          + T23*(gb(2)*gd(3) + gb(3)*gd(2))
     &          + T24*(gb(2)*gd(4) + gb(4)*gd(2))
     &          + T34*(gb(3)*gd(4) + gb(4)*gd(3))

      TV23 = TV23 - T12*(gb(1)*gc(2) + gb(2)*gc(1))
     &          - T13*(gb(1)*gc(3) + gb(3)*gc(1))
     &          - T14*(gb(1)*gc(4) + gb(4)*gc(1))
     &          + T23*(gb(2)*gc(3) + gb(3)*gc(2))
     &          + T24*(gb(2)*gc(4) + gb(4)*gc(2))
     &          + T34*(gb(3)*gc(4) + gb(4)*gc(3))
     	

      vab = ga(1)*gb(1)-ga(2)*gb(2)-ga(3)*gb(3)-ga(4)*gb(4)
      vac = ga(1)*gc(1)-ga(2)*gc(2)-ga(3)*gc(3)-ga(4)*gc(4)
      vad = ga(1)*gd(1)-ga(2)*gd(2)-ga(3)*gd(3)-ga(4)*gd(4)
      vbc = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      vbd = gb(1)*gd(1)-gb(2)*gd(2)-gb(3)*gd(3)-gb(4)*gd(4)
      vcd = gc(1)*gd(1)-gc(2)*gd(2)-gc(3)*gd(3)-gc(4)*gd(4)

      dvertx = -TV13*vbd-TV24*vac+TV23*vad+TV14*vbc 
     &+vbd*vac*T00-vad*vbc*T00

      vertex = -dvertx * g*g*g2

      return
      end

      subroutine ggggcf(g1,g2,g3,g4,g, vertex)
c
c this subroutine computes an amplitude of the four-point coupling of   
c gluoins. 
c                                         
c input:                                                                
c       complex g1(0:3)        : first gluon                         g1 
c       complex g2(0:3)        : second gluon                        g2
c       complex g3(0:3)        : third gluon                         g3 
c       complex g4(0:3)        : fourth gluon                        g4
c       real    g              : coupling constant            
c                                                                       
c output:                                                               
c       complex vertex         : amplitude          gamma(g1,g2,g3,g4)
c
      implicit none
      double complex g1(6),g2(6),g3(6),g4(6),vertex,g
      double complex dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),
     &           dvertx,v12,v13,v14,v23,v24,v34

      double precision r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )

      dv1(0)=dcmplx(g1(1))
      dv1(1)=dcmplx(g1(2))
      dv1(2)=dcmplx(g1(3))
      dv1(3)=dcmplx(g1(4))
      dv2(0)=dcmplx(g2(1))
      dv2(1)=dcmplx(g2(2))
      dv2(2)=dcmplx(g2(3))
      dv2(3)=dcmplx(g2(4))
      dv3(0)=dcmplx(g3(1))
      dv3(1)=dcmplx(g3(2))
      dv3(2)=dcmplx(g3(3))
      dv3(3)=dcmplx(g3(4))
      dv4(0)=dcmplx(g4(1)) 
      dv4(1)=dcmplx(g4(2))
      dv4(2)=dcmplx(g4(3))
      dv4(3)=dcmplx(g4(4))

      v12= dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13= dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14= dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23= dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24= dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34= dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

      dvertx =( -v14*v23 + 2.d0*v13*v24 - v12*v34)

      vertex = dcmplx( dvertx ) * (g*g)

      return
      end


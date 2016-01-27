      subroutine jgggcf(w1,w2,w3,g, jggg)
c
c This subroutine computes an off-shell gluon current
c from the four-point gluon coupling.
c
c input:
c       complex w1(6)          : first  gluon                        w1
c       complex w2(6)          : second gluon                        w2
c       complex w3(6)          : third  gluon                        w3
c       real    g              : coupling constant
c
c output:
c       complex jggg(6)        : gluon current             j^mu(w':w1,w2,w3)
c
      implicit none
      double complex w1(6),w2(6),w3(6),jggg(6)
      double complex dw1(0:3),dw2(0:3),dw3(0:3)
      double complex jj(0:3),dv,w32,w13,w21,dg2,g
      double precision p1(0:3),p2(0:3),p3(0:3),q(0:3),q2

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )

      jggg(5) = w1(5)+w2(5)+w3(5)
      jggg(6) = w1(6)+w2(6)+w3(6)

      dw1(0) = dcmplx(w1(1))
      dw1(1) = dcmplx(w1(2))
      dw1(2) = dcmplx(w1(3))
      dw1(3) = dcmplx(w1(4))
      dw2(0) = dcmplx(w2(1))
      dw2(1) = dcmplx(w2(2))
      dw2(2) = dcmplx(w2(3))
      dw2(3) = dcmplx(w2(4))
      dw3(0) = dcmplx(w3(1))
      dw3(1) = dcmplx(w3(2))
      dw3(2) = dcmplx(w3(3))
      dw3(3) = dcmplx(w3(4))
      q(0) = -dble(jggg(5))
      q(1) = -dble(jggg(6))
      q(2) = -dimag(jggg(6))
      q(3) = -dimag(jggg(5))

      q2 = q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)

      dg2 = g*g

      dv = rOne/dcmplx( q2 )

      w32 = dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)
      w13 = dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)
      w21 = dw2(0)*dw1(0)-dw2(1)*dw1(1)-dw2(2)*dw1(2)-dw2(3)*dw1(3)

      jj(0) = dg2*( -dw1(0)*w32 + 2d0*dw2(0)*w13 - dw3(0) *w21 )
      jj(1) = dg2*( -dw1(1)*w32 + 2d0*dw2(1)*w13 - dw3(1) *w21 )
      jj(2) = dg2*( -dw1(2)*w32 + 2d0*dw2(2)*w13 - dw3(2) *w21 )
      jj(3) = dg2*( -dw1(3)*w32 + 2d0*dw2(3)*w13 - dw3(3) *w21 )

      jggg(1) = dcmplx( jj(0)*dv )
      jggg(2) = dcmplx( jj(1)*dv )
      jggg(3) = dcmplx( jj(2)*dv )
      jggg(4) = dcmplx( jj(3)*dv )
c
      return
      end

      subroutine GLUON4_0(w1,w2,w3,w4,g,vertex)

      implicit none

      integer i
      double complex w1(18),w2(18),w3(18),w4(18)
      double complex z,w12(18),w23(18),vertex
      double complex g,gg,ci
      parameter(ci=(0d0,1d0))

      vertex=(0d0,0d0)
      gg = -ci*g

      call VVV1_1(w1,w2, gg,0d0,0d0, w12)
      call VVV1_1(w2,w3, gg,0d0,0d0, w23)

      call VVV1_0(w12,w3,w4,gg,z)
      vertex = vertex+z
      call VVV1_0(w1,w23,w4,gg,z)
      vertex = vertex+z
      call ggggcf(w1,w2,w3,w4,g,z)
      vertex = vertex+z
    
      return
      end

      subroutine gluon4_1(w1,w2,w3,g,mass,width, jgluon4)
 
      implicit none

      integer i
      double precision mass, width
      double complex w(6,12,12),w1(6),w2(6),w3(6)
      double complex wx(6,55),jgluon4(6)     
      double complex g,gg,ci
      parameter(ci=(0d0,1d0))

       gg = -ci*g

      do i=1,6

         w(i,1,1) = w1(i)
         w(i,2,2) = w2(i)
         w(i,3,3) = w3(i)
      enddo

      do i=1,2
         call VVV1_2(w(1,i,i),w(1,i+1,i+1), gg,0d0,0d0,w(1,i,i+1))     
      enddo

      call VVV1_2(w(1,1,2),w(1,3,3),gg,0d0,0d0,wx(1,1))              
      call VVV1_2(w(1,1,1)  ,w(1,2,3),gg,0d0,0d0,wx(1,2))              
      call jgggcf(w(1,1,1)  ,w(1,2,2),w(1,3,3), g,wx(1,3))
      call sumw(wx,3,jgluon4)     
      
      return            
      end


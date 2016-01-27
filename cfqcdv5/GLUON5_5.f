      subroutine gluon5_5(w1,w2,w3,w4,g,mass,width, jgluon5) 

      implicit none

      integer i
      double precision mass,width
      double complex w(6,12,12),w1(6),w2(6),w3(6),w4(6)
      double complex wx(6,55),jgluon5(6)   
      double complex g,gg,ci
      parameter(ci=(0d0,1d0))

       gg = -ci*g
      
      do i=1,6
         jgluon5(i) = (0d0,0d0)
         w(i,1,1) = w1(i)
         w(i,2,2) = w2(i)
         w(i,3,3) = w3(i)
         w(i,4,4) = w4(i)
      enddo
     
***************************************************************************

      do i=1,3
         call VVV1_2(w(1,i,i),w(1,i+1,i+1),gg,0d0,0d0,w(1,i,i+1))     
      enddo

***************************************************************************

       do i=1,2 
          call VVV1_2(w(1,i,i+1),w(1,i+2,i+2),gg,0d0,0d0,wx(1,1))              
          call VVV1_2(w(1,i,i),w(1,i+1,i+2),gg,0d0,0d0,wx(1,2))              
          call jgggcf(w(1,i,i),w(1,i+1,i+1),w(1,i+2,i+2),g,wx(1,3))             
          call sumw(wx,3,w(1,i,i+2))
       enddo
      
***************************************************************************
 
       call VVV1_2(w(1,1,3),w(1,4,4),gg,0d0,0d0,wx(1,1))                
       call VVV1_2(w(1,1,2),w(1,3,4),gg,0d0,0d0,wx(1,2))                
       call VVV1_2(w(1,1,1),w(1,2,4),gg,0d0,0d0,wx(1,3))                
       call jgggcf(w(1,1,2),w(1,3,3),w(1,4,4),g,wx(1,4))               
       call jgggcf(w(1,1,1),w(1,2,3),w(1,4,4),g,wx(1,5))               
       call jgggcf(w(1,1,1),w(1,2,2),w(1,3,4),g,wx(1,6))               
       call sumw(wx,6,jgluon5)
 
      return
      end
      

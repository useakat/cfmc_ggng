      subroutine GLUON5_0(w1,w2,w3,w4,w5,g,vertex) 

      implicit none    

      integer i
      double complex w(18,12,12),w1(18),w2(18),w3(18),w4(18),w5(18)
      double complex z(55),wx(18,55),vertex    
      double complex g,gg,ci
      parameter(ci=(0d0,1d0))
      
      vertex=(0d0,0d0)
      gg = -ci*g

      do i=1,18
         w(i,1,1) = w1(i)
         w(i,2,2) = w2(i)
         w(i,3,3) = w3(i)
         w(i,4,4) = w4(i)
         w(i,5,5) = w5(i)
      enddo

***************************************************************************

      do i=1,3
         call VVV1_1(w(1,i,i),w(1,i+1,i+1), gg,0d0,0d0,w(1,i,i+1))     
      enddo

***************************************************************************
         
       do i=1,2  
          call VVV1_1(w(1,i,i+1),w(1,i+2,i+2),gg, 0d0, 0d0, wx(1,1))              
          call VVV1_1(w(1,i,i)  ,w(1,i+1,i+2),gg, 0d0, 0d0, wx(1,2))              
          call jgggcf(w(1,i,i)  ,w(1,i+1,i+1),w(1,i+2,i+2), g,wx(1,3))             
          call sumw(wx,3,w(1,i,i+2))
       enddo   
      
***************************************************************************
      
       call VVV1_0(w(1,1,3),w(1,4,4),w(1,5,5),gg,z(1))      
       call VVV1_0(w(1,1,2),w(1,3,4),w(1,5,5),gg,z(2))      
       call VVV1_0(w(1,1,1),w(1,2,4),w(1,5,5),gg,z(3))      
       call ggggcf(w(1,1,2),w(1,3,3),w(1,4,4),w(1,5,5),g,z(4))     
       call ggggcf(w(1,1,1),w(1,2,3),w(1,4,4),w(1,5,5),g,z(5))     
       call ggggcf(w(1,1,1),w(1,2,2),w(1,3,4),w(1,5,5),g,z(6))     
       do i=1,6
          vertex = vertex+z(i)
       enddo   
       
       return 
       end



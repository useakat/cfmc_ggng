      subroutine gluon10_0(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,g,vertex)

      implicit none
     
      integer i,k, ngluons
      double complex w(18,12,12),w1(18),w2(18),w3(18),w4(18),w5(18)
      double complex w6(18),w7(18),w8(18),w9(18),w10(18)
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
         w(i,6,6) = w6(i)
         w(i,7,7) = w7(i)
         w(i,8,8) = w8(i)
         w(i,9,9) = w9(i)
         w(i,10,10) = w10(i)
      enddo

      ngluons = 10
      
      DO i=1,ngluons-2
         CALL VVV1_1(W(1,I,I),W(1,I+1,I+1), GG,0d0,0d0,W(1,I,I+1))
      ENDDO

***************************************************************************

       DO i=1,ngluons-3
          CALL VVV1_1(W(1,I,I+1),W(1,I+2,I+2),GG,0d0,0d0,WX(1,1))              
          CALL VVV1_1(W(1,I,I)  ,W(1,I+1,I+2),GG,0d0,0d0,WX(1,2))              
          CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+2), G,WX(1,3))              
          CALL SUMW(WX,3,W(1,I,I+2))
       ENDDO
      
***************************************************************************

      DO i=1,ngluons-4
         CALL VVV1_1(W(1,I,I+2),W(1,I+3,I+3),GG,0d0,0d0,WX(1,1))                
         CALL VVV1_1(W(1,I,I+1),W(1,I+2,I+3),GG,0d0,0d0,WX(1,2))                
         CALL VVV1_1(W(1,I,I)  ,W(1,I+1,I+3),GG,0d0,0d0,WX(1,3))                
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+3),G,WX(1,4))                
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+3),G,WX(1,5))                
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+3),G,WX(1,6))                
         CALL SUMW(WX,6,W(1,I,I+3))
      ENDDO
      
***************************************************************************

      DO i=1,ngluons-5
C 3-point
         CALL VVV1_1(W(1,I,I+3),W(1,I+4,I+4),GG,0d0,0d0,WX(1, 1))               
         CALL VVV1_1(W(1,I,I+2),W(1,I+3,I+4),GG,0d0,0d0,WX(1, 2))               
         CALL VVV1_1(W(1,I,I+1),W(1,I+2,I+4),GG,0d0,0d0,WX(1, 3))               
         CALL VVV1_1(W(1,I,I)  ,W(1,I+1,I+4),GG,0d0,0d0,WX(1, 4))               
C 4-point
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+3), W(1,I+4,I+4), G,WX(1, 5))             
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2), W(1,I+3,I+4), G,WX(1, 6))             
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+3), W(1,I+4,I+4), G,WX(1, 7))             
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1), W(1,I+2,I+4), G,WX(1, 8))             
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2), W(1,I+3,I+4), G,WX(1, 9))             
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+3), W(1,I+4,I+4), G,WX(1,10))             
         CALL SUMW(WX,10,W(1,I,I+4))
      ENDDO
      

***************************************************************************

      DO i=1,ngluons-6  
C 3-point
         CALL VVV1_1(W(1,I,I+4),W(1,I+5,I+5),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_1(W(1,I,I+3),W(1,I+4,I+5),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_1(W(1,I,I+2),W(1,I+3,I+5),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_1(W(1,I,I+1),W(1,I+2,I+5),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_1(W(1,I,I)  ,W(1,I+1,I+5),GG,0d0,0d0,WX(1, 5))              
C 4-point
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+4),  W(1,I+5,I+5), G,WX(1, 6))            
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+4),  W(1,I+5,I+5), G,WX(1, 7))            
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+4),  W(1,I+5,I+5), G,WX(1, 8))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+4),  W(1,I+5,I+5), G,WX(1, 9))            
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+3),  W(1,I+4,I+5), G,WX(1,10))            
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+3),  W(1,I+4,I+5), G,WX(1,11))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+3),  W(1,I+4,I+5), G,WX(1,12))            
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),  W(1,I+3,I+5), G,WX(1,13))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),  W(1,I+3,I+5), G,WX(1,14))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),  W(1,I+2,I+5), G,WX(1,15))            
         CALL SUMW(WX,15,W(1,I,I+5))
      ENDDO
     
***************************************************************************
 
      DO i=1,ngluons-7 
         
C 3-point
         CALL VVV1_1(W(1,I,I+5),W(1,I+6,I+6),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_1(W(1,I,I+4),W(1,I+5,I+6),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_1(W(1,I,I+3),W(1,I+4,I+6),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_1(W(1,I,I+2),W(1,I+3,I+6),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_1(W(1,I,I+1),W(1,I+2,I+6),GG,0d0,0d0,WX(1, 5))              
         CALL VVV1_1(W(1,I,I)  ,W(1,I+1,I+6),GG,0d0,0d0,WX(1, 6))              

C 4-point
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+5),  W(1,I+6,I+6), G,WX(1, 7))            
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+5),  W(1,I+6,I+6), G,WX(1, 8))            
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+5),  W(1,I+6,I+6), G,WX(1, 9))            
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+5),  W(1,I+6,I+6), G,WX(1,10))            
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+5),  W(1,I+6,I+6), G,WX(1,11))            
         
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+4),  W(1,I+5,I+6), G,WX(1,12))            
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+4),  W(1,I+5,I+6), G,WX(1,13))            
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+4),  W(1,I+5,I+6), G,WX(1,14))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+4),  W(1,I+5,I+6), G,WX(1,15))            
         
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+3),  W(1,I+4,I+6), G,WX(1,16))            
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+3),  W(1,I+4,I+6), G,WX(1,17))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+3),  W(1,I+4,I+6), G,WX(1,18))            
         
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),  W(1,I+3,I+6), G,WX(1,19))            
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),  W(1,I+3,I+6), G,WX(1,20))            
         
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),  W(1,I+2,I+6), G,WX(1,21))            
         
         CALL SUMW(WX,21,W(1,I,I+6))
      ENDDO
      


***************************************************************************
 
      DO i=1,ngluons-8 

C 3-point
         CALL VVV1_1(W(1,I,I+6),W(1,I+7,I+7),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_1(W(1,I,I+5),W(1,I+6,I+7),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_1(W(1,I,I+4),W(1,I+5,I+7),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_1(W(1,I,I+3),W(1,I+4,I+7),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_1(W(1,I,I+2),W(1,I+3,I+7),GG,0d0,0d0,WX(1, 5))              
         CALL VVV1_1(W(1,I,I+1),W(1,I+2,I+7),GG,0d0,0d0,WX(1, 6))              
         CALL VVV1_1(W(1,I,I)  ,W(1,I+1,I+7),GG,0d0,0d0,WX(1, 7))              

C 4-point
         CALL JGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+7), G,WX(1, 8))              
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+7), G,WX(1, 9))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+7), G,WX(1,10))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+7), G,WX(1,11))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+7), G,WX(1,12))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+7), G,WX(1,13))              
         
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+7), G,WX(1,14))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+7), G,WX(1,15))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+7), G,WX(1,16))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+7), G,WX(1,17))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+5),W(1,I+6,I+7), G,WX(1,18))              
         
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+7), G,WX(1,19))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+7), G,WX(1,20))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+7), G,WX(1,21))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+7), G,WX(1,22))              
         
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+7), G,WX(1,23))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+7), G,WX(1,24))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+7), G,WX(1,25))              
         
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+7), G,WX(1,26))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+7), G,WX(1,27))              
         
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+7), G,WX(1,28))              
         
         CALL SUMW(WX,28,W(1,I,I+7))
      ENDDO
      


***************************************************************************
      
         i=1
      
C 3-point
         CALL VVV1_0(W(1,I,I+7),W(1,I+8,I+8),W(1,I+9,I+9),GG,Z( 1))     
         CALL VVV1_0(W(1,I,I+6),W(1,I+7,I+8),W(1,I+9,I+9),GG,Z( 2))     
         CALL VVV1_0(W(1,I,I+5),W(1,I+6,I+8),W(1,I+9,I+9),GG,Z( 3))     
         CALL VVV1_0(W(1,I,I+4),W(1,I+5,I+8),W(1,I+9,I+9),GG,Z( 4))     
         CALL VVV1_0(W(1,I,I+3),W(1,I+4,I+8),W(1,I+9,I+9),GG,Z( 5))     
         CALL VVV1_0(W(1,I,I+2),W(1,I+3,I+8),W(1,I+9,I+9),GG,Z( 6))     
         CALL VVV1_0(W(1,I,I+1),W(1,I+2,I+8),W(1,I+9,I+9),GG,Z( 7))     
         CALL VVV1_0(W(1,I,I)  ,W(1,I+1,I+8),W(1,I+9,I+9),GG,Z( 8))     

C 4-point
         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z( 9))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z(10))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z(11))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z(12))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z(13))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z(14))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+7),W(1,I+8,I+8),W(1,I+9,I+9),
     &        G,Z(15))     
         
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+8),W(1,I+9,I+9),
     &        G,Z(16))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+8),W(1,I+9,I+9),
     &        G,Z(17))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+8),W(1,I+9,I+9),
     &        G,Z(18))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+8),W(1,I+9,I+9),
     &        G,Z(19))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+8),W(1,I+9,I+9),
     &        G,Z(20))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+8),W(1,I+9,I+9),
     &        G,Z(21))     
         
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+8),W(1,I+9,I+9),
     &        G,Z(22))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+8),W(1,I+9,I+9),
     &        G,Z(23))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+8),W(1,I+9,I+9),
     &        G,Z(24))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+8),W(1,I+9,I+9),
     &        G,Z(25))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+5),W(1,I+6,I+8),W(1,I+9,I+9),
     &        G,Z(26))     
         
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+8),W(1,I+9,I+9),
     &        G,Z(27))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+8),W(1,I+9,I+9),
     &        G,Z(28))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+8),W(1,I+9,I+9),
     &        G,Z(29))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+8),W(1,I+9,I+9),
     &        G,Z(30))     
         
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+8),W(1,I+9,I+9),
     &        G,Z(31))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+8),W(1,I+9,I+9),
     &        G,Z(32))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+8),W(1,I+9,I+9),
     &        G,Z(33))     
         
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+8),W(1,I+9,I+9),
     &        G,Z(34))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+8),W(1,I+9,I+9),
     &        G,Z(35))     
         
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+8),W(1,I+9,I+9),
     &        G,Z(36))     
         
         DO K=1,36
            vertex=vertex+Z(K)
         ENDDO  
       
      END

      subroutine gluon8_8(w1,w2,w3,w4,w5,w6,w7,g,mass,width, jgluo8)

      implicit none
     
      integer i, ngluons
      double precision mass,width
      double complex w(6,12,12),w1(6),w2(6),w3(6),w4(6),w5(6)
      double complex w6(6),w7(6),w8(6)
      double complex z(55),wx(6,55),jgluo8(6)
      double complex g,gg,ci
      parameter(ci=(0d0,1d0))
      
      gg = -ci*g

       do i = 1,6
         jgluo8(i) = (0d0,0d0)
      enddo

      do i=1,6
         w(i,1,1) = w1(i)
         w(i,2,2) = w2(i)
         w(i,3,3) = w3(i)
         w(i,4,4) = w4(i)
         w(i,5,5) = w5(i)
         w(i,6,6) = w6(i)
         w(i,7,7) = w7(i)
         w(i,8,8) = w8(i)
      enddo

      ngluons = 8
      
      DO i=1,ngluons-2
         CALL VVV1_2(W(1,I,I),W(1,I+1,I+1), GG,0d0,0d0,W(1,I,I+1))
      ENDDO

***************************************************************************

       DO i=1,ngluons-3
          CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+2),GG,0d0,0d0,WX(1,1))              
          CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+2),GG,0d0,0d0,WX(1,2))              
          CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+2), G,WX(1,3))              
          CALL SUMW(WX,3,W(1,I,I+2))
       ENDDO
      
***************************************************************************

      DO i=1,ngluons-4
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+3),GG,0d0,0d0,WX(1,1))                
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+3),GG,0d0,0d0,WX(1,2))                
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+3),GG,0d0,0d0,WX(1,3))                
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+3),G,WX(1,4))                
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+3),G,WX(1,5))                
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+3),G,WX(1,6))                
         CALL SUMW(WX,6,W(1,I,I+3))
      ENDDO
      
***************************************************************************

      DO i=1,ngluons-5
C 3-point
         CALL VVV1_2(W(1,I,I+3),W(1,I+4,I+4),GG,0d0,0d0,WX(1, 1))               
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+4),GG,0d0,0d0,WX(1, 2))               
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+4),GG,0d0,0d0,WX(1, 3))               
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+4),GG,0d0,0d0,WX(1, 4))               
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
         CALL VVV1_2(W(1,I,I+4),W(1,I+5,I+5),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_2(W(1,I,I+3),W(1,I+4,I+5),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+5),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+5),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+5),GG,0d0,0d0,WX(1, 5))              
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
         CALL VVV1_2(W(1,I,I+5),W(1,I+6,I+6),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_2(W(1,I,I+4),W(1,I+5,I+6),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_2(W(1,I,I+3),W(1,I+4,I+6),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+6),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+6),GG,0d0,0d0,WX(1, 5))              
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+6),GG,0d0,0d0,WX(1, 6))              

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
         
         CALL SUMW(WX,21,jgluo8)
      ENDDO
      
      RETURN
      END

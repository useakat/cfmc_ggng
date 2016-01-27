      subroutine gluon7_7(w1,w2,w3,w4,w5,w6,g,jgluo7)

      implicit none
     
      integer i, ngluons
      double complex w(6,12,12),w1(6),w2(6),w3(6),w4(6),w5(6)
      double complex w6(6),w7(6),w8(6)
      double complex z(55),wx(6,55),jgluo7(6)
      double complex g,gg,ci
      parameter(ci=(0d0,1d0))
      
      vertex=(0d0,0d0)
      gg = -ci*g

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

      IF(ngluons.eq.8) then
         i=1

C 3-point
         CALL VVV1_0(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+7),GG,Z( 1))    
         CALL VVV1_0(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+7),GG,Z( 2))    
         CALL VVV1_0(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+7),GG,Z( 3))    
         CALL VVV1_0(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+7),GG,Z( 4))    
         CALL VVV1_0(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+7),GG,Z( 5))    
         CALL VVV1_0(W(1,I,I)  ,W(1,I+1,I+6),W(1,I+7,I+7),GG,Z( 6))    
         
C 4-point
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+5),  W(1,I+6,I+6),W(1,I+7,I+7), G,Z( 7))  
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+5),  W(1,I+6,I+6),W(1,I+7,I+7), G,Z( 8))  
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+5),  W(1,I+6,I+6),W(1,I+7,I+7), G,Z( 9))  
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+5),  W(1,I+6,I+6),W(1,I+7,I+7), G,Z(10))  
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+5),  W(1,I+6,I+6),W(1,I+7,I+7), G,Z(11))  
         
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+4),  W(1,I+5,I+6),W(1,I+7,I+7), G,Z(12))  
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+4),  W(1,I+5,I+6),W(1,I+7,I+7), G,Z(13))  
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+4),  W(1,I+5,I+6),W(1,I+7,I+7), G,Z(14))  
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+4),  W(1,I+5,I+6),W(1,I+7,I+7), G,Z(15))  
         
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+3),  W(1,I+4,I+6),W(1,I+7,I+7), G,Z(16))  
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+3),  W(1,I+4,I+6),W(1,I+7,I+7), G,Z(17))  
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+3),  W(1,I+4,I+6),W(1,I+7,I+7), G,Z(18))  
         
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+2),  W(1,I+3,I+6),W(1,I+7,I+7), G,Z(19))  
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+2),  W(1,I+3,I+6),W(1,I+7,I+7), G,Z(20))  
         
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+1),  W(1,I+2,I+6),W(1,I+7,I+7), G,Z(21))  
         DO K=1,21
            AMP_NGLUONS=AMP_NGLUONS+Z(K)
         ENDDO  
         RETURN
      ENDIF
 
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
         
         CALL SUMW(WX,21,W(1,I,I+6))
      ENDDO
      


***************************************************************************

      IF(ngluons.eq.9) then
         i=1

C 3-point
         CALL VVV1_0(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+8),GG,Z( 1))     
         CALL VVV1_0(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+8),GG,Z( 2))     
         CALL VVV1_0(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+8),GG,Z( 3))     
         CALL VVV1_0(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+8),GG,Z( 4))     
         CALL VVV1_0(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+8),GG,Z( 5))     
         CALL VVV1_0(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+8),GG,Z( 6))     
         CALL VVV1_0(W(1,I,I)  ,W(1,I+1,I+7),W(1,I+8,I+8),GG,Z( 7))     
         
C 4-point
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+7),W(1,I+8,I+8),G,Z( 8))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+7),W(1,I+8,I+8),G,Z( 9))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+7),W(1,I+8,I+8),G,Z(10))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+7),W(1,I+8,I+8),G,Z(11))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+7),W(1,I+8,I+8),G,Z(12))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+7),W(1,I+8,I+8),G,Z(13))     
         
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+7),W(1,I+8,I+8),G,Z(14))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+7),W(1,I+8,I+8),G,Z(15))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+7),W(1,I+8,I+8),G,Z(16))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+7),W(1,I+8,I+8),G,Z(17))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+5),W(1,I+6,I+7),W(1,I+8,I+8),G,Z(18))     
         
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+7),W(1,I+8,I+8),G,Z(19))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+7),W(1,I+8,I+8),G,Z(20))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+7),W(1,I+8,I+8),G,Z(21))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+7),W(1,I+8,I+8),G,Z(22))     
         
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+7),W(1,I+8,I+8),G,Z(23))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+7),W(1,I+8,I+8),G,Z(24))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+7),W(1,I+8,I+8),G,Z(25))     
         
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+7),W(1,I+8,I+8),G,Z(26))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+7),W(1,I+8,I+8),G,Z(27))     
         
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+7),W(1,I+8,I+8),G,Z(28))     
 
         DO K=1,28
            AMP_NGLUONS=AMP_NGLUONS+Z(K)
         ENDDO  
         RETURN
      ENDIF

 
      DO i=1,ngluons-8 

C 3-point
         CALL VVV1_2(W(1,I,I+6),W(1,I+7,I+7),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_2(W(1,I,I+5),W(1,I+6,I+7),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_2(W(1,I,I+4),W(1,I+5,I+7),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_2(W(1,I,I+3),W(1,I+4,I+7),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+7),GG,0d0,0d0,WX(1, 5))              
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+7),GG,0d0,0d0,WX(1, 6))              
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+7),GG,0d0,0d0,WX(1, 7))              

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

      IF(ngluons.eq.10) then
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
         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z( 9))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z(10))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z(11))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z(12))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z(13))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z(14))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+7),W(1,I+8,I+8),W(1,I+9,I+9),G,Z(15))     
         
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+8),W(1,I+9,I+9),G,Z(16))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+8),W(1,I+9,I+9),G,Z(17))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+8),W(1,I+9,I+9),G,Z(18))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+8),W(1,I+9,I+9),G,Z(19))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+8),W(1,I+9,I+9),G,Z(20))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+8),W(1,I+9,I+9),G,Z(21))     
         
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+8),W(1,I+9,I+9),G,Z(22))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+8),W(1,I+9,I+9),G,Z(23))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+8),W(1,I+9,I+9),G,Z(24))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+8),W(1,I+9,I+9),G,Z(25))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+5),W(1,I+6,I+8),W(1,I+9,I+9),G,Z(26))     
         
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+8),W(1,I+9,I+9),G,Z(27))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+8),W(1,I+9,I+9),G,Z(28))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+8),W(1,I+9,I+9),G,Z(29))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+8),W(1,I+9,I+9),G,Z(30))     
         
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+8),W(1,I+9,I+9),G,Z(31))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+8),W(1,I+9,I+9),G,Z(32))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+8),W(1,I+9,I+9),G,Z(33))     
         
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+8),W(1,I+9,I+9),G,Z(34))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+8),W(1,I+9,I+9),G,Z(35))     
         
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+8),W(1,I+9,I+9),G,Z(36))     
         
         DO K=1,36
            AMP_NGLUONS=AMP_NGLUONS+Z(K)
         ENDDO  
         RETURN

      ENDIF


      DO i=1,ngluons-9 

C 3-point
         CALL VVV1_2(W(1,I,I+7),W(1,I+8,I+8),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_2(W(1,I,I+6),W(1,I+7,I+8),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_2(W(1,I,I+5),W(1,I+6,I+8),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_2(W(1,I,I+4),W(1,I+5,I+8),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_2(W(1,I,I+3),W(1,I+4,I+8),GG,0d0,0d0,WX(1, 5))              
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+8),GG,0d0,0d0,WX(1, 6))              
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+8),GG,0d0,0d0,WX(1, 7))              
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+8),GG,0d0,0d0,WX(1, 8))              

C 4-point
         CALL JGGGCF(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+8), G,WX(1, 9))              
         CALL JGGGCF(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+8), G,WX(1,10))              
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+8), G,WX(1,11))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+8), G,WX(1,12))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+8), G,WX(1,13))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+8), G,WX(1,14))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+7),W(1,I+8,I+8), G,WX(1,15))              

         CALL JGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+8), G,WX(1,16))              
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+8), G,WX(1,17))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+8), G,WX(1,18))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+8), G,WX(1,19))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+8), G,WX(1,20))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+8), G,WX(1,21))              
         
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+8), G,WX(1,22))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+8), G,WX(1,23))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+8), G,WX(1,24))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+8), G,WX(1,25))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+5),W(1,I+6,I+8), G,WX(1,26))              
         
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+8), G,WX(1,27))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+8), G,WX(1,28))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+8), G,WX(1,29))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+8), G,WX(1,30))              
         
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+8), G,WX(1,31))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+8), G,WX(1,32))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+8), G,WX(1,33))              
         
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+8), G,WX(1,34))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+8), G,WX(1,35))              
         
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+8), G,WX(1,36))              
         
         CALL SUMW(WX,36,W(1,I,I+8))
      ENDDO

***************************************************************************

      IF(ngluons.eq.11) then
         i=1
      
C 3-point
         CALL VVV1_0(W(1,I,I+8),W(1,I+9,I+9),W(1,I+10,I+10),GG,Z( 1))     
         CALL VVV1_0(W(1,I,I+7),W(1,I+8,I+9),W(1,I+10,I+10),GG,Z( 2))     
         CALL VVV1_0(W(1,I,I+6),W(1,I+7,I+9),W(1,I+10,I+10),GG,Z( 3))     
         CALL VVV1_0(W(1,I,I+5),W(1,I+6,I+9),W(1,I+10,I+10),GG,Z( 4))     
         CALL VVV1_0(W(1,I,I+4),W(1,I+5,I+9),W(1,I+10,I+10),GG,Z( 5))     
         CALL VVV1_0(W(1,I,I+3),W(1,I+4,I+9),W(1,I+10,I+10),GG,Z( 6))     
         CALL VVV1_0(W(1,I,I+2),W(1,I+3,I+9),W(1,I+10,I+10),GG,Z( 7))     
         CALL VVV1_0(W(1,I,I+1),W(1,I+2,I+9),W(1,I+10,I+10),GG,Z( 8))     
         CALL VVV1_0(W(1,I,I)  ,W(1,I+1,I+9),W(1,I+10,I+10),GG,Z( 9))     

C 4-point

         CALL GGGGCF(W(1,I,I+7),W(1,I+8,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(10))     
         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(11))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(12))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(13))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(14))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(15))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(16))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+8),W(1,I+9,I+9),W(1,I+10,I+10),G,Z(17))     

         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(18))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(19))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(20))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(21))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(22))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(23))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+7),W(1,I+8,I+9),W(1,I+10,I+10),G,Z(24))     
         
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+9),W(1,I+10,I+10),G,Z(25))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+9),W(1,I+10,I+10),G,Z(26))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+9),W(1,I+10,I+10),G,Z(27))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+9),W(1,I+10,I+10),G,Z(28))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+9),W(1,I+10,I+10),G,Z(29))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+9),W(1,I+10,I+10),G,Z(30))     
         
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+9),W(1,I+10,I+10),G,Z(31))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+9),W(1,I+10,I+10),G,Z(32))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+9),W(1,I+10,I+10),G,Z(33))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+9),W(1,I+10,I+10),G,Z(34))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+5),W(1,I+6,I+9),W(1,I+10,I+10),G,Z(35))     
         
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+9),W(1,I+10,I+10),G,Z(36))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+9),W(1,I+10,I+10),G,Z(37))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+9),W(1,I+10,I+10),G,Z(38))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+9),W(1,I+10,I+10),G,Z(39))     
         
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+9),W(1,I+10,I+10),G,Z(40))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+9),W(1,I+10,I+10),G,Z(41))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+9),W(1,I+10,I+10),G,Z(42))     
         
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+9),W(1,I+10,I+10),G,Z(43))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+9),W(1,I+10,I+10),G,Z(44))     
         
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+9),W(1,I+10,I+10),G,Z(45))     
         
         DO K=1,45
            AMP_NGLUONS=AMP_NGLUONS+Z(K)
         ENDDO  

      RETURN
      ENDIF

***************************************************************************

      DO i=1,ngluons-10

C 3-point
         CALL VVV1_2(W(1,I,I+8),W(1,I+9,I+9),GG,0d0,0d0,WX(1, 1))              
         CALL VVV1_2(W(1,I,I+7),W(1,I+8,I+9),GG,0d0,0d0,WX(1, 2))              
         CALL VVV1_2(W(1,I,I+6),W(1,I+7,I+9),GG,0d0,0d0,WX(1, 3))              
         CALL VVV1_2(W(1,I,I+5),W(1,I+6,I+9),GG,0d0,0d0,WX(1, 4))              
         CALL VVV1_2(W(1,I,I+4),W(1,I+5,I+9),GG,0d0,0d0,WX(1, 5))              
         CALL VVV1_2(W(1,I,I+3),W(1,I+4,I+9),GG,0d0,0d0,WX(1, 6))              
         CALL VVV1_2(W(1,I,I+2),W(1,I+3,I+9),GG,0d0,0d0,WX(1, 7))              
         CALL VVV1_2(W(1,I,I+1),W(1,I+2,I+9),GG,0d0,0d0,WX(1, 8))              
         CALL VVV1_2(W(1,I,I)  ,W(1,I+1,I+9),GG,0d0,0d0,WX(1, 9))              

C 4-point
         CALL JGGGCF(W(1,I,I+7),W(1,I+8,I+8),W(1,I+9,I+9), G,WX(1,10))              
         CALL JGGGCF(W(1,I,I+6),W(1,I+7,I+8),W(1,I+9,I+9), G,WX(1,11))              
         CALL JGGGCF(W(1,I,I+5),W(1,I+6,I+8),W(1,I+9,I+9), G,WX(1,12))              
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+8),W(1,I+9,I+9), G,WX(1,13))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+8),W(1,I+9,I+9), G,WX(1,14))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+8),W(1,I+9,I+9), G,WX(1,15))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+8),W(1,I+9,I+9), G,WX(1,16))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+8),W(1,I+9,I+9), G,WX(1,17))              

         CALL JGGGCF(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+9), G,WX(1,18))              
         CALL JGGGCF(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+9), G,WX(1,19))              
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+9), G,WX(1,20))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+9), G,WX(1,21))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+9), G,WX(1,22))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+9), G,WX(1,23))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+7),W(1,I+8,I+9), G,WX(1,24))              
         
         CALL JGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+9), G,WX(1,25))              
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+9), G,WX(1,26))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+9), G,WX(1,27))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+9), G,WX(1,28))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+9), G,WX(1,29))              
         CALL JGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+9), G,WX(1,30))              
         
         CALL JGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+9), G,WX(1,31))              
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+9), G,WX(1,32))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+9), G,WX(1,33))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+9), G,WX(1,34))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+5),W(1,I+6,I+9), G,WX(1,35))              
         
         CALL JGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+9), G,WX(1,36))              
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+9), G,WX(1,37))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+9), G,WX(1,38))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+9), G,WX(1,39))              
         
         CALL JGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+9), G,WX(1,40))              
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+9), G,WX(1,41))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+9), G,WX(1,42))              
         
         CALL JGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+9), G,WX(1,43))              
         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+9), G,WX(1,44))              

         CALL JGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+9), G,WX(1,45))              

         CALL SUMW(WX,45,W(1,I,I+9))
      ENDDO

***************************************************************************

      IF(ngluons.eq.12) then
         i=1
      
C 3-point
         CALL VVV1_0(W(1,I,I+9),W(1,I+10,I+10),W(1,I+11,I+11),GG,Z( 1))     
         CALL VVV1_0(W(1,I,I+8),W(1,I+9 ,I+10),W(1,I+11,I+11),GG,Z( 2))     
         CALL VVV1_0(W(1,I,I+7),W(1,I+8 ,I+10),W(1,I+11,I+11),GG,Z( 3))     
         CALL VVV1_0(W(1,I,I+6),W(1,I+7 ,I+10),W(1,I+11,I+11),GG,Z( 4))     
         CALL VVV1_0(W(1,I,I+5),W(1,I+6 ,I+10),W(1,I+11,I+11),GG,Z( 5))     
         CALL VVV1_0(W(1,I,I+4),W(1,I+5 ,I+10),W(1,I+11,I+11),GG,Z( 6))     
         CALL VVV1_0(W(1,I,I+3),W(1,I+4 ,I+10),W(1,I+11,I+11),GG,Z( 7))     
         CALL VVV1_0(W(1,I,I+2),W(1,I+3 ,I+10),W(1,I+11,I+11),GG,Z( 8))     
         CALL VVV1_0(W(1,I,I+1),W(1,I+2 ,I+10),W(1,I+11,I+11),GG,Z( 9))     
         CALL VVV1_0(W(1,I,I)  ,W(1,I+1 ,I+10),W(1,I+11,I+11),GG,Z(10))     

C 4-point

         CALL GGGGCF(W(1,I,I+8),W(1,I+9,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(11))     
         CALL GGGGCF(W(1,I,I+7),W(1,I+8,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(12))     
         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(13))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(14))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(15))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(16))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(17))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(18))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+9),W(1,I+10,I+10),W(1,I+11,I+11),G,Z(19))     

         CALL GGGGCF(W(1,I,I+7),W(1,I+8,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(20))     
         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(21))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(22))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(23))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(24))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(25))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(26))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+8),W(1,I+9,I+10),W(1,I+11,I+11),G,Z(27))     
         
         CALL GGGGCF(W(1,I,I+6),W(1,I+7,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(28))     
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(29))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(30))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(31))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(32))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(33))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+7),W(1,I+8,I+10),W(1,I+11,I+11),G,Z(34))     
         
         CALL GGGGCF(W(1,I,I+5),W(1,I+6,I+6),W(1,I+7,I+10),W(1,I+11,I+11),G,Z(35))     
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+6),W(1,I+7,I+10),W(1,I+11,I+11),G,Z(36))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+6),W(1,I+7,I+10),W(1,I+11,I+11),G,Z(37))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+6),W(1,I+7,I+10),W(1,I+11,I+11),G,Z(38))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+6),W(1,I+7,I+10),W(1,I+11,I+11),G,Z(39))     
         CALL GGGGCF(W(1,I,I  ),W(1,I+1,I+6),W(1,I+7,I+10),W(1,I+11,I+11),G,Z(40))     
         
         CALL GGGGCF(W(1,I,I+4),W(1,I+5,I+5),W(1,I+6,I+10),W(1,I+11,I+11),G,Z(41))     
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+5),W(1,I+6,I+10),W(1,I+11,I+11),G,Z(42))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+5),W(1,I+6,I+10),W(1,I+11,I+11),G,Z(43))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+5),W(1,I+6,I+10),W(1,I+11,I+11),G,Z(44))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+5),W(1,I+6,I+10),W(1,I+11,I+11),G,Z(45))     
         
         CALL GGGGCF(W(1,I,I+3),W(1,I+4,I+4),W(1,I+5,I+10),W(1,I+11,I+11),G,Z(46))     
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+4),W(1,I+5,I+10),W(1,I+11,I+11),G,Z(47))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+4),W(1,I+5,I+10),W(1,I+11,I+11),G,Z(48))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+4),W(1,I+5,I+10),W(1,I+11,I+11),G,Z(49))     
         
         CALL GGGGCF(W(1,I,I+2),W(1,I+3,I+3),W(1,I+4,I+10),W(1,I+11,I+11),G,Z(50))     
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+3),W(1,I+4,I+10),W(1,I+11,I+11),G,Z(51))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+3),W(1,I+4,I+10),W(1,I+11,I+11),G,Z(52))     
         
         CALL GGGGCF(W(1,I,I+1),W(1,I+2,I+2),W(1,I+3,I+10),W(1,I+11,I+11),G,Z(53))     
         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+2),W(1,I+3,I+10),W(1,I+11,I+11),G,Z(54))     

         CALL GGGGCF(W(1,I,I)  ,W(1,I+1,I+1),W(1,I+2,I+10),W(1,I+11,I+11),G,Z(55))     
         
         DO K=1,55
            AMP_NGLUONS=AMP_NGLUONS+Z(K)
         ENDDO  

      RETURN
      ENDIF
      
      
      END

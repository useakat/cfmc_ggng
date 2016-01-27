C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV3_1(V2, V3, V4, COUP, M1, W1, V1)
      IMPLICIT NONE
      DOUBLE COMPLEX V1(*)
      DOUBLE COMPLEX V2(*)
      DOUBLE COMPLEX V3(*)
      DOUBLE COMPLEX V4(*)
      DOUBLE COMPLEX COUP
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M1, W1
      DOUBLE COMPLEX OM1
      DOUBLE PRECISION P1(0:3)

      V1(5)= V2(5)+V3(5)+V4(5)
      V1(6)= V2(6)+V3(6)+V4(6)
      P1(0) = - DBLE(V1(5))
      P1(1) = - DBLE(V1(6))
      P1(2) = - DIMAG(V1(6))
      P1(3) = - DIMAG(V1(5))
      OM1 = 0D0
      IF (M1 .NE. 0D0) OM1=1D0/M1**2

      DENOM =1D0/(( (M1*( -M1+(0, 1)*W1))+( (P1(0)**2)-(P1(1)**2)
     $ -(P1(2)**2)-(P1(3)**2))))
      V1(1)= COUP*DENOM*( (OM1*( (V4(1)*( (V3(1)*( (0, 1)*(V2(2)
     $ *P1(1))+(0, 1)*(V2(3)*P1(2))+(0, 1)*(V2(4)*P1(3))))+(P1(0)
     $ *( (0, -1)*(V3(2)*V2(2))+(0, -1)*(V3(3)*V2(3))+(0, -1)*(V3(4)
     $ *V2(4))))))+( (V4(2)*( (V3(2)*( (0, 1)*(V2(1)*P1(0))+(0, 
     $ -1)*(V2(3)*P1(2))+(0, -1)*(V2(4)*P1(3))))+(P1(1)*( (0, 
     $ -1)*(V3(1)*V2(1))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)*V2(4))))))
     $ +( (V4(3)*( (V3(3)*( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))
     $ +(0, -1)*(V2(4)*P1(3))))+(P1(2)*( (0, -1)*(V3(1)*V2(1))
     $ +(0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(4)*V2(4))))))+(V4(4)*( (V3(4)
     $ *( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))+(0, -1)*(V2(3)
     $ *P1(2))))+(P1(3)*( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(2)*V2(2))
     $ +(0, 1)*(V3(3)*V2(3)))))))))*P1(0))+( (V2(1)*( (0, -1)*(V4(2)
     $ *V3(2))+(0, -1)*(V4(3)*V3(3))+(0, -1)*(V4(4)*V3(4))))+(V4(1)
     $ *( (0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)
     $ *V2(4))))))
      V1(2)= COUP*DENOM*( (OM1*( (V4(1)*( (V3(1)*( (0, 1)*(V2(2)
     $ *P1(1))+(0, 1)*(V2(3)*P1(2))+(0, 1)*(V2(4)*P1(3))))+(P1(0)
     $ *( (0, -1)*(V3(2)*V2(2))+(0, -1)*(V3(3)*V2(3))+(0, -1)*(V3(4)
     $ *V2(4))))))+( (V4(2)*( (V3(2)*( (0, 1)*(V2(1)*P1(0))+(0, 
     $ -1)*(V2(3)*P1(2))+(0, -1)*(V2(4)*P1(3))))+(P1(1)*( (0, 
     $ -1)*(V3(1)*V2(1))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)*V2(4))))))
     $ +( (V4(3)*( (V3(3)*( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))
     $ +(0, -1)*(V2(4)*P1(3))))+(P1(2)*( (0, -1)*(V3(1)*V2(1))
     $ +(0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(4)*V2(4))))))+(V4(4)*( (V3(4)
     $ *( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))+(0, -1)*(V2(3)
     $ *P1(2))))+(P1(3)*( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(2)*V2(2))
     $ +(0, 1)*(V3(3)*V2(3)))))))))*P1(1))+( (V2(2)*( (0, 1)*(V4(1)
     $ *V3(1))+(0, -1)*(V4(3)*V3(3))+(0, -1)*(V4(4)*V3(4))))+(V4(2)
     $ *( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)
     $ *V2(4))))))
      V1(3)= COUP*DENOM*( (OM1*( (V4(1)*( (V3(1)*( (0, 1)*(V2(2)
     $ *P1(1))+(0, 1)*(V2(3)*P1(2))+(0, 1)*(V2(4)*P1(3))))+(P1(0)
     $ *( (0, -1)*(V3(2)*V2(2))+(0, -1)*(V3(3)*V2(3))+(0, -1)*(V3(4)
     $ *V2(4))))))+( (V4(2)*( (V3(2)*( (0, 1)*(V2(1)*P1(0))+(0, 
     $ -1)*(V2(3)*P1(2))+(0, -1)*(V2(4)*P1(3))))+(P1(1)*( (0, 
     $ -1)*(V3(1)*V2(1))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)*V2(4))))))
     $ +( (V4(3)*( (V3(3)*( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))
     $ +(0, -1)*(V2(4)*P1(3))))+(P1(2)*( (0, -1)*(V3(1)*V2(1))
     $ +(0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(4)*V2(4))))))+(V4(4)*( (V3(4)
     $ *( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))+(0, -1)*(V2(3)
     $ *P1(2))))+(P1(3)*( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(2)*V2(2))
     $ +(0, 1)*(V3(3)*V2(3)))))))))*P1(2))+( (V2(3)*( (0, 1)*(V4(1)
     $ *V3(1))+(0, -1)*(V4(2)*V3(2))+(0, -1)*(V4(4)*V3(4))))+(V4(3)
     $ *( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(4)
     $ *V2(4))))))
      V1(4)= COUP*DENOM*( (OM1*( (V4(1)*( (V3(1)*( (0, 1)*(V2(2)
     $ *P1(1))+(0, 1)*(V2(3)*P1(2))+(0, 1)*(V2(4)*P1(3))))+(P1(0)
     $ *( (0, -1)*(V3(2)*V2(2))+(0, -1)*(V3(3)*V2(3))+(0, -1)*(V3(4)
     $ *V2(4))))))+( (V4(2)*( (V3(2)*( (0, 1)*(V2(1)*P1(0))+(0, 
     $ -1)*(V2(3)*P1(2))+(0, -1)*(V2(4)*P1(3))))+(P1(1)*( (0, 
     $ -1)*(V3(1)*V2(1))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)*V2(4))))))
     $ +( (V4(3)*( (V3(3)*( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))
     $ +(0, -1)*(V2(4)*P1(3))))+(P1(2)*( (0, -1)*(V3(1)*V2(1))
     $ +(0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(4)*V2(4))))))+(V4(4)*( (V3(4)
     $ *( (0, 1)*(V2(1)*P1(0))+(0, -1)*(V2(2)*P1(1))+(0, -1)*(V2(3)
     $ *P1(2))))+(P1(3)*( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(2)*V2(2))
     $ +(0, 1)*(V3(3)*V2(3)))))))))*P1(3))+( (V2(4)*( (0, 1)*(V4(1)
     $ *V3(1))+(0, -1)*(V4(2)*V3(2))+(0, -1)*(V4(3)*V3(3))))+(V4(4)
     $ *( (0, -1)*(V3(1)*V2(1))+(0, 1)*(V3(2)*V2(2))+(0, 1)*(V3(3)
     $ *V2(3))))))
      END



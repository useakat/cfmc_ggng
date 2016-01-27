C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) +
C      P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)
C     
      SUBROUTINE VVV1_3(V1, V2, C, M3, W3, V3)
      IMPLICIT NONE
      DOUBLE COMPLEX V1(6)
      DOUBLE COMPLEX V2(6)
      DOUBLE COMPLEX V3(6)
      DOUBLE COMPLEX C
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M3, W3
      DOUBLE COMPLEX OM3
      DOUBLE PRECISION P2(0:3),P3(0:3),P1(0:3)

      V3(5)= V1(5)+V2(5)
      V3(6)= V1(6)+V2(6)
      P2(0) =  DBLE(V2(5))
      P2(1) =  DBLE(V2(6))
      P2(2) =  DIMAG(V2(6))
      P2(3) =  DIMAG(V2(5))
      P3(0) = - DBLE(V3(5))
      P3(1) = - DBLE(V3(6))
      P3(2) = - DIMAG(V3(6))
      P3(3) = - DIMAG(V3(5))
      P1(0) =  DBLE(V1(5))
      P1(1) =  DBLE(V1(6))
      P1(2) =  DIMAG(V1(6))
      P1(3) =  DIMAG(V1(5))
      OM3 = 0D0
      IF (M3 .NE. 0D0) OM3=1D0/DCMPLX(M3**2,-W3*M3)

      DENOM =1D0/(( (M3*( -M3+(0, 1)*W3))+( (P3(0)**2)-(P3(1)**2)
     $ -(P3(2)**2)-(P3(3)**2))))
      V3(1)= C*DENOM*( (OM3*( (V2(1)*( (V1(1)*( (P3(1)*( (0, -1)*P1(1)
     $ +(0, 1)*P2(1)))+( (P3(2)*( (0, -1)*P1(2)+(0, 1)*P2(2)))
     $ +(P3(3)*( (0, -1)*P1(3)+(0, 1)*P2(3))))))+( (P1(0)*( (0, 1)
     $ *(V1(2)*P3(1))+(0, 1)*(V1(3)*P3(2))+(0, 1)*(V1(4)*P3(3))))
     $ +(P3(0)*( (0, -1)*(V1(2)*P2(1))+(0, -1)*(V1(3)*P2(2))+(0, 
     $ -1)*(V1(4)*P2(3)))))))+( (V2(2)*( (V1(2)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(2)*( (0, 1)*P1(2)+(0, -1)
     $ *P2(2)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(1)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(3)*P3(2))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(1)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(3)*P2(2))
     $ +(0, 1)*(V1(4)*P2(3)))))))+( (V2(3)*( (V1(3)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(2)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(2)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(4)*P2(3)))))))+(V2(4)*( (V1(4)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(2)*( (0, 1)*P1(2)+(0, -1)*P2(2))))))+( (P1(3)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(3)
     $ *P3(2))))+(P3(3)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(3)*P2(2))))))))))*P3(0))+( (V1(1)*( (V2(2)*( (0, 
     $ -1)*P1(1)+(0, 1)*P3(1)))+( (V2(3)*( (0, -1)*P1(2)+(0, 1)
     $ *P3(2)))+(V2(4)*( (0, -1)*P1(3)+(0, 1)*P3(3))))))+( (V2(1)
     $ *( (V1(2)*( (0, 1)*P2(1)+(0, -1)*P3(1)))+( (V1(3)*( (0, 1)
     $ *P2(2)+(0, -1)*P3(2)))+(V1(4)*( (0, 1)*P2(3)+(0, -1)*P3(3))))))
     $ +( (P1(0)*( (0, 1)*(V2(2)*V1(2))+(0, 1)*(V2(3)*V1(3))+(0, 1)
     $ *(V2(4)*V1(4))))+(P2(0)*( (0, -1)*(V2(2)*V1(2))+(0, -1)*(V2(3)
     $ *V1(3))+(0, -1)*(V2(4)*V1(4))))))))
      V3(2)= C*DENOM*( (OM3*( (V2(1)*( (V1(1)*( (P3(1)*( (0, -1)*P1(1)
     $ +(0, 1)*P2(1)))+( (P3(2)*( (0, -1)*P1(2)+(0, 1)*P2(2)))
     $ +(P3(3)*( (0, -1)*P1(3)+(0, 1)*P2(3))))))+( (P1(0)*( (0, 1)
     $ *(V1(2)*P3(1))+(0, 1)*(V1(3)*P3(2))+(0, 1)*(V1(4)*P3(3))))
     $ +(P3(0)*( (0, -1)*(V1(2)*P2(1))+(0, -1)*(V1(3)*P2(2))+(0, 
     $ -1)*(V1(4)*P2(3)))))))+( (V2(2)*( (V1(2)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(2)*( (0, 1)*P1(2)+(0, -1)
     $ *P2(2)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(1)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(3)*P3(2))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(1)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(3)*P2(2))
     $ +(0, 1)*(V1(4)*P2(3)))))))+( (V2(3)*( (V1(3)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(2)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(2)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(4)*P2(3)))))))+(V2(4)*( (V1(4)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(2)*( (0, 1)*P1(2)+(0, -1)*P2(2))))))+( (P1(3)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(3)
     $ *P3(2))))+(P3(3)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(3)*P2(2))))))))))*P3(1))+( (V1(2)*( (V2(1)*( (0, 1)
     $ *P1(0)+(0, -1)*P3(0)))+( (V2(3)*( (0, -1)*P1(2)+(0, 1)*P3(2)))
     $ +(V2(4)*( (0, -1)*P1(3)+(0, 1)*P3(3))))))+( (V2(2)*( (V1(1)
     $ *( (0, -1)*P2(0)+(0, 1)*P3(0)))+( (V1(3)*( (0, 1)*P2(2)
     $ +(0, -1)*P3(2)))+(V1(4)*( (0, 1)*P2(3)+(0, -1)*P3(3))))))
     $ +( (P1(1)*( (0, -1)*(V2(1)*V1(1))+(0, 1)*(V2(3)*V1(3))
     $ +(0, 1)*(V2(4)*V1(4))))+(P2(1)*( (0, 1)*(V2(1)*V1(1))+(0, 
     $ -1)*(V2(3)*V1(3))+(0, -1)*(V2(4)*V1(4))))))))
      V3(3)= C*DENOM*( (OM3*( (V2(1)*( (V1(1)*( (P3(1)*( (0, -1)*P1(1)
     $ +(0, 1)*P2(1)))+( (P3(2)*( (0, -1)*P1(2)+(0, 1)*P2(2)))
     $ +(P3(3)*( (0, -1)*P1(3)+(0, 1)*P2(3))))))+( (P1(0)*( (0, 1)
     $ *(V1(2)*P3(1))+(0, 1)*(V1(3)*P3(2))+(0, 1)*(V1(4)*P3(3))))
     $ +(P3(0)*( (0, -1)*(V1(2)*P2(1))+(0, -1)*(V1(3)*P2(2))+(0, 
     $ -1)*(V1(4)*P2(3)))))))+( (V2(2)*( (V1(2)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(2)*( (0, 1)*P1(2)+(0, -1)
     $ *P2(2)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(1)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(3)*P3(2))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(1)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(3)*P2(2))
     $ +(0, 1)*(V1(4)*P2(3)))))))+( (V2(3)*( (V1(3)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(2)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(2)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(4)*P2(3)))))))+(V2(4)*( (V1(4)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(2)*( (0, 1)*P1(2)+(0, -1)*P2(2))))))+( (P1(3)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(3)
     $ *P3(2))))+(P3(3)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(3)*P2(2))))))))))*P3(2))+( (V1(3)*( (V2(1)*( (0, 1)
     $ *P1(0)+(0, -1)*P3(0)))+( (V2(2)*( (0, -1)*P1(1)+(0, 1)*P3(1)))
     $ +(V2(4)*( (0, -1)*P1(3)+(0, 1)*P3(3))))))+( (V2(3)*( (V1(1)
     $ *( (0, -1)*P2(0)+(0, 1)*P3(0)))+( (V1(2)*( (0, 1)*P2(1)
     $ +(0, -1)*P3(1)))+(V1(4)*( (0, 1)*P2(3)+(0, -1)*P3(3))))))
     $ +( (P1(2)*( (0, -1)*(V2(1)*V1(1))+(0, 1)*(V2(2)*V1(2))
     $ +(0, 1)*(V2(4)*V1(4))))+(P2(2)*( (0, 1)*(V2(1)*V1(1))+(0, 
     $ -1)*(V2(2)*V1(2))+(0, -1)*(V2(4)*V1(4))))))))
      V3(4)= C*DENOM*( (OM3*( (V2(1)*( (V1(1)*( (P3(1)*( (0, -1)*P1(1)
     $ +(0, 1)*P2(1)))+( (P3(2)*( (0, -1)*P1(2)+(0, 1)*P2(2)))
     $ +(P3(3)*( (0, -1)*P1(3)+(0, 1)*P2(3))))))+( (P1(0)*( (0, 1)
     $ *(V1(2)*P3(1))+(0, 1)*(V1(3)*P3(2))+(0, 1)*(V1(4)*P3(3))))
     $ +(P3(0)*( (0, -1)*(V1(2)*P2(1))+(0, -1)*(V1(3)*P2(2))+(0, 
     $ -1)*(V1(4)*P2(3)))))))+( (V2(2)*( (V1(2)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(2)*( (0, 1)*P1(2)+(0, -1)
     $ *P2(2)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(1)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(3)*P3(2))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(1)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(3)*P2(2))
     $ +(0, 1)*(V1(4)*P2(3)))))))+( (V2(3)*( (V1(3)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(3)*( (0, 1)*P1(3)+(0, -1)*P2(3))))))+( (P1(2)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(4)
     $ *P3(3))))+(P3(2)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(4)*P2(3)))))))+(V2(4)*( (V1(4)*( (P3(0)*( (0, 
     $ -1)*P1(0)+(0, 1)*P2(0)))+( (P3(1)*( (0, 1)*P1(1)+(0, -1)
     $ *P2(1)))+(P3(2)*( (0, 1)*P1(2)+(0, -1)*P2(2))))))+( (P1(3)
     $ *( (0, 1)*(V1(1)*P3(0))+(0, -1)*(V1(2)*P3(1))+(0, -1)*(V1(3)
     $ *P3(2))))+(P3(3)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(3)*P2(2))))))))))*P3(3))+( (V1(4)*( (V2(1)*( (0, 1)
     $ *P1(0)+(0, -1)*P3(0)))+( (V2(2)*( (0, -1)*P1(1)+(0, 1)*P3(1)))
     $ +(V2(3)*( (0, -1)*P1(2)+(0, 1)*P3(2))))))+( (V2(4)*( (V1(1)
     $ *( (0, -1)*P2(0)+(0, 1)*P3(0)))+( (V1(2)*( (0, 1)*P2(1)
     $ +(0, -1)*P3(1)))+(V1(3)*( (0, 1)*P2(2)+(0, -1)*P3(2))))))
     $ +( (P1(3)*( (0, -1)*(V2(1)*V1(1))+(0, 1)*(V2(2)*V1(2))
     $ +(0, 1)*(V2(3)*V1(3))))+(P2(3)*( (0, 1)*(V2(1)*V1(1))+(0, 
     $ -1)*(V2(2)*V1(2))+(0, -1)*(V2(3)*V1(3))))))))

      END
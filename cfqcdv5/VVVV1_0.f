C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)
C     
      SUBROUTINE VVVV1_0(V1,V2,V3,V4,COUP,VERTEX)
      IMPLICIT NONE
      DOUBLE COMPLEX V1(*)
      DOUBLE COMPLEX V2(*)
      DOUBLE COMPLEX V3(*)
      DOUBLE COMPLEX V4(*)
      DOUBLE COMPLEX COUP
      DOUBLE COMPLEX VERTEX


      VERTEX = COUP*( (V4(1)*( (V2(1)*( (0, -1)*(V3(2)*V1(2))
     $ +(0, -1)*(V3(3)*V1(3))+(0, -1)*(V3(4)*V1(4))))+(V1(1)*( (0, 1)
     $ *(V3(2)*V2(2))+(0, 1)*(V3(3)*V2(3))+(0, 1)*(V3(4)*V2(4))))))
     $ +( (V4(2)*( (V2(2)*( (0, -1)*(V3(1)*V1(1))+(0, 1)*(V3(3)*V1(3))
     $ +(0, 1)*(V3(4)*V1(4))))+(V1(2)*( (0, 1)*(V3(1)*V2(1))+(0, 
     $ -1)*(V3(3)*V2(3))+(0, -1)*(V3(4)*V2(4))))))+( (V4(3)*( (V2(3)
     $ *( (0, -1)*(V3(1)*V1(1))+(0, 1)*(V3(2)*V1(2))+(0, 1)*(V3(4)
     $ *V1(4))))+(V1(3)*( (0, 1)*(V3(1)*V2(1))+(0, -1)*(V3(2)*V2(2))
     $ +(0, -1)*(V3(4)*V2(4))))))+(V4(4)*( (V2(4)*( (0, -1)*(V3(1)
     $ *V1(1))+(0, 1)*(V3(2)*V1(2))+(0, 1)*(V3(3)*V1(3))))+(V1(4)
     $ *( (0, 1)*(V3(1)*V2(1))+(0, -1)*(V3(2)*V2(2))+(0, -1)*(V3(3)
     $ *V2(3)))))))))
      END


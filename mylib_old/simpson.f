        SUBROUTINE SIMP(F,X0,X1,ANS,ERROR,NMAX,NFLAG)                   
C--------THIS SUBROUTIN INTEGRATES A FUNCTION F FROM X0 TO X1 AND RETURN
C--------THE INTEGRAL AS ANS BY USING THE SIMPSON'S FORMULA.  ALLOWED   
C--------ERROR AND A MAXIMUM NUMBER OF ITERATION NMAX SHOULD BE         
C--------SPECIFIED. (THE MAXIMUM DIVISION NUMBER IS 2*2**NMAX).         
C--------NFLAG=1 IF THE INTEGRAL DOES NOT CONVERGE.                     
C--------NFLAG=0 IF THE INTEGRAL CONVERGES.                             
        IMPLICIT REAL*8 (A-H,O-Z)
        NFLAG=0
        N   =2
        H   =(X1-X0)/2.
       
        SUM0=F(X1)+F(X0)
       
30        SUM1=F(X0+H)
        
        SUM2=0.
        S   =H/3.*(SUM0+4.*SUM1)
        
       DO 10 L=1,NMAX
        H   =H/2.
        N   =N*2
        SUM2=SUM1+SUM2
        SUM1=0.
        DO 20 I=1,N,2
20      SUM1=SUM1+F(X0+H*I)
        ANS =H/3.*(SUM0+4.*SUM1+2.*SUM2)
        IF(ABS(ANS-S).LE.ABS(ANS*ERROR)) RETURN
        S=ANS
10      CONTINUE
        NFLAG=1
        RETURN
        END

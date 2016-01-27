c====================================================
      SUBROUTINE DK2(IY,I1,I2,ZZ,WPS2,Jacob2,IFLG2) 
cc  I2 carries fz I1 carries 1-fz
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
c      include 'cparam.inc'
      INTEGER  IFLG2
      DIMENSION zz(2)
      REAL*8 Jacob2
      COMMON /PARTCL/P(5,50),M(3,50)
      REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
      COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
      REAL*8 PTcut2,AYEcut2,AYEcut3
      COMMON /cut2/PTcut2,AYEcut2,AYEcut3
      REAL*8 lx,lol,uuz,fz,lamd

      WPS2=0.D0 
      IFLG2=0.d0 
      MY=M(1,IY)
      MYSQ=M(2,IY)
      M1SQ=M(2,I1)
      M2SQ=M(2,I2)
      DIFF = MY-M(1,I1)-M(1,I2)
      XLA=((MYSQ-M2SQ-M1SQ)**2-4.*M1SQ*M2SQ)
      lamd=dsqrt(XLA)/MYSQ
      IF(DIFF.LE.0.) GO TO 105
      IF(XLA.GT.0)GO TO 110
105   IFLG2=1 
      RETURN
110   XLA=SQRT(XLA)

      EXCM=(MYSQ+M2SQ-M1SQ)/(2.*MY) 
      PXCM=XLA/(2.*MY)
c.......!!!!!!!!!!!!!!!!!!
      lx=PTcut/P(4,IY)
c      lx=PTcut/14000d0
      lol=dlog(lx/(1.0d0-lx))
      uuz=lol+2*abs(lol)*ZZ(1)
      fz=dexp(uuz)/(1+dexp(uuz))
      CV=(2*fz-1.0d0-M2SQ/MYSQ+M1SQ/MYSQ)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/lamd   
      IF(abs(CV).GT.1.) then 
      IFLG2=1
      RETURN
      ENDIF  

      WPS2=6.0d0*Pi*abs(lol)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/(2.d0*pi)**(3*2.d0-4.0d0)
      Jacob2=1d0/6d0/(1d0/fz+1d0/(1d0-fz))

      P(1,I1)=0.
      P(2,I1)=0.
      P(3,I1)=-PXCM 
      P(4,I1)=(MYSQ+M1SQ-M2SQ)/(2.*MY)
      P(1,I2)=0.
      P(2,I2)=0.
      P(3,I2)=PXCM
      P(4,I2)=(MYSQ+M2SQ-M1SQ)/(2.*MY)
      SV=SQRT(ABS(1.-CV*CV))
      CALL ROT12(I1,1,3,SV,CV)
      CALL ROT12(I2,1,3,SV,CV)
      PH = 2.0*PI*ZZ(2) 
      SPH=SIN(PH)
       CPH=COS(PH)
      CALL ROT12(I1,2,1,SPH,CPH)
      CALL ROT12(I2,2,1,SPH,CPH)
C
C     BOOST/ROTATE TO LAB FRAME
C
      ETAY=P(5,IY)/MY
      IF(ETAY.LE.1.E-4)GO TO 200
      GAMMAY=P(4,IY)/MY
      CALL BOOSTZ(I1,GAMMAY,ETAY)
      CALL BOOSTZ(I2,GAMMAY,ETAY)
      CV=P(3,IY)/P(5,IY)
      SV=1.-CV*CV
      IF(SV.GT.-.001)GO TO 130
      WPS2=0. 
      IFLG2=2 
      RETURN
130   SV=SQRT(ABS(SV))
      CALL ROT12(I1,1,3,SV,CV)
      CALL ROT12(I2,1,3,SV,CV)
      PTY=SQRT(P(1,IY)**2+P(2,IY)**2)
      IF(PTY.LE.1.E-4)GO TO 200
      CPH=P(1,IY)/PTY
        SPH=P(2,IY)/PTY
      CALL ROT12(I1,2,1,SPH,CPH)
      CALL ROT12(I2,2,1,SPH,CPH)
200   CALL PSET(I1)
      CALL PSET(I2)
      RETURN
      END


c$$$      SUBROUTINE DK2_2(IY,I1,I2,ZZ,WPS2,Jacob2,IFLG2) 
c$$$cc  I2 carries fz I1 carries 1-fz
c$$$      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
c$$$      INTEGER  IFLG2
c$$$      DIMENSION zz(2)
c$$$      REAL*8 Jacob2
c$$$      COMMON /PARTCL/P(5,50),M(3,50)
c$$$      REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
c$$$      COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
c$$$      REAL*8 PTcut2,AYEcut2,AYEcut3
c$$$      COMMON /cut2/PTcut2,AYEcut2,AYEcut3
c$$$      REAL*8 lx,lol,uuz,fz,lamd
c$$$
c$$$      WPS2=0.D0 
c$$$      IFLG2=0.d0 
c$$$      MY=M(1,IY)
c$$$      MYSQ=M(2,IY)
c$$$      M1SQ=M(2,I1)
c$$$      M2SQ=M(2,I2)
c$$$      DIFF = MY-M(1,I1)-M(1,I2)
c$$$      XLA=((MYSQ-M2SQ-M1SQ)**2-4.*M1SQ*M2SQ)
c$$$      lamd=dsqrt(XLA)/MYSQ
c$$$      IF(DIFF.LE.0.) GO TO 105
c$$$      IF(XLA.GT.0)GO TO 110
c$$$105   IFLG2=1 
c$$$      RETURN
c$$$110   XLA=SQRT(XLA)
c$$$
c$$$      EXCM=(MYSQ+M2SQ-M1SQ)/(2.*MY) 
c$$$      PXCM=XLA/(2.*MY)
c$$$c.......!!!!!!!!!!!!!!!!!!
c$$$      lx=PTcut/P(4,IY)
c$$$      lol=dlog(lx/(1.0d0-lx))
c$$$      uuz=lol+2*abs(lol)*ZZ(1)
c$$$      fz=dexp(uuz)/(1+dexp(uuz))
c$$$      CV=(2*fz-1.0d0-M2SQ/MYSQ+M1SQ/MYSQ)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/lamd   
c$$$      IF(abs(CV).GT.1.) then 
c$$$      IFLG2=1
c$$$      RETURN
c$$$      ENDIF  
c$$$
c$$$      WPS2=6.0d0*Pi*abs(lol)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/(2.d0*pi)**(3*2.d0-4.0d0)
c$$$      Jacob2=1d0/6d0/(1d0/fz+1d0/(1d0-fz))
c$$$
c$$$      P(1,I1)=0.
c$$$      P(2,I1)=0.
c$$$      P(3,I1)=-PXCM 
c$$$      P(4,I1)=(MYSQ+M1SQ-M2SQ)/(2.*MY)
c$$$      P(1,I2)=0.
c$$$      P(2,I2)=0.
c$$$      P(3,I2)=PXCM
c$$$      P(4,I2)=(MYSQ+M2SQ-M1SQ)/(2.*MY)
c$$$      SV=SQRT(ABS(1.-CV*CV))
c$$$      CALL ROT12(I1,1,3,SV,CV)
c$$$      CALL ROT12(I2,1,3,SV,CV)
c$$$      PH = 2.0*PI*ZZ(2) 
c$$$      SPH=SIN(PH)
c$$$       CPH=COS(PH)
c$$$      CALL ROT12(I1,2,1,SPH,CPH)
c$$$      CALL ROT12(I2,2,1,SPH,CPH)
c$$$C
c$$$C     BOOST/ROTATE TO LAB FRAME
c$$$C
c$$$      ETAY=P(5,IY)/MY
c$$$      IF(ETAY.LE.1.E-4)GO TO 200
c$$$      GAMMAY=P(4,IY)/MY
c$$$      CALL BOOSTZ(I1,GAMMAY,ETAY)
c$$$      CALL BOOSTZ(I2,GAMMAY,ETAY)
c$$$      CV=P(3,IY)/P(5,IY)
c$$$      SV=1.-CV*CV
c$$$      IF(SV.GT.-.001)GO TO 130
c$$$      WPS2=0. 
c$$$      IFLG2=2 
c$$$      RETURN
c$$$130   SV=SQRT(ABS(SV))
c$$$      CALL ROT12(I1,1,3,SV,CV)
c$$$      CALL ROT12(I2,1,3,SV,CV)
c$$$      PTY=SQRT(P(1,IY)**2+P(2,IY)**2)
c$$$      IF(PTY.LE.1.E-4)GO TO 200
c$$$      CPH=P(1,IY)/PTY
c$$$        SPH=P(2,IY)/PTY
c$$$      CALL ROT12(I1,2,1,SPH,CPH)
c$$$      CALL ROT12(I2,2,1,SPH,CPH)
c$$$200   CALL PSET(I1)
c$$$      CALL PSET(I2)
c$$$      RETURN
c$$$      END

c==============================================================
C-------------------------------THREE-BODY DECAY: AUXILIARY SUBROUTINES
C----------------------------------------------------------------------
C ................. BOOSTN -  general boost ...........................
C
      SUBROUTINE BOOSTN(P,R,Q)
C
C
C     The four vector P is assumed to be given in the rest frame of R,
C     which must be a timelike vector.
C     output Q is the vector P boosted to the frame in which R is given.
C                                              Compare Jackson, p.517
C                                              D. Zeppenfeld (28.6.1985)
C
      REAL*4 P(5),R(4),Q(5)
      REAL*4 BETA(3), X, Y, GAMMA
      INTEGER I

      X = 0D0
      Y = 0D0
      DO I = 1,3
         BETA(I) = R(I)/R(4)
         X = X + BETA(I)**2
         Y = Y + BETA(I)*P(I)
      ENDDO
      IF (X.GT.1D-16.AND.X.LT.(1D0-1D-12)) THEN
         GAMMA = 1D0/DSQRT(1D0-X)
         DO I = 1,3
            Q(I) = P(I)+BETA(I)*(Y*(GAMMA-1D0)/X + GAMMA*P(4))
         ENDDO
         Q(4) = GAMMA*(P(4) + Y)
      ELSE
         DO I = 1,4
            Q(I) = P(I)
         ENDDO
         IF(X.GE.(1D0-1D-12))
     *      WRITE(6,1000) R,R(4)**2-R(1)**2-R(2)**2-R(3)**2
      ENDIF

      Q(5) = SQRT(Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3))
 1000 FORMAT (' The reference vector ',4G12.3,' is not timelike.'/
     1        ' R**2 = ',G12.3)
      END
C............................................Lorentz boost along Z axis
       SUBROUTINE BOOSTZ(I,GAMMA,ETA)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
         DOUBLE PRECISION M
         COMMON /PARTCL/ P(5,50),M(3,50)

         TMP    = GAMMA*P(3,I) + ETA*P(4,I)
         P(4,I) = GAMMA*P(4,I) + ETA*P(3,I)
         P(3,I) = TMP
         P(5,I) = SQRT(P(1,I)**2 + P(2,I)**2 + P(3,I)**2)
       RETURN 
       END
C............................................Lorentz boost along X axis
       SUBROUTINE BOOST1(I,GAMMA,ETA)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M 
         COMMON /PARTCL/ P(5,50),M(3,50)

         TMP    = GAMMA*P(1,I) + ETA*P(4,I)
         P(4,I) = GAMMA*P(4,I) + ETA*P(1,I)
         P(1,I) = TMP
         P(5,I) = SQRT(P(1,I)*P(1,I) + P(2,I)*P(2,I) + P(3,I)*P(3,I))
       RETURN
       END
C............................................Lorentz boost along Y axis
       SUBROUTINE BOOSTY(I,GAMMA,ETA)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M 
         COMMON /PARTCL/ P(5,50),M(3,50)

         TMP    = GAMMA*P(2,I) + ETA*P(4,I)
         P(4,I) = GAMMA*P(4,I) + ETA*P(2,I)
         P(2,I) = TMP
         P(5,I) = SQRT(P(1,I)*P(1,I) + P(2,I)*P(2,I) + P(3,I)*P(3,I))
       RETURN
       END
C......................................................Spatial rotation
       SUBROUTINE ROT12(I,K1,K2,S,C)

         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M 
         COMMON /PARTCL/ P(5,50),M(3,50)

         TMP     = C*P(K1,I) + S*P(K2,I)
         P(K2,I) = C*P(K2,I) - S*P(K1,I)
         P(K1,I) = TMP
       RETURN
       END
C....................................Store magnitude of momentum vector
       SUBROUTINE PSET(I)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M(3,50),P(5,50),TMP
         COMMON /PARTCL/ P,M

         TMP    = P(1,I)**2 + P(2,I)**2 + P(3,I)**2
         P(5,I) = SQRT(TMP)
       RETURN
       END
C=====================TRANSFER momentums for HELAS===================#

         SUBROUTINE TRANSFER(I,Q)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         COMMON/PARTCL/P(5,50),M(3,50)
         DIMENSION Q(0:3)
                
          Q(0) = P(4,I) 
          Q(1) = P(1,I) 
          Q(2) = P(2,I)
          Q(3) = P(3,I) 

         RETURN
         END

         SUBROUTINE TRANSFER4(Q,p1,p2,p3,p4)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         DOUBLE PRECISION Q(0:3,6),p1(0:3),p2(0:3),p3(0:3),p4(0:3)
           Q(0,1)=p1(0) 
           Q(1,1)=p1(1) 
           Q(2,1)=p1(2) 
           Q(3,1)=p1(3) 
           Q(0,2)=p2(0) 
           Q(1,2)=p2(1) 
           Q(2,2)=p2(2) 
           Q(3,2)=p2(3) 
           Q(0,3)=p3(0) 
           Q(1,3)=p3(1) 
           Q(2,3)=p3(2) 
           Q(3,3)=p3(3) 
           Q(0,4)=p4(0) 
           Q(1,4)=p4(1) 
           Q(2,4)=p4(2) 
           Q(3,4)=p4(3) 
         RETURN
         END	  	 

         SUBROUTINE TRANSFER5(Q,p1,p2,p3,p4,p5)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         DOUBLE PRECISION Q(0:3,6),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
           Q(0,1)=p1(0) 
           Q(1,1)=p1(1) 
           Q(2,1)=p1(2) 
           Q(3,1)=p1(3) 
           Q(0,2)=p2(0) 
           Q(1,2)=p2(1) 
           Q(2,2)=p2(2) 
           Q(3,2)=p2(3) 
           Q(0,3)=p3(0) 
           Q(1,3)=p3(1) 
           Q(2,3)=p3(2) 
           Q(3,3)=p3(3) 
           Q(0,4)=p4(0) 
           Q(1,4)=p4(1) 
           Q(2,4)=p4(2) 
           Q(3,4)=p4(3) 
           Q(0,5)=p5(0) 
           Q(1,5)=p5(1) 
           Q(2,5)=p5(2) 
           Q(3,5)=p5(3) 
         RETURN
         END	  	 

         SUBROUTINE TRANSFER6(Q,p1,p2,p3,p4,p5,p6)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         real*8 Q(0:3,6),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
           Q(0,1)=p1(0) 
           Q(1,1)=p1(1) 
           Q(2,1)=p1(2) 
           Q(3,1)=p1(3) 
           Q(0,2)=p2(0) 
           Q(1,2)=p2(1) 
           Q(2,2)=p2(2) 
           Q(3,2)=p2(3) 
           Q(0,3)=p3(0) 
           Q(1,3)=p3(1) 
           Q(2,3)=p3(2) 
           Q(3,3)=p3(3) 
           Q(0,4)=p4(0) 
           Q(1,4)=p4(1) 
           Q(2,4)=p4(2) 
           Q(3,4)=p4(3) 
           Q(0,5)=p5(0) 
           Q(1,5)=p5(1) 
           Q(2,5)=p5(2) 
           Q(3,5)=p5(3) 
           Q(0,6)=p6(0) 
           Q(1,6)=p6(1) 
           Q(2,6)=p6(2) 
           Q(3,6)=p6(3) 
         RETURN
         END	  	 

C=============================Rapidity================================#
      FUNCTION RPID(P)
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
      double precision P(0:3)

      RPID = 0.5d0*dlog((p(0)+p(3))/(p(0)-p(3)))

      return  
      end
C=========================== phi=========================#
         SUBROUTINE COLAT2(I,PH)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         COMMON/PARTCL/P(5,50),M(3,50)
         PI = DACOS(-1.d0)
    
            IF(P(1,I).NE.0.d0) GO TO 88
                 IF(P(2,I).GT.0.d0)  PH = PI/2.d0
                 IF(P(2,I).EQ.0.d0)  PH = 0.d0
                 IF(P(2,I).LT.0.d0)  PH = 3*PI/2.d0
                 RETURN

88            IF(P(1,I).LT.0.d0)THEN
                 PH = DATAN(P(2,I)/P(1,I))+PI
             ELSE
              IF(P(1,I).GT.0.d0.AND.P(2,I).LT.0.d0)THEN
                 PH = DATAN(P(2,I)/P(1,I))+2.*PI
             ELSE
                 PH= DATAN(P(2,I)/P(1,I))
             ENDIF
             ENDIF
             RETURN
             END 
C=========================== phi=========================#
         SUBROUTINE COLAT(P,PH)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         double precision P(0:3)
         PI = DACOS(-1.d0)
    
            IF(P(1).NE.0.d0) GO TO 88
                 IF(P(2).GT.0.d0)  PH = PI/2.d0
                 IF(P(2).EQ.0.d0)  PH = 0.d0
                 IF(P(2).LT.0.d0)  PH = 3*PI/2.d0
                 RETURN

88            IF(P(1).LT.0.d0)THEN
                 PH = DATAN(P(2)/P(1))+PI
             ELSE
              IF(P(1).GT.0.d0.AND.P(2).LT.0.d0)THEN
                 PH = DATAN(P(2)/P(1))+2.*PI
             ELSE
                 PH= DATAN(P(2)/P(1))
             ENDIF
             ENDIF
             RETURN
             END 
C=======================PT=========================#
c
c         FUNCTION PT(P)
c         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
c         double precision P(0:3)
c
c          PT = dsqrt(p(1)**2 + p(2)**2)
c            RETURN
c            END 
c
C=============================Rapidity================================#

      SUBROUTINE RAPIDITY(L,AYE)
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
      COMMON/PARTCL/P(5,50),M(3,50)

      aye = 0.5d0*dlog((p(4,l)+p(3,l))/(p(4,l)-p(3,l)))
                                !aye = abs(aye)
      return  
      end

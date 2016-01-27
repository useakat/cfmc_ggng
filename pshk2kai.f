c----------------------------------------------------------------------------------------------
       SUBROUTINE pshk0(SS,n,rdn,x1,x2,WPSn,Jacob,IFLG)
       DOUBLE PRECISION P(5,50),M(2,50)
       COMMON /PARTCL/P,M
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       REAL*8 PTcut,AYEcut,QCut,Pi
       COMMON /cut/PTcut,AYEcut,QCut,Pi
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG
       DOUBLE PRECISION rdn(30),KEPT2(20),KPT4(20)
       DOUBLE PRECISION spt1,spt2,E00,PZ0,cutpt,cutpt2,cuteta,cuteta2,cuteta3
       DOUBLE PRECISION SS, sh, x1,x2,pipt2,WPSn,pau2,Jacob
       
       WPSn   = 0.0D0
       IFLG = 0 
       pipt2=1.0d0
       pau2=(dsqrt(SS)/2d0)**2
       cutpt=PTcut
       cutpt2=PTcut2
       cuteta=AYEcut
       cuteta2=AYEcut2
       cuteta3=AYEcut3*0.8d0
       spt1=0.0d0 
       spt2=0.0d0 

       KETA(n)=-cuteta3 +2*cuteta3*rdn(3*n-5)
       KEPT2(n)=dlog(cutpt2**2) +(dlog(pau2) -dlog(cutpt2**2))*rdn(3*n-4)
       KPT2(n)=dexp(KEPT2(n))
       KPT(n)=dsqrt(KPT2(n))
       KPhi(n)=2.0d0*PI*rdn(3*n-3)
       spt1=spt1-KPT(n)*dcos(KPhi(n))
       spt2=spt2-KPT(n)*dsin(KPhi(n))

       if(n.gt.2) then   
          do I=1,n-2
             KETA(i)=-cuteta+2.0d0*cuteta*rdn(3*i-2)
             KEPT2(i)=dlog(cutpt**2)+(dlog(pau2)-dlog(cutpt**2))*rdn(3*i-1)
             KPT2(i)=dexp(KEPT2(i))
             KPT(i)=dsqrt(KPT2(i))
             KPhi(i)=2.0d0*PI*rdn(3*i)
             spt1=spt1-KPT(i)*dcos(KPhi(i))
             spt2=spt2-KPT(i)*dsin(KPhi(i))
          enddo
       endif

       KETA(n-1)=-cuteta2+2.0d0*cuteta2*rdn(3*n-2)
       KPT2(n-1)=spt1**2+spt2**2
       KPT(n-1)=dsqrt(spt1**2+spt2**2)

       IF(spt1.eq.0.d0) THEN 
          IF(spt2.GT.0.d0)   KPhi(n-1) = PI/2.d0
          IF(spt2.EQ.0.d0)  KPhi(n-1) = 0.d0
          IF(spt2.LT.0.d0)    KPhi(n-1) = 3.*PI/2.d0
       ELSEIF(spt1.LT.0.d0.AND.spt2.GT.0.d0) THEN
          KPhi(n-1) = DATAN(spt2/spt1)+PI
       ELSEIF(spt1.LT.0.d0.AND.spt2.LT.0.d0)THEN
          KPhi(n-1) = DATAN(spt2/spt1)+PI
       ELSEIF(spt1.GT.0.d0.AND.spt2.LT.0.d0)THEN
          KPhi(n-1) = DATAN(spt2/spt1)+2.0d0*PI
       ELSE
          KPhi(n-1)= DATAN(spt2/spt1)
       ENDIF

       if (n.gt.2) then
          do I=1,n-2
             pipt2=pipt2*KPT2(i)
          enddo
          pipt2=pipt2*KPT2(n)
       else
          pipt2 = KPT2(n)
       endif

       E00=0.0d0
       PZ0=0.0d0
       do J=1,n
          P(4,J+2)=dsqrt(KPT2(J)+M(2,J+2))*(dexp(KETA(J))+dexp(-KETA(J)))/2.0d0
          P(3,J+2)=dsqrt(KPT2(J)+M(2,J+2))*(dexp(KETA(J))-dexp(-KETA(J)))/2.0d0
          P(1,J+2)=KPT(J)*dcos(KPhi(j))
          P(2,J+2)=KPT(J)*dsin(KPhi(j)) 
          P(5,J+2)=dsqrt(P(1,J+2)**2+P(2,J+2)**2+P(3,J+2)**2) 
          E00=E00+P(4,J+2)
          PZ0=PZ0+P(3,J+2) 
       enddo

       x1=(E00+PZ0)/dsqrt(SS) 
       x2=(E00-PZ0)/dsqrt(SS) 
       sh=ss*x1*x2
       If(x1.gt.1.0d0.or.x2.gt.1.0d0) then
          IFLG=1
       endif
       P(4,1)=x1*dsqrt(ss)/2.0d0
       P(1,1)=0.d0
       P(2,1)=0.d0
       P(3,1)=x1*dsqrt(ss)/2.0d0
       P(5,1)=x1*dsqrt(ss)/2.0d0
       P(4,2)=x2*dsqrt(ss)/2.0d0
       P(1,2)=0.d0
       P(2,2)=0.d0
       P(3,2)=-x2*dsqrt(ss)/2.0d0
       P(5,2)=x2*dsqrt(ss)/2.0d0

       WPSn=(2.0d0*cuteta)**(n-2)*(2.0d0*cuteta2)*(2*cuteta3)
     # *(dlog(pau2)-dlog(cutpt**2))**(n-2)*(dlog(pau2) -dlog(cutpt2**2))*(2.*Pi)**(n-1)/4**n*2.0d0
       WPSn=WPSn*0.38937966D+9*(2.0d0*Pi)**(4.0-3.*n)/sh/ss
       Jacob=pipt2

c!!!!! In case we call pshk2 or pshk3  before pshk
        PTcut2=PTcut 
        AYEcut2=AYEcut
        AYEcut3=AYEcut
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       RETURN                                                                   
       END 


c==============================================================
       SUBROUTINE pshk2kai(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)

       DOUBLE PRECISION P(5,50),M(2,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi
       COMMON /cut/PTcut,AYEcut,QCut,Pi
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION rdn(30),rdn2(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mini2
       real*8 emax,ewgt

       integer zimp,logimp1,mimp
       common /imps/ zimp,logimp1,mimp

       WPSn   = 0.0D0
       IFLG = 0 
       IFLG2 = 0 

       do J=1,2
          rdn2(J)=rdn(3*n-1-J)
       enddo

       mtemp1=M(1,n+1)
       mtemp2=M(2,n+1)                
       mini2=Qcut**2
       if(n.gt.2) then
          mmx=dlog(mini2)+(dlog(SS) -dlog(mini2))*rdn(3*(n-1)-1)     
          M(2,n+1)=dexp(mmx)
          M(1,n+1)=dsqrt(M(2,n+1))
          WPSn1=(dlog(SS) -dlog(mini2))/(2*pi)
          Jacobn1=dexp(mmx)

          emax = dsqrt(ss)
          AYEcut3 = dlog((emax +dsqrt(emax**2 -M(2,n+1)))/M(1,n+1))
         
          CALL pshk0(SS,n-1,rdn,x1,x2,WPSn2,Jacobn2,IFLG2)
          IF(IFLG2.ne.0) then
             IFLG = 1
             return
          endif

       elseif(n.eq.2) then
          x1=Qcut**2/SS+(1.0d0-Qcut**2/SS)*rdn(1)
          x2=Qcut**2/SS/x1+(1.0d0-Qcut**2/SS/x1)*rdn(2)
          mmx=x1*x2*SS
          M(2,n)=mmx
          M(1,n)=dsqrt(M(2,n))
          P(1,1)=0.0d0
          P(2,1)=0.0d0
          P(3,1)=x1*dsqrt(SS)/2
          P(4,1)=x1*dsqrt(SS)/2
          P(5,1)=x1*dsqrt(SS)/2
          P(1,2)=0.0d0
          P(2,2)=0.0d0
          P(3,2)=-x2*dsqrt(SS)/2
          P(4,2)=x2*dsqrt(SS)/2
          P(5,2)=x2*dsqrt(SS)/2
          P(1,3)=0.0d0
          P(2,3)=0.0d0
          P(3,3)=(x1-x2)*dsqrt(SS)/2
          P(4,3)=(x1+x2)*dsqrt(SS)/2
          P(5,3)=abs(x1-x2)*dsqrt(SS)/2
          WPSn1=(1.0d0-Qcut**2/SS)*(1.0d0-Qcut**2/SS/x1)
          WPSn2=0.38937966D+9/2.0d0/(x1*x2*SS) 
          Jacobn1=1.0d0
          Jacobn2=1.0d0
       else
          print *,"ERROR, n<2 in pshk2 subroutine"
       endif
       
       M(2,50)=M(2,n+1)
       M(1,50)=M(1,n+1)
       P(1,50)=P(1,n+1)
       P(2,50)=P(2,n+1)
       P(3,50)=P(3,n+1)
       P(4,50)=P(4,n+1)
       P(5,50)=P(5,n+1)
       M(1,n+1)=mtemp1
       M(2,n+1)=mtemp2

       CALL DK2(50,n+1,n+2,rdn2,WPSn3,Jacobn3,IFLG2) 
       
       IF(IFLG2.ne.0) then
          IFLG=1
          RETURN
       endif
       
       KPT2(n-1)=P(1,n+1)**2+P(2,n+1)**2
       KPT2(n)=(P(1,n+2)**2+P(2,n+2)**2)
       KPT(n-1)=dsqrt(P(1,n+1)**2+P(2,n+1)**2)
       KPT(n)=dsqrt(P(1,n+2)**2+P(2,n+2)**2)

       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n+2,KPHI(n))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n+2,KETA(n))

       WPSn=WPSn1*WPSn2*WPSn3
       Jacobn=Jacobn1*Jacobn2*Jacobn3

       RETURN                                                                   
       END 

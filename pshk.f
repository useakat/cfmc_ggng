c----------------------------------------------------------------------------------------------
       SUBROUTINE pshk(SS,n,rdn,x1,x2,WPSn,Jacob,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG
       DOUBLE PRECISION rdn(30),KEPT2(20)
       DOUBLE PRECISION spt1,spt2,E00,PZ0,cutpt,cutpt2,cuteta,cuteta2
       DOUBLE PRECISION cuteta3
       DOUBLE PRECISION SS,sh,x1,x2,pipt2,WPSn,pau2,Jacob
       double precision jjacob(3*n-2),wwgt(3*n-2)
       real*8 ipw1,ipw2
       real*8 ptmaxl,ptminl

       IFLG = 0 
       spt1 = 0d0 
       spt2 = 0d0 

       ipw1 = 1d0
       ipw2 = 1.5d0

       pau2 = (dsqrt(SS)/2d0)**2

       do i = 1,n-1
c          if (i.ge.3) then
c             ipw1 = 0d0
c             ipw2 = 1d0
c          endif
          if ((M(3,i+2).eq.1).or.(M(1,i+2).gt.0d0)) then
             if (M(3,i+2).eq.1) then 
                cutpt = Qcut
c                cutpt = M(1,i+2)
             elseif (M(1,i+2).gt.0d0) then
                cutpt = M(1,i+2)
             endif
             if (ipw1.eq.1) then
                ptminl = dlog(cutpt**2)
                ptmaxl = dlog(pau2)
                KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                KPT2(i) = dexp(KEPT2(i)) -cutpt**2
                KPT(i) = dsqrt(KPT2(i))
                jjacob(3*i-2) = dexp(KEPT2(i))
                wwgt(3*i-2) = ptmaxl -ptminl
             elseif (ipw1.gt.1) then
                ptminl = (cutpt**2)**(1d0-ipw1)
                ptmaxl = pau2**(1d0-ipw1)
                KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                KPT2(i) = KEPT2(i)**(1d0/(1d0-ipw1)) -cutpt**2
                KPT(i) = dsqrt(KPT2(i))
                jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(1d0-ipw1)
                wwgt(3*i-2) = ptmaxl -ptminl
             endif
             if (M(3,i+2).eq.1) then
                cuteta = dacosh( dsqrt( ss/( M(2,i+2) +KPT2(i) ) ) )
             elseif (M(1,i+2).gt.0d0) then
                cuteta = min(ayecut
     &               ,dacosh( dsqrt( ss/( M(2,i+2) +KPT2(i) ) ) ))
             endif
          elseif( M(3,i+2).eq.0 ) then
                cutpt = ptcut
             if (ipw2.eq.1) then
                ptminl = dlog(cutpt**2)
                ptmaxl = dlog(pau2)
                KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                KPT2(i) = dexp(KEPT2(i))
                KPT(i) = dsqrt(KPT2(i))
                jjacob(3*i-2) = dexp(KEPT2(i))
                wwgt(3*i-2) = ptmaxl -ptminl
             elseif (ipw2.gt.1) then
                ptminl = (cutpt**2)**(1d0-ipw2)
                ptmaxl = pau2**(1d0-ipw2)
                KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                KPT2(i) = KEPT2(i)**(1d0/(1d0-ipw2))
                KPT(i) = dsqrt(KPT2(i))
c                jjacob(3*i-2) = KEPT2(i)**(ipw2/(1d0-ipw2))/(1d0-ipw2)
c                wwgt(3*i-2) = ptmaxl -ptminl
                jjacob(3*i-2) = KEPT2(i)**(ipw2/(1d0-ipw2))/(ipw2-1d0)
                wwgt(3*i-2) = -(ptmaxl -ptminl)
             endif
             cuteta = ayecut
          else
             iflg = 1
             return
          endif

          KETA(i) = -cuteta +2*cuteta*rdn(3*i-1)
          jjacob(3*i-1) = 1d0
          wwgt(3*i-1) = 2*cuteta
          
          KPhi(i) = 2*pi*rdn(3*i)
          jjacob(3*i) = 1d0
          wwgt(3*i) = 2*pi

          spt1 = spt1 -KPT(i)*dcos(KPhi(i))
          spt2 = spt2 -KPT(i)*dsin(KPhi(i))
       enddo

       KPT2(n) = spt1**2 +spt2**2
       KPT(n) = dsqrt(KPT2(n))
       if (M(3,n+2).eq.1) then
c          cuteta2 = dacosh( dsqrt( ss/( M(2,n+2) +KPT2(n) ) ) )
          cuteta2 = ayecut
       elseif( M(3,n+2).eq.0 ) then
c          cuteta2 = min( 
c     &         ayecut, dacosh( dsqrt( ss/( M(2,n+2) +KPT2(n) ) ) ) )
          cuteta2 = ayecut
       else
          iflg = 1
          return
       endif
       KETA(n) = -cuteta2 +2*cuteta2*rdn(3*n-2)
       jjacob(3*n-2) = 1d0
       wwgt(3*n-2) = 2*cuteta2

       IF(spt1.eq.0d0) THEN 
          IF(spt2.GT.0d0) KPhi(n) = PI/2d0
          IF(spt2.EQ.0d0) KPhi(n) = 0d0
          IF(spt2.LT.0d0) KPhi(n) = 3*PI/2d0
       ELSEIF(spt1.LT.0d0.AND.spt2.GT.0d0) THEN
          KPhi(n) = DATAN(spt2/spt1)+PI
       ELSEIF(spt1.LT.0d0.AND.spt2.LT.0d0)THEN
          KPhi(n) = DATAN(spt2/spt1)+PI
       ELSEIF(spt1.GT.0d0.AND.spt2.LT.0d0)THEN
          KPhi(n) = DATAN(spt2/spt1)+2d0*PI
       ELSE
          KPhi(n) = DATAN(spt2/spt1)
       ENDIF

       E00 = 0d0
       PZ0 = 0d0
       do J = 1,n
          P(4,J+2) = dsqrt( KPT2(J) +M(2,J+2) )*( dexp(KETA(J)) 
     &         +dexp(-KETA(J)) )/2d0
          P(3,J+2) = dsqrt( KPT2(J) +M(2,J+2) )*( dexp(KETA(J)) 
     &         -dexp(-KETA(J)) )/2d0
          P(1,J+2) = KPT(J)*dcos(KPhi(j))
          P(2,J+2) = KPT(J)*dsin(KPhi(j)) 
          P(5,J+2) = dsqrt( P(1,J+2)**2 +P(2,J+2)**2 +P(3,J+2)**2 ) 
          E00 = E00 +P(4,J+2)
          PZ0 = PZ0 +P(3,J+2) 
       enddo

       x1 = ( E00 +PZ0 )/dsqrt(SS) 
       x2 = ( E00 -PZ0 )/dsqrt(SS) 
       sh = ss*x1*x2
       If((x1.gt.1d0).or.(x2.gt.1d0)) then
          IFLG = 1
       endif

       P(4,1) = x1*dsqrt(ss)/2d0
       P(1,1) = 0d0
       P(2,1) = 0d0
       P(3,1) = x1*dsqrt(ss)/2d0
       P(5,1) = x1*dsqrt(ss)/2d0
       P(4,2) = x2*dsqrt(ss)/2d0
       P(1,2) = 0d0
       P(2,2) = 0d0
       P(3,2) = -x2*dsqrt(ss)/2d0
       P(5,2) = x2*dsqrt(ss)/2d0

       wgt = 1d0
       jacob = 1d0
       do i = 1,3*n-2
          wgt = wgt*wwgt(i)
          jacob = jacob*jjacob(i)
       enddo
       wpsn = wgt*2*pi/ss/(4*pi)**(2*n-2)/(2*pi)**(n-1)/(2*sh)
     &      *0.38937966d9

       RETURN                                                                   
       END 

c==============================================================
       SUBROUTINE pshk2(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION rdn(30),rdn2(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mini2,max2
       double precision K(0:3),KK(0:3,2),mmass(2)
       real*8 ipw

       WPSn = 0d0
       IFLG = 0 
       IFLG2 = 0 

       ipw = 1.1d0

       do J = 1,2
          rdn2(J) = rdn(3*n-1-J)
       enddo

       mtemp1 = M(1,n)
       mtemp2 = M(2,n)       
       if(n.gt.2) then
          if (ipw.eq.1) then
             mini2 = dlog(Qcut**2)
             max2 = dlog(Qmax**2)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n) = 1
             M(2,n) = dexp(mmx)
             M(1,n) = dsqrt(M(2,n))
             WPSn1 = ( max2 -mini2 )/2d0/Pi
             Jacobn1 = dexp(mmx)
         elseif (ipw.gt.1) then
             mini2 = (Qmax**2)**(1d0-ipw)
             max2 = (Qcut**2)**(1d0-ipw)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n) = 1
             M(2,n) = mmx**(1d0/(1d0-ipw))
             M(1,n) = dsqrt(M(2,n))
             WPSn1 = ( max2 -mini2 )/2d0/Pi
             Jacobn1 = M(2,n)**ipw/(ipw-1d0)
         endif
c     AYEcut2=dlog(2.0d0*dsqrt(SS)/M(1,n+1))
          AYEcut2 = dlog((dsqrt(SS) +dsqrt(SS -M(2,n)))/M(1,n))
          CALL pshk(SS,n-1,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          IF(IFLG.ne.0) return
c       elseif(n.eq.2) then
c          x1=Qcut**2/SS+(1.0d0-Qcut**2/SS)*rdn(1)
c          x2=Qcut**2/SS/x1+(1.0d0-Qcut**2/SS/x1)*rdn(2)
c          mmx=x1*x2*SS
c          M(2,n)=mmx
c          M(1,n)=dsqrt(M(2,n))
c          P(1,1)=0.0d0
c          P(2,1)=0.0d0
c          P(3,1)=x1*dsqrt(SS)/2
c          P(4,1)=x1*dsqrt(SS)/2
c          P(5,1)=x1*dsqrt(SS)/2
c          P(1,2)=0.0d0
c          P(2,2)=0.0d0
c          P(3,2)=-x2*dsqrt(SS)/2
c          P(4,2)=x2*dsqrt(SS)/2
c          P(5,2)=x2*dsqrt(SS)/2
c          P(1,3)=0.0d0
c          P(2,3)=0.0d0
c          P(3,3)=(x1-x2)*dsqrt(SS)/2
c          P(4,3)=(x1+x2)*dsqrt(SS)/2
c          P(5,3)=abs(x1-x2)*dsqrt(SS)/2
c          WPSn1=(1.0d0-Qcut**2/SS)*(1.0d0-Qcut**2/SS/x1)
c          WPSn2=0.38937966D+9/2.0d0/(x1*x2*SS) 
c          Jacobn1=1.0d0
c          Jacobn2=1.0d0
       else
          write(99,*) "ERROR: n<2 (pshk2)"
          stop
       endif

       do i = 1,5
          P(i,n+2) = P(i,n+1)
       enddo

       M(2,50) = M(2,n)
       M(1,50) = M(1,n)
       do i = 1,5
          P(i,50) = P(i,n)
       enddo
       M(1,n) = mtemp1
       M(2,n) = mtemp2
       CALL DK2(50,n,n+1,rdn2,WPSn3,Jacobn3,IFLG2) 
c       mmass(1) = M(1,n+1)
c       mmass(2) = M(1,n+2)
c       call ps2bd_in_pshk(rdn2,P(1,50),mmass,P(1,n+1),wpsn3,jacobn3
c     &      ,iflg2)
       IF(IFLG2.ne.0) then
          IFLG=1
          RETURN
       endif

       KPT2(n-2) = P(1,n)**2 +P(2,n)**2
       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT(n-2) = dsqrt(KPT2(n-2))
       KPT(n-1) = dsqrt(KPT2(n-1))
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n+1,KETA(n-1))

       WPSn = WPSn1*WPSn2*WPSn3
       Jacobn = Jacobn1*Jacobn2*Jacobn3

       RETURN                                                                   
       END 

c==============================================================
       SUBROUTINE pshk2_2(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2,IFLG3
       DOUBLE PRECISION rdn(30),rdn2(2),rdn3(2)
       DOUBLE PRECISION SS,x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4
       DOUBLE PRECISION mtemp1_1,mtemp1_2,mtemp2_1,mtemp2_2,mmx,mmx1
       double precision mmx2,mini2,max2,mmass(4)

       WPSn = 0d0
       IFLG = 0 
       IFLG2 = 0 

       do J = 1,2
          rdn2(J) = rdn(3*n-1-J)
          rdn3(J) = rdn(3*n-3-J)
       enddo

       mtemp1_1 = M(1,n)
       mtemp1_2 = M(2,n)                
       mtemp2_1 = M(1,n-1)
       mtemp2_2 = M(2,n-1)                
       mini2 = Qcut**2
       max2 = Qmax**2

       if(n.ge.4) then
          mmx1 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-3)     
          mmx2 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-4)     
c          mmx2=dlog(8d0**2)+(dlog(16d0**2)-dlog(8d0**2))*rdn(3*(n-1)-4)     
          M(3,n) = 1
          M(2,n) = dexp(mmx1)
          M(1,n) = dsqrt(M(2,n))
          M(3,n-1) = 1
          M(2,n-1) = dexp(mmx2)
          M(1,n-1) = dsqrt(M(2,n-1))

          WPSn1 = (( dlog(max2) -dlog(mini2) )/2d0/pi)**2
c          WPSn1=((dlog(max2)- dlog(mini2))/2.0d0/Pi)*((dlog(16d0**2)- dlog(8d0**2))/2.0d0/Pi)
          Jacobn1 = dexp(mmx1)*dexp(mmx2)

c          AYEcut2 = dlog(( dsqrt(SS) +dsqrt( SS -M(2,n) ) )/M(1,n))
c          AYEcut3 = dlog(( dsqrt(SS) +dsqrt( SS -M(2,n-1) ) )/M(1,n-1))

          CALL pshk(SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          IF(IFLG.ne.0) then
             return
          endif
       else
          print *,"ERROR, n < 4 in pshk2_2 subroutine"
          stop
       endif

       M(1,50) = M(1,n)
       M(2,50) = M(2,n)
       P(1,50) = P(1,n)
       P(2,50) = P(2,n)
       P(3,50) = P(3,n)
       P(4,50) = P(4,n)
       P(5,50) = P(5,n)
       M(1,n) = mtemp1_1
       M(2,n) = mtemp1_2

       M(1,49) = M(1,n-1)
       M(2,49) = M(2,n-1)
       P(1,49) = P(1,n-1)
       P(2,49) = P(2,n-1)
       P(3,49) = P(3,n-1)
       P(4,49) = P(4,n-1)
       P(5,49) = P(5,n-1)
       M(1,n-1) = mtemp2_1
       M(2,n-1) = mtemp2_2

c       CALL DK2(50,n+1,n+2,rdn2,WPSn3,Jacobn3,IFLG2) 
       mmass(1) = M(1,n+1)
       mmass(2) = M(1,n+2)
       call ps2bd_in_pshk(rdn2,P(1,50),mmass,P(1,n+1),wpsn3,jacobn3
     &      ,iflg2)
       IF(IFLG2.ne.0) then
          IFLG = 1
          RETURN
       endif

c       CALL DK2(49,n-1,n,rdn3,WPSn4,Jacobn4,IFLG3) 
       mmass(3) = M(1,n-1)
       mmass(4) = M(1,n)
       call ps2bd_in_pshk(rdn3,P(1,49),mmass(3),P(1,n-1),wpsn4,jacobn4
     &      ,iflg3)
       IF(IFLG3.ne.0) then
          IFLG = 1
          RETURN
       endif

       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT2(n) = P(1,n+2)**2 +P(2,n+2)**2
       KPT(n-1) = dsqrt(KPT2(n-1))
       KPT(n) = dsqrt(KPT2(n))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n+2,KPHI(n))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n+2,KETA(n))

       WPSn=WPSn1*WPSn2*WPSn3*WPSn4
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4

       RETURN                                                                   
       END 

c=========================================================
       SUBROUTINE pshk3(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION  rdn(30),rdn31(2),rdn32(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4,WPSn5
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4,Jacobn5
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mmy,mini2

       WPSn   = 0.0D0
       IFLG = 0 
       IFLG2 = 0 

       do J=1,2
          rdn31(J)=rdn(3*n-1-J)
       enddo
       do J=1,2
          rdn32(J)=rdn(3*n-3-J)
       enddo
 
       mtemp1=M(1,n)
       mtemp2=M(2,n)                
       mini2=Qcut**2
       if(n.gt.3) then
          mmx=dlog(mini2)+(dlog((dsqrt(SS))**2)- dlog(mini2))*rdn(3*(n-1)-4)   
          M(2,n)=dexp(mmx)
          M(1,n)=dsqrt(M(2,n))
          WPSn1=(dlog((dsqrt(SS))**2)- dlog(mini2))/2.0d0/Pi
          Jacobn1=dexp(mmx)
          AYEcut2=dlog(2.0d0*dsqrt(SS)/M(1,n))
          CALL pshk(SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
       elseif(n.eq.3) then
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
          P(1,n)=0.0d0
          P(2,n)=0.0d0
          P(3,n)=(x1-x2)*dsqrt(SS)/2
          P(4,n)=(x1+x2)*dsqrt(SS)/2
          P(5,n)=abs(x1-x2)*dsqrt(SS)/2
          WPSn1=(1.0d0-Qcut**2/SS)*(1.0d0-Qcut**2/SS/x1)
          WPSn2=0.38937966D+9/2.0d0/(x1*x2*SS)
          Jacobn1=1.0d0
          Jacobn2=1.0d0
       else
          print *,"ERROR, n<3 in pshk3 subroutine"
          stop
       endif
       
       M(2,50)=M(2,n)
       M(1,50)=M(1,n)
       P(1,50)=P(1,n)
       P(2,50)=P(2,n)
       P(3,50)=P(3,n)
       P(4,50)=P(4,n)
       P(5,50)=P(5,n)
       M(1,n)=mtemp1
       M(2,n)=mtemp2
       
       mmy=dlog(mini2)+(dlog(M(2,50))- dlog(mini2))*rdn(3*(n-1)-3)   

       M(2,49)=dexp(mmy)
       M(1,49)=dsqrt(M(2,49))
       WPSn3=(dlog(M(2,50))- dlog(mini2))/2.0d0/Pi
       Jacobn3=dexp(mmy)

       CALL DK2(50,n,49,rdn31,WPSn4,Jacobn4,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif
       CALL DK2(49,n+1,n+2,rdn32,WPSn5,Jacobn5,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif

       KPT2(n-2)=P(1,n)**2+P(2,n)**2
       KPT2(n-1)=(P(1,n+1)**2+P(2,n+1)**2)
       KPT2(n)=(P(1,n+2)**2+P(2,n+2)**2)
       KPT(n-2)=dsqrt(P(1,n)**2+P(2,n)**2)
       KPT(n-1)=dsqrt(P(1,n+1)**2+P(2,n+1)**2)
       KPT(n)=dsqrt(P(1,n+2)**2+P(2,n+2)**2)
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n+2,KPHI(n))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n+2,KETA(n))
       WPSn=WPSn1*WPSn2*WPSn3*WPSn4*WPSn5
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4*Jacobn5
       RETURN                                                                   
       END 
c=========================================================
       SUBROUTINE pshk41(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION  rdn(30),rdn31(2),rdn32(2),rdn33(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4,WPSn5,WPSn6,WPSn7
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4,Jacobn5,Jacobn6,Jacobn7
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mmy,mmz,mini2


       WPSn   = 0.0D0
       IFLG = 0 
       IFLG2 = 0 

       do J=1,2
       rdn31(J)=rdn(3*n-1-J)
       rdn32(J)=rdn(3*n-3-J)
       rdn33(J)=rdn(3*n-5-J) 
       enddo
 
       mtemp1=M(1,n-1)
       mtemp2=M(2,n-1)            
       mini2=Qcut**2

       if(n.gt.4) then
       mmx=dlog(mini2)+(dlog((dsqrt(SS))**2)- dlog(mini2))*rdn(3*n-10)   
       M(2,n-1)=dexp(mmx)
       M(1,n-1)=dsqrt(M(2,n-1))
       WPSn1=(dlog((dsqrt(SS))**2)- dlog(mini2))/2.0d0/Pi
       Jacobn1=dexp(mmx)
       AYEcut2=dlog(2.0d0*dsqrt(SS)/M(1,n-1))
       CALL pshk(SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
       elseif(n.eq.4) then
       x1=Qcut**2/SS+(1.0d0-Qcut**2/SS)*rdn(1)
       x2=Qcut**2/SS/x1+(1.0d0-Qcut**2/SS/x1)*rdn(2)
       mmx=x1*x2*SS
       M(2,n-1)=mmx
       M(1,n-1)=dsqrt(M(2,n-1))
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
       P(1,n-1)=0.0d0
       P(2,n-1)=0.0d0
       P(3,n-1)=(x1-x2)*dsqrt(SS)/2
       P(4,n-1)=(x1+x2)*dsqrt(SS)/2
       P(5,n-1)=abs(x1-x2)*dsqrt(SS)/2
       WPSn1=(1.0d0-Qcut**2/SS)*(1.0d0-Qcut**2/SS/x1)
       WPSn2=0.38937966D+9/2.0d0/(x1*x2*SS)
       Jacobn1=1.0d0
       Jacobn2=1.0d0
       else
          print *,"ERROR, n<4 in pshk41 subroutine"
          stop
       endif

       M(2,50)=M(2,n-1)
       M(1,50)=M(1,n-1)
       P(1,50)=P(1,n-1)
       P(2,50)=P(2,n-1)
       P(3,50)=P(3,n-1)
       P(4,50)=P(4,n-1)
       P(5,50)=P(5,n-1)
       M(1,n-1)=mtemp1
       M(2,n-1)=mtemp2
       
       mmy=dlog(mini2)+(dlog(M(2,50))- dlog(mini2))*rdn(3*n-9)   
       M(2,49)=dexp(mmy)
       M(1,49)=dsqrt(M(2,49))
       WPSn3=(dlog(M(2,50))- dlog(mini2))/2.0d0/Pi
       Jacobn3=dexp(mmy)
       mmz=dlog(mini2)+(dlog(M(2,49))- dlog(mini2))*rdn(3*n-8)   
       M(2,48)=dexp(mmz)
       M(1,48)=dsqrt(M(2,48))
       WPSn4=(dlog(M(2,49))- dlog(mini2))/2.0d0/Pi
       Jacobn4=dexp(mmz)
       CALL DK2(50,n-1,49,rdn31,WPSn5,Jacobn5,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif
       CALL DK2(49,n,48,rdn32,WPSn6,Jacobn6,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif
       CALL DK2(48,n+1,n+2,rdn33,WPSn7,Jacobn7,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif

       KPT2(n-3)=P(1,n-1)**2+P(2,n-1)**2
       KPT2(n-2)=P(1,n)**2+P(2,n)**2
       KPT2(n-1)=(P(1,n+1)**2+P(2,n+1)**2)
       KPT2(n)=(P(1,n+2)**2+P(2,n+2)**2)
       KPT(n-3)=dsqrt(P(1,n-1)**2+P(2,n-1)**2)
       KPT(n-2)=dsqrt(P(1,n)**2+P(2,n)**2)
       KPT(n-1)=dsqrt(P(1,n+1)**2+P(2,n+1)**2)
       KPT(n)=dsqrt(P(1,n+2)**2+P(2,n+2)**2)
       CALL COLAT2(n-1,KPHI(n-3))
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n+2,KPHI(n))
       CALL RAPIDITY(n-1,KETA(n-3))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n+2,KETA(n))
       WPSn=WPSn1*WPSn2*WPSn3*WPSn4*WPSn5*WPSn6*WPSn7
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4*Jacobn5*Jacobn6*Jacobn7
       RETURN                                                                   
       END 
c=========================================================
       SUBROUTINE pshk42(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION  rdn(30),rdn31(2),rdn32(2),rdn33(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4,WPSn5,WPSn6,WPSn7
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4,Jacobn5,Jacobn6,Jacobn7
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mmy,mmz,mini2


       WPSn   = 0.0D0
       IFLG = 0 
       IFLG2 = 0 

       do J=1,2
       rdn31(J)=rdn(3*n-1-J)
       rdn32(J)=rdn(3*n-3-J)
       rdn33(J)=rdn(3*n-5-J) 
       enddo
 
       mtemp1=M(1,n-1)
       mtemp2=M(2,n-1)            
       mini2=Qcut**2

       if(n.gt.4) then
       mmx=dlog(mini2)+(dlog((dsqrt(SS))**2)- dlog(mini2))*rdn(3*n-10)   
       M(2,n-1)=dexp(mmx)
       M(1,n-1)=dsqrt(M(2,n-1))
       WPSn1=(dlog((dsqrt(SS))**2)- dlog(mini2))/2.0d0/Pi
       Jacobn1=dexp(mmx)
       AYEcut2=dlog(2.0d0*dsqrt(SS)/M(1,n-1))
       CALL pshk(SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
       elseif(n.eq.4) then
       x1=Qcut**2/SS+(1.0d0-Qcut**2/SS)*rdn(1)
       x2=Qcut**2/SS/x1+(1.0d0-Qcut**2/SS/x1)*rdn(2)
       mmx=x1*x2*SS
       M(2,n-1)=mmx
       M(1,n-1)=dsqrt(M(2,n-1))
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
       P(1,n-1)=0.0d0
       P(2,n-1)=0.0d0
       P(3,n-1)=(x1-x2)*dsqrt(SS)/2
       P(4,n-1)=(x1+x2)*dsqrt(SS)/2
       P(5,n-1)=abs(x1-x2)*dsqrt(SS)/2
       WPSn1=(1.0d0-Qcut**2/SS)*(1.0d0-Qcut**2/SS/x1)
       WPSn2=0.38937966D+9/2.0d0/(x1*x2*SS)
       Jacobn1=1.0d0
       Jacobn2=1.0d0
       else
          print *,"ERROR, n<4 in pshk41 subroutine"
          stop
       endif

       M(2,50)=M(2,n-1)
       M(1,50)=M(1,n-1)
       P(1,50)=P(1,n-1)
       P(2,50)=P(2,n-1)
       P(3,50)=P(3,n-1)
       P(4,50)=P(4,n-1)
       P(5,50)=P(5,n-1)
       M(1,n-1)=mtemp1
       M(2,n-1)=mtemp2
       
       mmy=dlog(mini2)+(dlog(M(2,50))- dlog(mini2))*rdn(3*n-9)   
       M(2,49)=dexp(mmy)
       M(1,49)=dsqrt(M(2,49))
       WPSn3=(dlog(M(2,50))- dlog(mini2))/2.0d0/Pi
       Jacobn3=dexp(mmy)
       mmz=dlog(mini2)+((dlog(M(2,50)-M(2,49)))-dlog(mini2))*rdn(3*n-8)   
       M(2,48)=dexp(mmz)
       M(1,48)=dsqrt(M(2,48))
       WPSn4=((dlog(M(2,50)-M(2,49)))-dlog(mini2))/2.0d0/Pi
       Jacobn4=dexp(mmz)

       CALL DK2(50,49,48,rdn31,WPSn5,Jacobn5,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif
       CALL DK2(49,n-1,n,rdn32,WPSn6,Jacobn6,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif
       CALL DK2(48,n+1,n+2,rdn33,WPSn7,Jacobn7,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif

       KPT2(n-3)=P(1,n-1)**2+P(2,n-1)**2
       KPT2(n-2)=P(1,n)**2+P(2,n)**2
       KPT2(n-1)=(P(1,n+1)**2+P(2,n+1)**2)
       KPT2(n)=(P(1,n+2)**2+P(2,n+2)**2)
       KPT(n-3)=dsqrt(P(1,n-1)**2+P(2,n-1)**2)
       KPT(n-2)=dsqrt(P(1,n)**2+P(2,n)**2)
       KPT(n-1)=dsqrt(P(1,n+1)**2+P(2,n+1)**2)
       KPT(n)=dsqrt(P(1,n+2)**2+P(2,n+2)**2)
       CALL COLAT2(n-1,KPHI(n-3))
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n+2,KPHI(n))
       CALL RAPIDITY(n-1,KETA(n-3))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n+2,KETA(n))
       WPSn=WPSn1*WPSn2*WPSn3*WPSn4*WPSn5*WPSn6*WPSn7
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4*Jacobn5*Jacobn6*Jacobn7
       RETURN                                                                   
       END 
c=========================================================
       SUBROUTINE pshk43(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION  rdn(30),rdn31(2),rdn32(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4,WPSn5
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4,Jacobn5
       DOUBLE PRECISION mtemp1,mtemp2,mtempb1,mtempb2,mmx,mmy,mini2


       WPSn   = 0.0D0
       IFLG = 0 
       IFLG2 = 0 

       do J=1,2
       rdn31(J)=rdn(3*n-1-J)
       rdn32(J)=rdn(3*n-3-J)
       enddo
 
       mtemp1=M(1,n-1)
       mtemp2=M(2,n-1)            
       mtempb1=M(1,n)
       mtempb2=M(2,n) 
       mini2=Qcut**2


       mmx=dlog(mini2)+(dlog((dsqrt(SS))**2)- dlog(mini2))*rdn(3*n-6)   
       M(2,n-1)=dexp(mmx)
       M(1,n-1)=dsqrt(M(2,n-1))
       WPSn1=(dlog((dsqrt(SS))**2)- dlog(mini2))/2.0d0/Pi
       Jacobn1=dexp(mmx)
       AYEcut2=dlog(2.0d0*dsqrt(SS)/M(1,n-1))
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       PTcut2=0.0001d0

       mmy=dlog(mini2)+(dlog((dsqrt(SS))**2)- dlog(mini2))*rdn(3*n-7)   
       M(2,n)=dexp(mmy)
       M(1,n)=dsqrt(M(2,n))
       WPSn2=(dlog((dsqrt(SS))**2)- dlog(mini2))/2.0d0/Pi
       Jacobn2=dexp(mmy)
       AYEcut3=dlog(2.0d0*dsqrt(SS)/M(1,n))

       CALL pshk0(SS,n-2,rdn,x1,x2,WPSn3,Jacobn3,IFLG)

       M(2,49)=M(2,n-1)
       M(1,49)=M(1,n-1)
       P(1,49)=P(1,n-1)
       P(2,49)=P(2,n-1)
       P(3,49)=P(3,n-1)
       P(4,49)=P(4,n-1)
       P(5,49)=P(5,n-1)
       M(1,n-1)=mtemp1
       M(2,n-1)=mtemp2
       M(2,50)=M(2,n)
       M(1,50)=M(1,n)
       P(1,50)=P(1,n)
       P(2,50)=P(2,n)
       P(3,50)=P(3,n)
       P(4,50)=P(4,n)
       P(5,50)=P(5,n)
       M(1,n)=mtempb1
       M(2,n)=mtempb2


       CALL DK2(49,n-1,n,rdn31,WPSn4,Jacobn4,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif
       CALL DK2(50,n+1,n+2,rdn32,WPSn5,Jacobn5,IFLG2) 
       IF(IFLG2.ne.0) then
       IFLG=1
       RETURN
       endif

       KPT2(n-3)=P(1,n-1)**2+P(2,n-1)**2
       KPT2(n-2)=P(1,n)**2+P(2,n)**2
       KPT2(n-1)=(P(1,n+1)**2+P(2,n+1)**2)
       KPT2(n)=(P(1,n+2)**2+P(2,n+2)**2)
       KPT(n-3)=dsqrt(P(1,n-1)**2+P(2,n-1)**2)
       KPT(n-2)=dsqrt(P(1,n)**2+P(2,n)**2)
       KPT(n-1)=dsqrt(P(1,n+1)**2+P(2,n+1)**2)
       KPT(n)=dsqrt(P(1,n+2)**2+P(2,n+2)**2)
       CALL COLAT2(n-1,KPHI(n-3))
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n+2,KPHI(n))
       CALL RAPIDITY(n-1,KETA(n-3))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n+2,KETA(n))
       WPSn=WPSn1*WPSn2*WPSn3*WPSn4*WPSn5
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4*Jacobn5
       RETURN                                                                   
       END 
c====================================================
      SUBROUTINE DK2(IY,I1,I2,ZZ,WPS2,Jacob2,IFLG2) 
cc  I2 carries fz I1 carries 1-fz
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
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
      lol=dlog(lx/(1.0d0-lx))
      uuz=lol+2*abs(lol)*ZZ(1)
      fz=dexp(uuz)/(1+dexp(uuz))
      CV=(2*fz-1.0d0-M2SQ/MYSQ+M1SQ/MYSQ)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/lamd   
      IF(abs(CV).GT.1.) then 
      IFLG2=1
      RETURN
      ENDIF  

      WPS2=6.0d0*Pi*abs(lol)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/(2.d0*pi)**(3*2.d0-4.0d0)
      Jacob2=1.0d0/6.0d0/(1.0d0/fz+1.0d0/(1.0d0-fz))

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

c----------------------------------------------------------------------------------------------
       SUBROUTINE pshkpos(ipos,delyflag,SS,n,rdn,x1,x2,WPSn,Jacob,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20),DKETA(20)
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG
       DOUBLE PRECISION rdn(30),KEPT2(20),KKPT2(20)
       DOUBLE PRECISION spt1,spt2,E00,PZ0,cutpt,cutpt2,cuteta(20)
       DOUBLE PRECISION SS,sh,x1,x2,pipt2,WPSn,pau2,Jacob
       double precision jjacob(3*n-2),wwgt(3*n-2)
       real*8 ipw1,ipw2
       real*8 ptmaxl,ptminl,rnd12,cutps,cutps_f
       integer idivps,ipos,delyflag
       real*8 M_tmp(3),KPT2_tmp,KPhi_tmp,cuteta_tmp

       IFLG = 0 
       spt1 = 0d0 
       spt2 = 0d0 
       idivps = 0
       cutps_f = 2d0
       ipw1 = 1d0
       ipw2 = 1.5d0

       pau2 = (dsqrt(SS)/2d0)**2

       if (ipos.gt.0) then ! Exchange the n-th and (n-ipos)-th momentum.
          do i = 1,3
             M_tmp(i) = M(i,2+n-ipos)
             M(i,2+n-ipos) = M(i,2+n)
             M(i,2+n) = M_tmp(i)
          enddo
       endif

       do i = 1,n-1
c          if (i.ge.3) then
c             ipw1 = 0d0
c             ipw2 = 1d0
c          endif
          if ((M(3,i+2).eq.1).or.(M(1,i+2).gt.0d0)) then
             if (M(3,i+2).eq.1) then 
                cutpt = Qcut
                if (idivps.eq.1) then
                   rnd12 = rdn(3*i-2)
                   cutps = cutps_f*Qcut
                   if (rnd12.lt.0.5d0) then
                      ipw1 = 2d0
                      ptminl = (cutpt**2)**(1d0-ipw1)
                      ptmaxl = (cutps**2)**(1d0-ipw1)
                      KEPT2(i) = ptminl +( ptmaxl -ptminl )*rnd12*2d0
                      KPT2(i) = (KEPT2(i))**(1d0/(1d0-ipw1)) -cutpt**2
                      KPT(i) = dsqrt(KPT2(i))
                      jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(ipw1-1d0)
                      wwgt(3*i-2) = -(ptmaxl -ptminl)*2d0
                   elseif (rnd12.ge.0.5d0) then
                      ipw1 = 1d0
                      ptminl = dlog(cutps**2)
                      ptmaxl = dlog(pau2)
                      KEPT2(i) = ptminl +( ptmaxl -ptminl )*(1d0-(1d0-rnd12)*2d0)
                      KPT2(i) = dexp(KEPT2(i)) -cutpt**2
                      KPT(i) = dsqrt(KPT2(i))
                      jjacob(3*i-2) = dexp(KEPT2(i))
                      wwgt(3*i-2) = (ptmaxl -ptminl)*2d0
                   endif
               else
                  ipw1 = 1d0
                  if (ipw1.eq.0) then
                     ptminl = cutpt**2
                     ptmaxl = pau2
                     KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                     KPT2(i) = KEPT2(i) -cutpt**2
                     KPT(i) = dsqrt(KPT2(i))
                     jjacob(3*i-2) = 1d0
                     wwgt(3*i-2) = ptmaxl -ptminl
                  elseif (ipw1.eq.1) then
                     ptminl = dlog(cutpt**2)
                     ptmaxl = dlog(pau2)
                     KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                     KPT2(i) = dexp(KEPT2(i)) -cutpt**2
                     KPT(i) = dsqrt(KPT2(i))
                     jjacob(3*i-2) = dexp(KEPT2(i))
                     wwgt(3*i-2) = ptmaxl -ptminl
                  elseif (ipw1.ne.1) then
                     ptminl = (cutpt**2)**(1d0-ipw1)
                     ptmaxl = pau2**(1d0-ipw1)
                     KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                     KPT2(i) = (KEPT2(i))**(1d0/(1d0-ipw1)) -cutpt**2
                     KPT(i) = dsqrt(KPT2(i))
                     if (ipw1.lt.1) then
                        jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(1d0-ipw1)
                        wwgt(3*i-2) = (ptmaxl -ptminl)
                     elseif (ipw1.gt.1) then
                        jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(ipw1-1d0)
                        wwgt(3*i-2) = -(ptmaxl -ptminl)
                     endif
                  endif
               endif
             elseif (M(1,i+2).gt.0d0) then
                cutpt = M(1,i+2)
                if (ipw1.eq.0) then
                   ptminl = cutpt**2
                   ptmaxl = pau2
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = KEPT2(i) -cutpt**2
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = 1d0
                   wwgt(3*i-2) = ptmaxl -ptminl
                elseif (ipw1.eq.1) then
                   ptminl = dlog(cutpt**2)
                   ptmaxl = dlog(pau2)
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = dexp(KEPT2(i)) -cutpt**2
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = dexp(KEPT2(i))
                   wwgt(3*i-2) = ptmaxl -ptminl
                elseif (ipw1.ne.1) then
                   ptminl = (cutpt**2)**(1d0-ipw1)
                   ptmaxl = pau2**(1d0-ipw1)
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = (KEPT2(i))**(1d0/(1d0-ipw1)) -cutpt**2
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(ipw1-1d0)
                   wwgt(3*i-2) = -(ptmaxl -ptminl)
                endif
             endif
             if (M(3,i+2).eq.1) then
c                cuteta = dacosh( dsqrt( ss/( M(2,i+2) +KPT2(i) ) ) )
                cuteta(i) = ayecut
             elseif (M(1,i+2).gt.0d0) then
                cuteta(i) = min(ayecut
c     &               ,dacosh( dsqrt( ss/( M(2,i+2) +KPT2(i) ) ) ))
     &               ,dacosh( dsqrt( ss/( M(2,i+2) ) ) ))
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
                jjacob(3*i-2) = KEPT2(i)**(ipw2/(1d0-ipw2))/(ipw2-1d0)
                wwgt(3*i-2) = -(ptmaxl -ptminl)
             endif
             cuteta(i) = ayecut
          else
             iflg = 1
             return
          endif

c          KETA(i) = -cuteta +2*cuteta*rdn(3*i-1)
c          jjacob(3*i-1) = 1d0
c          wwgt(3*i-1) = 2*cuteta
          
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
          cuteta(n) = ayecut
       elseif( M(3,n+2).eq.0 ) then
          if (M(1,n+2).eq.0) then
             cuteta(n) = ayecut
          else
             cuteta(n) = min( 
c     &         ayecut, dacosh( dsqrt( ss/( M(2,n+2) +KPT2(n) ) ) ) )
     &            ayecut, dacosh( dsqrt( ss/( M(2,n+2) ) ) ) )
          endif         
       else
          iflg = 1
          return
       endif
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

       if (ipos.gt.0) then  ! The (n-ipos)-th momentum is generated last.
          do i = 1,3
             M_tmp(i) = M(i,2+n-ipos)
             M(i,2+n-ipos) = M(i,2+n)
             M(i,2+n) = M_tmp(i)
          enddo
          KPT2_tmp = KPT2(n-ipos)
          KPT2(n-ipos) = KPT2(n)
          KPT2(n) = KPT2_tmp
          KPT(n-ipos) = dsqrt(KPT2(n-ipos))
          KPT(n) = dsqrt(KPT2(n))
          KPhi_tmp = KPhi(n-ipos)
          KPhi(n-ipos) = KPhi(n)
          KPhi(n) = KPhi_tmp
          cuteta_tmp = cuteta(n-ipos)
          cuteta(n-ipos) = cuteta(n)
          cuteta(n) = cuteta_tmp
       endif

       if (delyflag.eq.0) then
          do i = 1,n-1
             KETA(i) = -cuteta(i) +2*cuteta(i)*rdn(3*i-1)
             jjacob(3*i-1) = 1d0
             wwgt(3*i-1) = 2*cuteta(i)          
          enddo
          KETA(n) = -cuteta(n) +2*cuteta(n)*rdn(3*n-2)
          jjacob(3*n-2) = 1d0
          wwgt(3*n-2) = 2*cuteta(n)
       elseif (delyflag.eq.1) then
          KETA(1) = -cuteta(1) +2*cuteta(1)*rdn(3*1-1)
          jjacob(3*1-1) = 1d0
          wwgt(3*1-1) = 2*cuteta(1)
          do i = 2,n-1
             DKETA(i) = -cuteta(i)-KETA(i-1) +2*cuteta(i)*rdn(3*i-1)
             KETA(i) = KETA(i-1) +DKETA(i)
             jjacob(3*i-1) = 1d0
             wwgt(3*i-1) = 2*cuteta(i)
          enddo
          DKETA(n) = -cuteta(n)-KETA(n-1) +2*cuteta(n)*rdn(3*n-2)
          KETA(n) = KETA(n-1) +DKETA(n)
          jjacob(3*n-2) = 1d0
          wwgt(3*n-2) = 2*cuteta(n)
       endif

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
       SUBROUTINE pshkpos2(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
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

       ipw = 1d0

       do J = 1,2
          rdn2(J) = rdn(3*n-1-J)
       enddo

       mtemp1 = M(1,n+1)
       mtemp2 = M(2,n+1)       
       if(n.gt.2) then
          if (ipw.eq.1) then
             mini2 = dlog(Qcut**2)
             max2 = dlog(Qmax**2)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n+1) = 1
             M(2,n+1) = dexp(mmx)
             M(1,n+1) = dsqrt(M(2,n+1))
             WPSn1 = ( max2 -mini2 )/2d0/Pi
             Jacobn1 = dexp(mmx)
         elseif (ipw.gt.1) then
             mini2 = (Qmax**2)**(1d0-ipw)
             max2 = (Qcut**2)**(1d0-ipw)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n+1) = 1
             M(2,n+1) = mmx**(1d0/(1d0-ipw))
             M(1,n+1) = dsqrt(M(2,n+1)) 
             WPSn1 = -( max2 -mini2 )/2d0/Pi
             Jacobn1 = mmx**(ipw/(1d0-ipw))/(ipw-1d0)
         endif
c     AYEcut2=dlog(2.0d0*dsqrt(SS)/M(1,n+1))
          AYEcut2 = dlog((dsqrt(SS) +dsqrt(SS -M(2,n+1)))/M(1,n+1))
          CALL pshkpos(0,0,SS,n-1,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
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

c       do i = 1,5
c          P(i,n+2) = P(i,n+1)
c       enddo

       M(2,50) = M(2,n+1)
       M(1,50) = M(1,n+1)
       do i = 1,5
          P(i,50) = P(i,n+1)
       enddo
       M(1,n+1) = mtemp1
       M(2,n+1) = mtemp2
c       CALL DK2(50,n+1,n+2,rdn2,WPSn3,Jacobn3,IFLG2) 
       mmass(1) = M(1,n+1)
       mmass(2) = M(1,n+2)
       call ps2bd_in_pshk(rdn2,P(1,50),mmass,P(1,n+1),wpsn3,jacobn3
     &      ,iflg2)
       IF(IFLG2.ne.0) then
          IFLG=1
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

       WPSn = WPSn1*WPSn2*WPSn3
       Jacobn = Jacobn1*Jacobn2*Jacobn3

       RETURN                                                                   
       END 

c==============================================================
       SUBROUTINE pshkpos2_2(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2,IFLG3
       DOUBLE PRECISION rdn(30),rdn2(2),rdn3(2)
       DOUBLE PRECISION SS,x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4
       DOUBLE PRECISION mtemp1_1,mtemp1_2,mtemp2_1,mtemp2_2,mmx,mmx1
       double precision mmx2,mini2,max2,mmass(4),ipw12,ipwm,mmxminl,mmxmaxl

       WPSn = 0d0
       IFLG = 0 
       IFLG2 = 0 
       ipwm = 1d0

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
          if (ipwm.eq.0) THEN
             mmx1 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-3)     
             mmx2 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-4)     
c     mmx2=dlog(8d0**2)+(dlog(16d0**2)-dlog(8d0**2))*rdn(3*(n-1)-4)     
             M(3,n) = 1
             M(2,n) = mmx1
             M(1,n) = dsqrt(M(2,n))
             M(3,n-1) = 1
             M(2,n-1) = mmx2
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (( max2 -mini2 )/2d0/pi)**2
c     WPSn1=((dlog(max2)- dlog(mini2))/2.0d0/Pi)*((dlog(16d0**2)- dlog(8d0**2))/2.0d0/Pi)
             Jacobn1 = 1d0
          elseif (ipwm.eq.1) then
             mmx1 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-3)     
             mmx2 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-4)     
c     mmx2=dlog(8d0**2)+(dlog(16d0**2)-dlog(8d0**2))*rdn(3*(n-1)-4)     
             M(3,n) = 1
             M(2,n) = dexp(mmx1)
             M(1,n) = dsqrt(M(2,n))
             M(3,n-1) = 1
             M(2,n-1) = dexp(mmx2)
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (( dlog(max2) -dlog(mini2) )/2d0/pi)**2
c     WPSn1=((dlog(max2)- dlog(mini2))/2.0d0/Pi)*((dlog(16d0**2)- dlog(8d0**2))/2.0d0/Pi)
             Jacobn1 = dexp(mmx1)*dexp(mmx2)
          elseif(ipwm.ne.1) then
             mmxminl = mini2**(1d0-ipwm)
             mmxmaxl = max2**(1d0-ipwm)
             mmx1 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-3)     
             M(2,n) = (mmx1)**(1d0/(1d0-ipwm))
             M(1,n) = dsqrt(M(2,n))
             M(3,n) = 1
             mmx2 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-4)     
             M(2,n-1) = (mmx2)**(1d0/(1d0-ipwm))
             M(1,n-1) = dsqrt(M(2,n-1))
             M(3,n-1) = 1
             Jacobn1 = mmx1**(ipwm/(1d0-ipwm))/(ipwm-1d0)
     &            *mmx2**(ipwm/(1d0-ipwm))/(ipwm-1d0)
             WPSn1 = ((mmxmaxl -mmxminl)/(2*pi))**2
          endif


c          AYEcut2 = dlog(( dsqrt(SS) +dsqrt( SS -M(2,n) ) )/M(1,n))
c          AYEcut3 = dlog(( dsqrt(SS) +dsqrt( SS -M(2,n-1) ) )/M(1,n-1))

c          CALL pshk(SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          CALL pshkpos(0,0,SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
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

       KPT2(n) = P(1,n+2)**2 +P(2,n+2)**2
       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT2(n-2) = P(1,n)**2 +P(2,n)**2
       KPT2(n-3) = P(1,n-1)**2 +P(2,n-1)**2
       KPT(n) = dsqrt(KPT2(n))
       KPT(n-1) = dsqrt(KPT2(n-1))
       KPT(n-2) = dsqrt(KPT2(n-2))
       KPT(n-3) = dsqrt(KPT2(n-3))
       CALL COLAT2(n+2,KPHI(n))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n-1,KPHI(n-3))
       CALL RAPIDITY(n+2,KETA(n))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n-1,KETA(n-3))

       WPSn=WPSn1*WPSn2*WPSn3*WPSn4
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4

       RETURN                                                                   
       END 

c==============================================================
       SUBROUTINE pshkpos2_2_0(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2,IFLG3
       DOUBLE PRECISION rdn(30),rdn2(2),rdn3(2)
       DOUBLE PRECISION SS,x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4
       DOUBLE PRECISION mtemp1_1,mtemp1_2,mtemp2_1,mtemp2_2,mmx,mmx1
       double precision mmx2,mini2,max2,mmass(4),ipw12,ipwm,mmxminl,mmxmaxl

       WPSn = 0d0
       IFLG = 0 
       IFLG2 = 0 
       ipwm = 1d0

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
          if (ipwm.eq.0) THEN
             mmx1 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-3)     
             mmx2 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-4)     
c     mmx2=dlog(8d0**2)+(dlog(16d0**2)-dlog(8d0**2))*rdn(3*(n-1)-4)     
             M(3,n) = 1
             M(2,n) = mmx1
             M(1,n) = dsqrt(M(2,n))
             M(3,n-1) = 1
             M(2,n-1) = mmx2
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (( max2 -mini2 )/2d0/pi)**2
c     WPSn1=((dlog(max2)- dlog(mini2))/2.0d0/Pi)*((dlog(16d0**2)- dlog(8d0**2))/2.0d0/Pi)
             Jacobn1 = 1d0
          elseif (ipwm.eq.1) then
             mmx1 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-3)     
             mmx2 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-4)     
c     mmx2=dlog(8d0**2)+(dlog(16d0**2)-dlog(8d0**2))*rdn(3*(n-1)-4)     
             M(3,n) = 1
             M(2,n) = dexp(mmx1)
             M(1,n) = dsqrt(M(2,n))
             M(3,n-1) = 1
             M(2,n-1) = dexp(mmx2)
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (( dlog(max2) -dlog(mini2) )/2d0/pi)**2
c     WPSn1=((dlog(max2)- dlog(mini2))/2.0d0/Pi)*((dlog(16d0**2)- dlog(8d0**2))/2.0d0/Pi)
             Jacobn1 = dexp(mmx1)*dexp(mmx2)
          elseif(ipwm.ne.1) then
             mmxminl = mini2**(1d0-ipwm)
             mmxmaxl = max2**(1d0-ipwm)
             mmx1 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-3)     
             M(2,n) = (mmx1)**(1d0/(1d0-ipwm))
             M(1,n) = dsqrt(M(2,n))
             M(3,n) = 1
             mmx2 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-4)     
             M(2,n-1) = (mmx2)**(1d0/(1d0-ipwm))
             M(1,n-1) = dsqrt(M(2,n-1))
             M(3,n-1) = 1
             Jacobn1 = mmx1**(ipwm/(1d0-ipwm))/(ipwm-1d0)
     &            *mmx2**(ipwm/(1d0-ipwm))/(ipwm-1d0)
             WPSn1 = ((mmxmaxl -mmxminl)/(2*pi))**2
          endif


c          AYEcut2 = dlog(( dsqrt(SS) +dsqrt( SS -M(2,n) ) )/M(1,n))
c          AYEcut3 = dlog(( dsqrt(SS) +dsqrt( SS -M(2,n-1) ) )/M(1,n-1))

c          CALL pshk(SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          CALL pshkpos(0,0,SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
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

       KPT2(n) = P(1,n+2)**2 +P(2,n+2)**2
       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT2(n-2) = P(1,n)**2 +P(2,n)**2
       KPT2(n-3) = P(1,n-1)**2 +P(2,n-1)**2
       KPT(n) = dsqrt(KPT2(n))
       KPT(n-1) = dsqrt(KPT2(n-1))
       KPT(n-2) = dsqrt(KPT2(n-2))
       KPT(n-3) = dsqrt(KPT2(n-3))
       CALL COLAT2(n+2,KPHI(n))
       CALL COLAT2(n+1,KPHI(n-1))
       CALL COLAT2(n,KPHI(n-2))
       CALL COLAT2(n-1,KPHI(n-3))
       CALL RAPIDITY(n+2,KETA(n))
       CALL RAPIDITY(n+1,KETA(n-1))
       CALL RAPIDITY(n,KETA(n-2))
       CALL RAPIDITY(n-1,KETA(n-3))

       WPSn=WPSn1*WPSn2*WPSn3*WPSn4
       Jacobn=Jacobn1*Jacobn2*Jacobn3*Jacobn4

       RETURN                                                                   
       END 

c=========================================================
       SUBROUTINE pshkpos3(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
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
c          CALL pshk(SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          CALL pshkpos(0,0,SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
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
       SUBROUTINE pshkpos41(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
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
c       CALL pshk(SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
       CALL pshkpos(0,0,SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
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
       SUBROUTINE pshkpos42(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
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
c       CALL pshk(SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
       CALL pshkpos(0,0,SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
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
       SUBROUTINE pshkpos43(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
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

c       CALL pshk0(SS,n-2,rdn,x1,x2,WPSn3,Jacobn3,IFLG)
       CALL pshkpos(0,0,SS,n-2,rdn,x1,x2,WPSn3,Jacobn3,IFLG)

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

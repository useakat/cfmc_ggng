c----------------------------------------------------------------------------------------------
       SUBROUTINE pshkpos(ipos,delyflag,SS,n,rdn,x1,x2,WPSn,Jacob,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20),DKETA(20)
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3,AYEcut4,ptmin_cut(20),ptmax_cut(20)
       COMMON /cut2/ptmin_cut,ptmax_cut,PTcut2,AYEcut2,AYEcut3,AYEcut4
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG
       DOUBLE PRECISION rdn(30),KEPT2(20),KKPT2(20)
       DOUBLE PRECISION spt1,spt2,E00,PZ0,cutpt,cutpt2,cuteta(20)
       DOUBLE PRECISION SS,sh,x1,x2,pipt2,WPSn,pau2,Jacob
       double precision jjacob(3*n-2),wwgt(3*n-2)
       real*8 ipw1,ipw2,cutpt_min2,cutpt_max2
       real*8 ptmaxl,ptminl,rnd12,cutps,cutps_f
       integer idivps,ipos,delyflag,idebug
       real*8 M_tmp(3),KPT2_tmp,KPhi_tmp,cuteta_tmp,meff2
       real*8 Amin(20),ETmin(20),pmeff2(20)

       idebug = 0 ! debug mode (Use only for test)
       IFLG = 0 
       spt1 = 0d0 
       spt2 = 0d0 
       cutps_f = 2d0
       ipw1 = 1d0
       ipw2 = 1.5d0

       pau2 = (dsqrt(SS)/2d0)**2

       do i = 1,n
          cuteta(i) = ayecut
c$$$          if ((M(3,i+2).eq.1).or.(M(1,i+2).gt.0d0)) then
c$$$             if (M(3,i+2).eq.1) then
c$$$c                cuteta = dacosh( dsqrt( ss/( M(2,i+2) +KPT2(i) ) ) )
c$$$                cuteta(i) = ayecut
c$$$             elseif (M(1,i+2).gt.0d0) then
c$$$                cuteta(i) = min(ayecut
c$$$c     &               ,dacosh( dsqrt( ss/( M(2,i+2) +KPT2(i) ) ) ))
c$$$c     &               ,dacosh( dsqrt( ss/( M(2,i+2) ) ) ))
c$$$             endif
c$$$          else
c$$$             cuteta(i) = ayecut
c$$$          else
c$$$             iflg = 1
c$$$             return
c$$$          endif
       enddo
       do i = 1,n-1
          KETA(i) = -cuteta(i) +2*cuteta(i)*rdn(3*i-1)
          jjacob(3*i-1) = 1d0
          wwgt(3*i-1) = 2*cuteta(i)          
       enddo
       KETA(n) = -cuteta(n) +2*cuteta(n)*rdn(3*n-2)
       jjacob(3*n-2) = 1d0
       wwgt(3*n-2) = 2*cuteta(n)

c$$$       do i = 1,n
c$$$          ETmin(i) = dsqrt(M(2,i+2) +ptmin_cut(i)**2)
c$$$       enddo
c$$$       do i = 1,n-1
c$$$          Amin(i) = 0d0
c$$$          do j = 1,n-1
c$$$             if (j.ne.i) then
c$$$                Amin(i) = Amin(i) +ETmin(j)*dexp(KETA(j)-KETA(i))
c$$$             endif
c$$$          enddo
c$$$          pmeff2(i) = 2*M(2,i+2)/( 1d0 +2*M(1,i+2)/Amin(i) )
c$$$       enddo 
       
       do i = 1,n-1
          if ((M(3,i+2).eq.1).or.(M(1,i+2).gt.0d0)) then
             if (M(3,i+2).eq.1) then 
c                meff2 = pmeff2(i)
                meff2 = M(2,i+2)
c                meff2 = Qcut**2
                cutpt_min2 = meff2 +ptmin_cut(i)**2
                cutpt_max2 = min(meff2 +ptmax_cut(i)**2,ss)
c                cutpt_max2 = meff2 +ptmax_cut(i)**2
                ipw1 = 2d0  ! for the pt-jacobian of a splitting leg
                if (ipw1.eq.0) then
                   ptminl = cutpt_min2
                   ptmaxl = cutpt_max2
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = KEPT2(i) -meff2
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = 1d0
                   wwgt(3*i-2) = ptmaxl -ptminl
                elseif (ipw1.eq.1d0) then
                   ptminl = dlog(cutpt_min2)
                   ptmaxl = dlog(cutpt_max2)
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = dexp(KEPT2(i)) -meff2
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = dexp(KEPT2(i))
                   wwgt(3*i-2) = ptmaxl -ptminl
                elseif (ipw1.ne.1) then
                   ptminl = (cutpt_min2)**(1d0-ipw1)
                   ptmaxl = (cutpt_max2)**(1d0-ipw1)
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = (KEPT2(i))**(1d0/(1d0-ipw1)) -meff2
                   KPT(i) = dsqrt(KPT2(i))
                   if (ipw1.lt.1) then
                      jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(1d0-ipw1)
                      wwgt(3*i-2) = (ptmaxl -ptminl)
                   elseif (ipw1.gt.1) then
                      jjacob(3*i-2) = KEPT2(i)**(ipw1/(1d0-ipw1))/(ipw1-1d0)
                      wwgt(3*i-2) = -(ptmaxl -ptminl)
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
          else
             iflg = 1
             return
          endif

          KPhi(i) = 2*pi*rdn(3*i)
          jjacob(3*i) = 1d0
          wwgt(3*i) = 2*pi

          spt1 = spt1 -KPT(i)*dcos(KPhi(i))
          spt2 = spt2 -KPT(i)*dsin(KPhi(i))
       enddo

       KPT2(n) = spt1**2 +spt2**2
       KPT(n) = dsqrt(KPT2(n))
       if ((KPT(n).ge.ptmin_cut(n)).and.(KPT(n).le.ptmax_cut(n))) then
          continue
       else
          if (idebug.eq.1) then
             write(*,*) "IFLG=1: KPT(n) is out of range"
          endif
          IFLG = 1
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

       P(4,1) = 0d0
       P(1,1) = 0d0
       P(2,1) = 0d0
       P(3,1) = 0d0
       P(5,1) = 0d0
       P(4,2) = 0d0
       P(1,2) = 0d0
       P(2,2) = 0d0
       P(3,2) = 0d0
       P(5,2) = 0d0
       x1 = ( E00 +PZ0 )/dsqrt(SS) 
       x2 = ( E00 -PZ0 )/dsqrt(SS) 
       sh = ss*x1*x2
       If((x1.gt.1d0).or.(x2.gt.1d0)) then
          if (idebug.eq.1) then
             write(*,*) "x1 or x2 > 1, iflg = 1"
          endif
          IFLG = 1
          return
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
       REAL*8 PTcut2,AYEcut2,AYEcut3,AYEcut4,ptmin_cut(20),ptmax_cut(20)
       COMMON /cut2/ptmin_cut,ptmax_cut,PTcut2,AYEcut2,AYEcut3,AYEcut4
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2,idebug
       DOUBLE PRECISION rdn(30),rdn2(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mini2,max2
       double precision K(0:3),KK(0:3,2),mmass(2),P12(5),PT12,eta12
       real*8 phi12,ET12,a1,a2,sinphi_max,sinphi_min,maxphi,minphi
       real*8 ipw,sqrts,beta,sinth,costh,costhcut_kt
       real*8 costhcut,mincosth,maxcosth,PPP(0:3,2),PP(0:3,2)
       real*8 phi,ETij,gamma,gambeta,D,pau,mmin,mmax,m12
       real*8 esbeta,es4sq
       external esbeta,es4sq

       idebug = 0 ! debug mode (Use only for test)
       WPSn = 0d0
       wpsn1 = 1d0
       wpsn2 = 1d0
       wpsn3 = 1d0
       IFLG = 0 
       IFLG2 = 0 
       ipw = 1d0  ! for the m-jacobian of a splitting leg

       pau = dsqrt(SS)/2d0
       mmin = Qcut
       mmax = Qmax

       do J = 1,2
          rdn2(J) = rdn(3*n-1-J)
       enddo

       mtemp1 = M(1,n+1)
       mtemp2 = M(2,n+1)       
       if(n.gt.2) then
          if (ipw.eq.1) then
             mini2 = dlog(mmin**2)
             max2 = dlog(mmax**2)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n+1) = 1
             M(2,n+1) = dexp(mmx)
             M(1,n+1) = dsqrt(M(2,n+1))
             WPSn1 = ( max2 -mini2 )/(2*pi)
             Jacobn1 = dexp(mmx)
         elseif (ipw.gt.1) then
             mini2 = (mmax**2)**(1d0-ipw)
             max2 = (mmin**2)**(1d0-ipw)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n+1) = 1
             M(2,n+1) = mmx**(1d0/(1d0-ipw))
             M(1,n+1) = dsqrt(M(2,n+1)) 
             WPSn1 = -( max2 -mini2 )/(2*pi)
             Jacobn1 = mmx**(ipw/(1d0-ipw))/(ipw-1d0)
         endif
         m12 = M(1,n+1)         ! invariant mass of the particle 1 and 2
         if ((m12.ge.mmin).and.(m12.le.mmax)) then
            continue
         else
            if (idebug.eq.1) then
               write(*,*) "IFLG=1: m12 is out of range"
               write(*,*) "m12 =",m12
            endif
            IFLG = 1
            return
         endif

         if (m12.lt.2*ptcut) then
            ptmin_cut(n-1) = 2*ptcut*dsqrt(1d0 -( m12/(2*ptcut) )**2) ! pT12 minimum cut 
            ptmax_cut(n-1) = min( ptcut*( m12**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
         else
            ptmin_cut(n-1) = 0d0 ! pT12 minimum cut 
            ptmax_cut(n-1) = pau ! pT12 maximum cut 
         endif
         AYEcut2 = dlog((dsqrt(SS) +dsqrt(SS -m12**2))/m12)

         CALL pshkpos(0,0,SS,n-1,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
         IF(IFLG.ne.0) return

       else
          if (idebug.eq.1) then
             write(99,*) "ERROR: n<2 (pshk2)"
          endif
          stop
       endif

       P12(1) = P(1,n+1)
       P12(2) = P(2,n+1)
       P12(3) = P(3,n+1)
       P12(4) = P(4,n+1)
       P12(5) = P(5,n+1)
       pt12 = KPT(n-1)

       M(1,n+1) = mtemp1
       M(2,n+1) = mtemp2
       mmass(1) = M(1,n+1)
       mmass(2) = M(1,n+2)
c       call ps2bd_in_pshk(rdn2,P12,mmass,P(1,n+1),wpsn3,jacobn3,iflg2)
       call ps2bd2_in_pshk(rdn2,P12,mmass,ptcut,rcut,P(1,n+1),wpsn3
     &      ,jacobn3,iflg2)
       IF(IFLG2.ne.0) then
          IFLG=1
          RETURN
       endif

       KPT2(n) = P(1,n+2)**2 +P(2,n+2)**2
       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT(n) = dsqrt(KPT2(n))
       KPT(n-1) = dsqrt(KPT2(n-1))

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
       REAL*8 PTcut2,AYEcut2,AYEcut3,AYEcut4,ptmin_cut(20),ptmax_cut(20)
       COMMON /cut2/ptmin_cut,ptmax_cut,PTcut2,AYEcut2,AYEcut3,AYEcut4
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2,IFLG3,idebug
       DOUBLE PRECISION rdn(30),rdn2(2),rdn3(2)
       DOUBLE PRECISION SS,x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4
       DOUBLE PRECISION mtemp1_1,mtemp1_2,mtemp2_1,mtemp2_2,mmx,mmx1
       double precision mmx2,mini2,max2,mmass(4),ipw12,ipwm,mmxminl,mmxmaxl
       real*8 m12_1,m12_2,mmin,mmax

       WPSn = 0d0
       IFLG = 0 
       IFLG2 = 0 
       ipwm = 1d0

       pau = dsqrt(SS)/2d0
       mmin = Qcut
       mmax = Qmax

       do J = 1,2
          rdn2(J) = rdn(3*n-1-J)
          rdn3(J) = rdn(3*n-3-J)
       enddo

       mtemp1_1 = M(1,n)
       mtemp1_2 = M(2,n)                
       mtemp2_1 = M(1,n-1)
       mtemp2_2 = M(2,n-1)                
       mini2 = mmin**2
       max2 = mmax**2

       if(n.ge.4) then
          if (ipwm.eq.0) THEN
             mmx1 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-3)     
             mmx2 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-4)     
             M(3,n) = 1
             M(2,n) = mmx1
             M(1,n) = dsqrt(M(2,n))
             M(3,n-1) = 1
             M(2,n-1) = mmx2
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (( max2 -mini2 )/(2*pi))**2
             Jacobn1 = 1d0
          elseif (ipwm.eq.1) then
             mmx1 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-3)     
             mmx2 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-4)     
             M(3,n) = 1
             M(2,n) = dexp(mmx1)
             M(1,n) = dsqrt(M(2,n))
             M(3,n-1) = 1
             M(2,n-1) = dexp(mmx2)
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (( dlog(max2) -dlog(mini2) )/(2*pi))**2
             Jacobn1 = dexp(mmx1)*dexp(mmx2)
          elseif(ipwm.ne.1) then
             mmxminl = mini2**(1d0-ipwm)
             mmxmaxl = max2**(1d0-ipwm)
             mmx1 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-3)     
             M(3,n) = 1
             M(2,n) = (mmx1)**(1d0/(1d0-ipwm))
             M(1,n) = dsqrt(M(2,n))
             mmx2 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-4)     
             M(3,n-1) = 1
             M(2,n-1) = mmx2**(1d0/(1d0-ipwm))
             M(1,n-1) = dsqrt(M(2,n-1))
             WPSn1 = (-(mmxmaxl -mmxminl)/(2*pi))**2
             Jacobn1 = mmx1**(ipwm/(1d0-ipwm))/(ipwm-1d0)
     &            *mmx2**(ipwm/(1d0-ipwm))/(ipwm-1d0)
          endif
         m12_1 = M(1,n)         ! invariant mass of the particle 1 and 2
         m12_2 = M(1,n-1)         ! invariant mass of the particle 1 and 2
         if ((m12_1.ge.mmin).and.(m12_1.le.mmax)) then
            continue
         else
            if (idebug.eq.1) then
               write(*,*) "IFLG=1: m12_1 is out of range"
               write(*,*) "m12_1 =",m12_1
            endif
            IFLG = 1
            return
         endif
         if ((m12_2.ge.mmin).and.(m12_2.le.mmax)) then
            continue
         else
            if (idebug.eq.1) then
               write(*,*) "IFLG=1: m12_2 is out of range"
               write(*,*) "m12_2 =",m12_2
            endif
            IFLG = 1
            return
         endif

         if (m12_1.lt.2*ptcut) then
            ptmin_cut(n-2) = 2*ptcut*dsqrt(1d0 -( m12_1/(2*ptcut) )**2) ! pT12 minimum cut 
            ptmax_cut(n-2) = min( ptcut*( m12_1**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
         else
            ptmin_cut(n-2) = 0d0 ! pT12 minimum cut 
            ptmax_cut(n-2) = pau ! pT12 maximum cut 
         endif
         if (m12_2.lt.2*ptcut) then
            ptmin_cut(n-3) = 2*ptcut*dsqrt(1d0 -( m12_2/(2*ptcut) )**2) ! pT12 minimum cut 
            ptmax_cut(n-3) = min( ptcut*( m12_2**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
         else
            ptmin_cut(n-3) = 0d0 ! pT12 minimum cut 
            ptmax_cut(n-3) = pau ! pT12 maximum cut 
         endif
         AYEcut2 = dlog((dsqrt(SS) +dsqrt(SS -m12_1**2))/m12_1)
         AYEcut3 = dlog((dsqrt(SS) +dsqrt(SS -m12_2**2))/m12_2)

          CALL pshkpos(0,0,SS,n-2,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          IF(IFLG.ne.0) then
             return
          endif
       else
          if (idebug.eq.1) then
             print *,"ERROR, n < 4 in pshk2_2 subroutine"
          endif
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

       mmass(1) = M(1,n+1)
       mmass(2) = M(1,n+2)
       call ps2bd2_in_pshk(rdn2,P(1,50),mmass,ptcut,rcut,P(1,n+1),wpsn3
     &      ,jacobn3,iflg2)
       IF(IFLG2.ne.0) then
          IFLG = 1
          RETURN
       endif

       mmass(3) = M(1,n-1)
       mmass(4) = M(1,n)
       call ps2bd2_in_pshk(rdn3,P(1,49),mmass(3),ptcut,rcut,P(1,n-1)
     &      ,wpsn4,jacobn4,iflg3)
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
c       CALL COLAT2(n+2,KPHI(n))
c       CALL COLAT2(n+1,KPHI(n-1))
c       CALL COLAT2(n,KPHI(n-2))
c       CALL COLAT2(n-1,KPHI(n-3))
c       CALL RAPIDITY(n+2,KETA(n))
c       CALL RAPIDITY(n+1,KETA(n-1))
c       CALL RAPIDITY(n,KETA(n-2))
c       CALL RAPIDITY(n-1,KETA(n-3))

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
       REAL*8 PTcut2,AYEcut2,AYEcut3,AYEcut4,ptmin_cut(20),ptmax_cut(20)
       COMMON /cut2/ptmin_cut,ptmax_cut,PTcut2,AYEcut2,AYEcut3,AYEcut4
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
       SUBROUTINE pshkpos2_3(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3,AYEcut4,ptmin_cut(20),ptmax_cut(20)
       COMMON /cut2/ptmin_cut,ptmax_cut,PTcut2,AYEcut2,AYEcut3,AYEcut4
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2,IFLG3,IFLG4,idebug
       DOUBLE PRECISION rdn(30),rdn2(2),rdn3(2),rdn4(2)
       DOUBLE PRECISION SS,x1,x2,WPSn,WPSn1,WPSn2,WPSn3,WPSn4,WPSn5
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3,Jacobn4,Jacobn5
       DOUBLE PRECISION mtemp1_1,mtemp1_2,mtemp2_1,mtemp2_2
       real*8 mtemp3_1,mtemp3_2,mmx,mmx1,mmx2,mmx3,mini2,max2
       double precision mmass(6),ipw12,ipwm,mmxminl,mmxmaxl
       real*8 m12_1,m12_2,mmin,mmax

       idebug = 0               ! debug mode (Use only for test)       
       WPSn = 0d0
       IFLG = 0 
       IFLG2 = 0 
       IFLG3 = 0 
       IFLG4 = 0 
       ipwm = 1d0
       
       pau = dsqrt(SS)/2d0
       mmin = Qcut
       mmax = Qmax

       do J = 1,2
          rdn2(J) = rdn(3*n-1-J)
          rdn3(J) = rdn(3*n-3-J)
          rdn4(J) = rdn(3*n-5-J)
       enddo
       
       mtemp1_1 = M(1,n-1)
       mtemp1_2 = M(2,n-1)                
       mtemp2_1 = M(1,n-2)
       mtemp2_2 = M(2,n-2)                
       mtemp3_1 = M(1,n-3)
       mtemp3_2 = M(2,n-3)                
       mini2 = mmin**2
       max2 = mmax**2
       
       if(n.ge.6) then
          if (ipwm.eq.0) THEN
             mmx1 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-3)     
             mmx2 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-4)     
             mmx3 = mini2 +( max2 -mini2 )*rdn(3*(n-1)-5)
             M(3,n-1) = 1
             M(2,n-1) = mmx1
             M(1,n-1) = dsqrt(M(2,n-1))
             M(3,n-2) = 1
             M(2,n-2) = mmx2
             M(1,n-2) = dsqrt(M(2,n-2))
             M(3,n-3) = 1
             M(2,n-3) = mmx3
             M(1,n-3) = dsqrt(M(2,n-3))
             WPSn1 = (( max2 -mini2 )/(2*pi))**3
             Jacobn1 = 1d0
          elseif (ipwm.eq.1) then
             mmx1 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-3)     
             mmx2 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-4)     
             mmx3 = dlog(mini2) +( dlog(max2) -dlog(mini2) )*rdn(3*(n-1)-5)
             M(3,n-1) = 1
             M(2,n-1) = dexp(mmx1)
             M(1,n-1) = dsqrt(M(2,n-1))
             M(3,n-2) = 1
             M(2,n-2) = dexp(mmx2)
             M(1,n-2) = dsqrt(M(2,n-2))
             M(3,n-3) = 1
             M(2,n-3) = dexp(mmx3)
             M(1,n-3) = dsqrt(M(2,n-3))
             WPSn1 = (( dlog(max2) -dlog(mini2) )/(2*pi))**3
             Jacobn1 = dexp(mmx1)*dexp(mmx2)*dexp(mmx3)
          elseif(ipwm.ne.1) then
             mmxminl = mini2**(1d0-ipwm)
             mmxmaxl = max2**(1d0-ipwm)
             mmx1 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-3)     
             M(3,n-1) = 1
             M(2,n-1) = mmx1**(1d0/(1d0-ipwm))
             M(1,n-1) = dsqrt(M(2,n))
             mmx2 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-4)     
             M(3,n-2) = 1
             M(2,n-2) = mmx2**(1d0/(1d0-ipwm))
             M(1,n-2) = dsqrt(M(2,n-1))
             mmx3 = mmxminl +( mmxmaxl -mmxminl )*rdn(3*(n-1)-5)     
             M(3,n-3) = 1
             M(2,n-3) = mmx3**(1d0/(1d0-ipwm))
             M(1,n-3) = dsqrt(M(2,n-2))
             WPSn1 = (-(mmxmaxl -mmxminl)/(2*pi))**3
             Jacobn1 = mmx1**(ipwm/(1d0-ipwm))/(ipwm-1d0)
     &            *mmx2**(ipwm/(1d0-ipwm))/(ipwm-1d0)
     &            *mmx3**(ipwm/(1d0-ipwm))/(ipwm-1d0)
          endif
          m12_1 = M(1,n-1)        ! invariant mass of the particle 1 and 2
          m12_2 = M(1,n-2)      ! invariant mass of the particle 1 and 2
          m12_3 = M(1,n-3)      ! invariant mass of the particle 1 and 2
          if ((m12_1.ge.mmin).and.(m12_1.le.mmax)) then
             continue
          else
             if (idebug.eq.1) then
                write(*,*) "IFLG=1: m12_1 is out of range"
                write(*,*) "m12_1 =",m12_1
             endif
             IFLG = 1
             return
          endif
          if ((m12_2.ge.mmin).and.(m12_2.le.mmax)) then
             continue
          else
             if (idebug.eq.1) then
                write(*,*) "IFLG=1: m12_2 is out of range"
                write(*,*) "m12_2 =",m12_2
             endif
             IFLG = 1
             return
          endif
          if ((m12_3.ge.mmin).and.(m12_3.le.mmax)) then
             continue
          else
             if (idebug.eq.1) then
                write(*,*) "IFLG=1: m12_2 is out of range"
                write(*,*) "m12_2 =",m12_2
             endif
             IFLG = 1
             return
          endif

          if (m12_1.lt.2*ptcut) then
             ptmin_cut(n-3) = 2*ptcut*dsqrt(1d0 -( m12_1/(2*ptcut) )**2) ! pT12 minimum cut 
             ptmax_cut(n-3) = min( ptcut*( m12_1**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
          else
             ptmin_cut(n-3) = 0d0 ! pT12 minimum cut 
             ptmax_cut(n-3) = pau ! pT12 maximum cut 
          endif
          if (m12_2.lt.2*ptcut) then
             ptmin_cut(n-4) = 2*ptcut*dsqrt(1d0 -( m12_2/(2*ptcut) )**2) ! pT12 minimum cut 
             ptmax_cut(n-4) = min( ptcut*( m12_2**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
          else
             ptmin_cut(n-4) = 0d0 ! pT12 minimum cut 
             ptmax_cut(n-4) = pau ! pT12 maximum cut 
          endif
          if (m12_3.lt.2*ptcut) then
             ptmin_cut(n-5) = 2*ptcut*dsqrt(1d0 -( m12_3/(2*ptcut) )**2) ! pT12 minimum cut 
             ptmax_cut(n-5) = min( ptcut*( m12_3**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
          else
             ptmin_cut(n-5) = 0d0 ! pT12 minimum cut 
             ptmax_cut(n-5) = pau ! pT12 maximum cut 
          endif
          AYEcut2 = dlog((dsqrt(SS) +dsqrt(SS -m12_1**2))/m12_1)
          AYEcut3 = dlog((dsqrt(SS) +dsqrt(SS -m12_2**2))/m12_2)
          AYEcut4 = dlog((dsqrt(SS) +dsqrt(SS -m12_3**2))/m12_3)

          CALL pshkpos(0,0,SS,n-3,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
          IF(IFLG.ne.0) then
             return
          endif
       else
c     print *,"ERROR, n < 4 in pshk2_2 subroutine"
          stop
       endif
       
       M(1,50) = M(1,n-1)
       M(2,50) = M(2,n-1)
       P(1,50) = P(1,n-1)
       P(2,50) = P(2,n-1)
       P(3,50) = P(3,n-1)
       P(4,50) = P(4,n-1)
       P(5,50) = P(5,n-1)
       M(1,n-1) = mtemp1_1
       M(2,n-1) = mtemp1_2
       
       M(1,49) = M(1,n-2)
       M(2,49) = M(2,n-2)
       P(1,49) = P(1,n-2)
       P(2,49) = P(2,n-2)
       P(3,49) = P(3,n-2)
       P(4,49) = P(4,n-2)
       P(5,49) = P(5,n-2)
       M(1,n-2) = mtemp2_1
       M(2,n-2) = mtemp2_2
       
       M(1,48) = M(1,n-3)
       M(2,48) = M(2,n-3)
       P(1,48) = P(1,n-3)
       P(2,48) = P(2,n-3)
       P(3,48) = P(3,n-3)
       P(4,48) = P(4,n-3)
       P(5,48) = P(5,n-3)
       M(1,n-3) = mtemp3_1
       M(2,n-3) = mtemp3_2
       
       mmass(1) = M(1,n+1)
       mmass(2) = M(1,n+2)
       call ps2bd2_in_pshk(rdn2,P(1,50),mmass,ptcut,rcut,P(1,n+1),wpsn3
     &      ,jacobn3,iflg2)
       IF(IFLG2.ne.0) then
          IFLG = 1
          RETURN
       endif
       
       mmass(3) = M(1,n-1)
       mmass(4) = M(1,n)
       call ps2bd2_in_pshk(rdn3,P(1,49),mmass(3),ptcut,rcut,P(1,n-1)
     &      ,wpsn4,jacobn4,iflg3)
       IF(IFLG3.ne.0) then
          IFLG = 1
          RETURN
       endif
       
       mmass(5) = M(1,n-3)
       mmass(6) = M(1,n-2)
       call ps2bd2_in_pshk(rdn4,P(1,48),mmass(5),ptcut,rcut,P(1,n-3)
     &      ,wpsn5,jacobn5,iflg4)
       IF(IFLG4.ne.0) then
          IFLG = 1
          RETURN
       endif
       
       KPT2(n) = P(1,n+2)**2 +P(2,n+2)**2
       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT2(n-2) = P(1,n)**2 +P(2,n)**2
       KPT2(n-3) = P(1,n-1)**2 +P(2,n-1)**2
       KPT2(n-4) = P(1,n-2)**2 +P(2,n-2)**2
       KPT2(n-5) = P(1,n-3)**2 +P(2,n-3)**2
       KPT(n) = dsqrt(KPT2(n))
       KPT(n-1) = dsqrt(KPT2(n-1))
       KPT(n-2) = dsqrt(KPT2(n-2))
       KPT(n-3) = dsqrt(KPT2(n-3))
       KPT(n-4) = dsqrt(KPT2(n-4))
       KPT(n-5) = dsqrt(KPT2(n-5))
       
       WPSn = WPSn1*WPSn2*WPSn3*WPSn4*WPSn5
       Jacobn = Jacobn1*Jacobn2*Jacobn3*Jacobn4*Jacobn5

       RETURN                                                                   
       END 
      

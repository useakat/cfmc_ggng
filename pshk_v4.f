c----------------------------------------------------------------------------------------------
       SUBROUTINE pshkpos(ipos,delyflag,SS,n,rdn,x1,x2,WPSn,Jacob,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20),DKETA(20)
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG
       DOUBLE PRECISION rdn(30),KEPT2(20),KKPT2(20)
       DOUBLE PRECISION spt1,spt2,E00,PZ0,cutpt,cutpt2,cuteta(20)
       DOUBLE PRECISION SS,sh,x1,x2,pipt2,WPSn,pau2,Jacob
       double precision jjacob(3*n-2),wwgt(3*n-2)
       real*8 ipw1,ipw2,cutpt_min2,cutpt_max2
       real*8 ptmaxl,ptminl,rnd12,cutps,cutps_f
       integer idivps,ipos,delyflag
       real*8 M_tmp(3),KPT2_tmp,KPhi_tmp,cuteta_tmp

       IFLG = 0 
       spt1 = 0d0 
       spt2 = 0d0 
       cutps_f = 2d0
       ipw1 = 1d0
       ipw2 = 1.5d0

       pau2 = (dsqrt(SS)/2d0)**2

c       if (ipos.gt.0) then ! Exchange the n-th and (n-ipos)-th momentum.
c          do i = 1,3
c             M_tmp(i) = M(i,2+n-ipos)
c             M(i,2+n-ipos) = M(i,2+n)
c             M(i,2+n) = M_tmp(i)
c          enddo
c       endif
       do i = 1,n-1
          if ((M(3,i+2).eq.1).or.(M(1,i+2).gt.0d0)) then
             if (M(3,i+2).eq.1) then 
                cutpt_min2 = M(2,i+2) +ptmin_cut**2
                cutpt_max2 = min(M(2,i+2)+ptmax_cut**2,pau2)
                ipw1 = 1d0
                if (ipw1.eq.0) then
                   ptminl = cutpt_min2
                   ptmaxl = cutpt_max2
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = KEPT2(i) -M(2,i+2)
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = 1d0
                   wwgt(3*i-2) = ptmaxl -ptminl
                elseif (ipw1.eq.1) then
                   ptminl = dlog(cutpt_min2)
                   ptmaxl = dlog(cutpt_max2)
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = dexp(KEPT2(i)) -M(2,i+2)
                   KPT(i) = dsqrt(KPT2(i))
                   jjacob(3*i-2) = dexp(KEPT2(i))
                   wwgt(3*i-2) = ptmaxl -ptminl
                elseif (ipw1.ne.1) then
                   ptminl = (cutpt_min2)**(1d0-ipw1)
                   ptmaxl = (cutpt_max2)**(1d0-ipw1)
                   KEPT2(i) = ptminl +( ptmaxl -ptminl )*rdn(3*i-2)
                   KPT2(i) = (KEPT2(i))**(1d0/(1d0-ipw1)) -M(2,i+2)
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
       if ((KPT(n).ge.ptmin_cut).and.(KPT(n).le.ptmax_cut)) then
          continue
       else
          write(*,*) "IFLG=1: KPT(n) is out of range"
          IFLG = 1
          return
       endif

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

c       if (ipos.gt.0) then  ! The (n-ipos)-th momentum is generated last.
c          do i = 1,3
c             M_tmp(i) = M(i,2+n-ipos)
c             M(i,2+n-ipos) = M(i,2+n)
c             M(i,2+n) = M_tmp(i)
c          enddo
c          KPT2_tmp = KPT2(n-ipos)
c          KPT2(n-ipos) = KPT2(n)
c          KPT2(n) = KPT2_tmp
c          KPT(n-ipos) = dsqrt(KPT2(n-ipos))
c          KPT(n) = dsqrt(KPT2(n))
c          KPhi_tmp = KPhi(n-ipos)
c          KPhi(n-ipos) = KPhi(n)
c          KPhi(n) = KPhi_tmp
c          cuteta_tmp = cuteta(n-ipos)
c          cuteta(n-ipos) = cuteta(n)
c          cuteta(n) = cuteta_tmp
c       endif

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
          write(*,*) "x1 or x2 > 1, iflg = 1"
          write(*,*)
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
c       include 'inc/const.inc'
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG,IFLG2
       DOUBLE PRECISION rdn(30),rdn2(2)
       DOUBLE PRECISION SS, x1,x2,WPSn,WPSn1,WPSn2,WPSn3
       DOUBLE PRECISION Jacobn,Jacobn1,Jacobn2,Jacobn3
       DOUBLE PRECISION mtemp1,mtemp2,mmx,mini2,max2
       double precision K(0:3),KK(0:3,2),mmass(2),P12(0:3),PT12,eta12
       real*8 phi12,ET12,a1,a2,sinphi_max,sinphi_min,maxphi,minphi
       real*8 ipw,sqrts,beta,sinth,costh,costhcut_kt
       real*8 costhcut,mincosth,maxcosth,PPP(0:3,2),PP(0:3,2)
       real*8 phi,ETij,gamma,gambeta,D,pau,mmin,mmax,m12
       real*8 esbeta,es4sq
       external esbeta,es4sq

       WPSn = 0d0
       wpsn1 = 1d0
       wpsn2 = 1d0
       wpsn3 = 1d0
       IFLG = 0 
       IFLG2 = 0 
       ipw = 1d0

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
             WPSn1 = ( max2 -mini2 )/2d0/Pi
             Jacobn1 = dexp(mmx)
         elseif (ipw.gt.1) then
             mini2 = (mmax**2)**(1d0-ipw)
             max2 = (mmin**2)**(1d0-ipw)
             mmx = mini2 +( max2 -mini2 )*rdn(3*(n-1)-1)     
             M(3,n+1) = 1
             M(2,n+1) = mmx**(1d0/(1d0-ipw))
             M(1,n+1) = dsqrt(M(2,n+1)) 
             WPSn1 = -( max2 -mini2 )/2d0/Pi
             Jacobn1 = mmx**(ipw/(1d0-ipw))/(ipw-1d0)
         endif
         m12 = M(1,n+1)         ! invariant mass of the particle 1 and 2
         if ((m12.ge.mmin).and.(m12.le.mmax)) then
            continue
         else
            write(*,*) "IFLG=1: m12 is out of range"
            write(*,*) "m12 =",m12
            IFLG = 1
            return
         endif

         if (m12.lt.2*ptcut) then
            ptmin_cut = 2*ptcut*dsqrt(1d0 -( m12/(2*ptcut) )**2) ! pT12 minimum cut 
            ptmax_cut = min( ptcut*( m12**2/mmin**2 +1d0 ), pau ) ! pT12 maximum cut
         else
            ptmin_cut = 0d0 ! pT12 minimum cut 
            ptmax_cut = pau ! pT12 maximum cut 
         endif
         AYEcut2 = dlog((dsqrt(SS) +dsqrt(SS -m12**2))/m12)

         CALL pshkpos(0,0,SS,n-1,rdn,x1,x2,WPSn2,Jacobn2,IFLG)
         IF(IFLG.ne.0) return

       else
          write(99,*) "ERROR: n<2 (pshk2)"
          stop
       endif

c$$$       write(*,*) P(4,1),P(1,1),P(2,1),P(3,1)
c$$$       write(*,*) P(4,2),P(1,2),P(2,2),P(3,2)
c$$$       write(*,*) P(4,3),P(1,3),P(2,3),P(3,3)
c$$$       write(*,*) P(4,4),P(1,4),P(2,4),P(3,4)
c$$$       write(*,*) P(4,3)**2 -P(1,3)**2 -P(2,3)**2 -P(3,3)**2
c$$$       write(*,*) P(4,4)**2 -P(1,4)**2 -P(2,4)**2 -P(3,4)**2
c$$$       write(*,*) M(2,4)

       P12(0) = P(4,n+1)
       P12(1) = P(1,n+1)
       P12(2) = P(2,n+1)
       P12(3) = P(3,n+1)
       pt12 = KPT(n-1)

       M(1,n+1) = mtemp1
       M(2,n+1) = mtemp2
       mmass(1) = M(1,n+1)
       mmass(2) = M(1,n+2)

       ET12 = dsqrt(m12**2 +pt12**2)
       gamma = ET12/m12
       gambeta = pt12/m12

       costh = -1d0 +2*rdn2(1)
       wpsn3 = 2d0
       jacob3 = 1d0

       phi = 2*pi*rdn2(2)
       wpsn3 = 2*pi
       jacob3 = 1d0

c     generation of costh*
c
c$$$       if (m12.lt.2*ptcut) then
c$$$          maxcosth = 1d0/pt12*(ET12 -2*ptcut)
c$$$          if (rcut.le.2*m12/pt12) then
c$$$             mincosth = 0d0
c$$$          else
c$$$             mincosth = pt12/ET12*dsqrt(1d0 -(2*m12/(pt12*rcut))**2)
c$$$          endif
c$$$          
c$$$          if (rdn2(1).lt.0.5) then
c$$$             costh = -maxcosth +(-mincosth +maxcosth)*2*rdn2(1)
c$$$          else
c$$$             costh = mincosth +(maxcosth -mincosth)*2*(rdn2(1)-0.5d0)
c$$$          endif
c$$$          wpsn3 = (maxcosth -mincosth)*2
c$$$          jacobn3 = 1d0
c$$$       else
c$$$          maxcosth = 1d0
c$$$          mincosth = -1d0
c$$$          costh = mincosth +(maxcosth -mincosth)*rdn2(1)
c$$$          wpsn3 = (maxcosth -mincosth)
c$$$          jacobn3 = 1d0
c$$$       endif
       if ((costh.ge.-1d0).and.(costh.le.1d0)) then
          continue
       else
          write(*,*) "IFLG=1: costh is out of range"
          write(*,*) "costh =",costh
          iflg = 1
          return
       endif
       sinth = dsqrt(1d0 -costh**2)

c$$$c     generation of phi*
c$$$c
c$$$c       write(*,*) "var",ptcut,m12,pt12,costh,sinth,gamma,gambeta/gamma
c$$$       if (m12.lt.2*ptcut) then
c$$$c      checking if the argument of dsqrt is positive
c$$$          a1 = 2*dsqrt(4*m12**4 +costh**2*pt12**2*ET12**2*rcut**4)
c$$$          a2 = rcut**2*(pt12**2 +costh**2*ET12**2)
c$$$          if (a1.lt.a2) then
c$$$             write(*,*) "IFLG=1: a1 < a2"
c$$$             write(*,*) a1,a2
c$$$             iflg = 1
c$$$             return
c$$$          endif
c$$$          
c$$$          sinphi_max = 1d0/(sinth*m12*rcut)*dsqrt(a1-a2)
c$$$          if (sinphi_max.le.1d0) then
c$$$             maxphi = dasin(sinphi_max)
c$$$          else
c$$$             maxphi = pi/2d0
c$$$          endif
c$$$ 
c$$$          if (costh.le.0d0) then
c$$$             if (ptcut.ge.0.5*dabs(pt12 +ET12*costh)) then
c$$$                sinphi_min = 1d0/sinth*dsqrt( (2*ptcut/m12)**2 
c$$$     &               -gamma**2*(gambeta/gamma +costh)**2 )
c$$$                minphi = dasin(sinphi_min)
c$$$             else
c$$$                minphi = 0d0
c$$$             endif
c$$$          else
c$$$             if (ptcut.ge.0.5d0*dabs(pt12 -ET12*costh)) then
c$$$                sinphi_min = 1d0/sinth*dsqrt( (2*ptcut/m12)**2 
c$$$     &               -gamma**2*(gambeta/gamma -costh)**2 )
c$$$                minphi = dasin(sinphi_min)
c$$$             else
c$$$                minphi = 0d0
c$$$             endif
c$$$          endif
c$$$
c$$$          if (rdn2(2).lt.0.25d0) then
c$$$             phi = minphi +4*rdn2(2)*(maxphi-minphi)
c$$$          elseif (rdn2(2).lt.0.5d0) then
c$$$             phi = pi -maxphi +4*(rdn2(2)-0.25d0)*(maxphi-minphi)
c$$$          elseif (rdn2(2).lt.0.75d0) then
c$$$             phi = pi +minphi +4*(rdn2(2)-0.5d0)*(maxphi-minphi)
c$$$          elseif (rdn2(2).le.1d0) then
c$$$             phi = 2*pi -maxphi +4*(rdn2(2)-0.75d0)*(maxphi-minphi)
c$$$          endif
c$$$          wpsn3 = wpsn3*4*(maxphi-minphi)
c$$$          jacobn3 = jacobn3*1d0
c$$$       else
c$$$          maxphi = 2*pi
c$$$          minphi = 0d0
c$$$          phi = minphi +(maxphi -minphi)*rdn2(2)
c$$$          wpsn3 = wpsn3*(maxphi-minphi)
c$$$          jacobn3 = jacobn3*1d0
c$$$       endif
       if ((phi.ge.0d0).and.(phi.le.2*pi)) then
          continue
       else
          write(*,*) "IFLG=1: phi is out of range"
          write(*,*) "phi =",phi
          iflg = 1
          return
       endif

c     construct 4-momentum of the particle 1 and 2 in the center of mass frame of the particle 12
c
       beta = esbeta(mmass(1)**2/m12**2,mmass(2)**2/m12**2)
c       call mom2cx(m12,mmass(1),mmass(2),costh,phi, PPP(0,1),PPP(0,2))
       PPP(0,1) = m12/2d0*(1d0 +(mmass(1)**2 -mmass(2)**2)/m12**2)
       PPP(1,1) = m12/2d0*beta*costh
       PPP(2,1) = m12/2d0*beta*sinth*dsin(phi)
       PPP(3,1) = -m12/2d0*beta*sinth*dcos(phi)
       PPP(0,2) = m12/2d0*(1d0 +(mmass(2)**2 -mmass(1)**2)/m12**2)
       PPP(1,2) = -PPP(1,1)
       PPP(2,2) = -PPP(2,1)
       PPP(3,2) = -PPP(3,1)

c       write(*,*) P12(0),P12(1),P12(2),P12(3)
c       write(*,*) dsqrt(P12(0)**2 -P12(1)**2 -P12(2)**2 -P12(3)**2)
c       write(*,*) m12
c       write(*,*) PPP(0,1)+PPP(0,2)
c       write(*,*) PPP(1,1)+PPP(1,2)
c       write(*,*) PPP(2,1)+PPP(2,2)
c       write(*,*) PPP(3,1)+PPP(3,2)
c       write(*,*) PPP(0,1)**2 -PPP(1,1)**2 -PPP(2,1)**2 -PPP(3,1)**2
c       write(*,*) PPP(0,2)**2 -PPP(1,2)**2 -PPP(2,2)**2 -PPP(3,2)**2

       call boostx(PPP(0,1),P12,PP(0,1))
       call boostx(PPP(0,2),P12,PP(0,2))
       
c$$$       write(*,*) P12(0),P12(1),P12(2),P12(3)
c$$$       write(*,*) dsqrt(P12(0)**2 -P12(1)**2 -P12(2)**2 -P12(3)**2)
c$$$       write(*,*) m12
c$$$       write(*,*) PP(0,1)+PP(0,2)
c$$$       write(*,*) PP(1,1)+PP(1,2)
c$$$       write(*,*) PP(2,1)+PP(2,2)
c$$$       write(*,*) PP(3,1)+PP(3,2)
c$$$       write(*,*) PP(0,1)**2 -PP(1,1)**2 -PP(2,1)**2 -PP(3,1)**2
c$$$       write(*,*) PP(0,2)**2 -PP(1,2)**2 -PP(2,2)**2 -PP(3,2)**2

       P(1,n+1) = PP(1,1)
       P(2,n+1) = PP(2,1)
       P(3,n+1) = PP(3,1)
       P(4,n+1) = PP(0,1)
       P(5,n+1) = dsqrt(PP(1,1)**2 +PP(2,1)**2 +PP(3,1)**2)
       P(1,n+2) = PP(1,2)
       P(2,n+2) = PP(2,2)
       P(3,n+2) = PP(3,2)
       P(4,n+2) = PP(0,2)
       P(5,n+2) = dsqrt(PP(1,2)**2 +PP(2,2)**2 +PP(3,2)**2)

       KPT2(n-1) = P(1,n+1)**2 +P(2,n+1)**2
       KPT2(n) = P(1,n+2)**2 +P(2,n+2)**2
       KPT(n-1) = dsqrt(KPT2(n-1))
       KPT(n) = dsqrt(KPT2(n))

       WPSn = WPSn1*WPSn2*WPSn3*beta/(8*pi)/(4*pi)
       Jacobn = Jacobn1*Jacobn2*Jacobn3

       RETURN                                                                   
       END 

c==============================================================
       SUBROUTINE pshkpos2_2(SS,n,rdn,x1,x2,WPSn,Jacobn,IFLG)
       DOUBLE PRECISION P(5,50),M(3,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,rcut,QCut,Pi,Qmax
       COMMON /cut/PTcut,AYEcut,rcut,QCut,Pi,Qmax
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
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
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
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
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
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
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
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
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
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
       REAL*8 PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
       COMMON /cut2/PTcut2,AYEcut2,AYEcut3,ptmin_cut,ptmax_cut
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

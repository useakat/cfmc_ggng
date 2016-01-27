      subroutine ptcut_chk(n,p,ptcut, flg)
      implicit none

      integer i,n
      real*8 p(0:3,n),ptcut
      integer flg

      real*8 pt
      external pt
     
      flg = 0
     
      do i = 1,n
         if( pt(p(0,i)).lt.ptcut ) flg = 1
      enddo

      end


      subroutine etacut_chk(n,p,etacut, flg)
      implicit none

      integer i,n
      real*8 p(0:3,n),etacut
      integer flg

      real*8 eta
      external eta
     
      flg = 0
     
      do i = 1,n
         if( abs(eta(p(0,i))).gt.etacut ) flg = 1
      enddo

      end


      subroutine ptijcut_chk(n,p,ptijcut, flg)
      implicit none

      integer i,j,n
      real*8 p(0:3,n),ptijcut
      integer flg

      real*8 pTij
      external pTij
     
      flg = 0
     
      do i = 1,n-1
         do j = i+1,n
            if (pTij(p(0,i),p(0,j)).lt.ptijcut) then
               flg = 1
               return
            endif
         enddo
      enddo

c      do i = 1,n-1
c         pt(i) = ppt(p(0,i))
c         do j = i+1,n
c            pt(j) = ppt(p(0,j))
c            dr = pdeltar(p(0,i),p(0,j))
c            ptij = min(pt(i),pt(j))*dr
c            if (ptij.lt.ptijcut) then
c               flg = 1
c               return
c            endif
c         enddo
c      enddo
      return
      end


      subroutine mijcut_chk(n,p,mijcut, flg)
c by Yoshitaro Takaesu -Aug/04/2011 @KEK
      implicit none
      
      integer i,j
      
      integer n
      real*8 p(0:3,n),mijcut
      integer flg

      real*8 pmij

      real*8 mij
      external mij

      flg = 0
      
      do i = 1,n-1
         do j = i+1,n
            pmij = mij(p(0,i),p(0,j))     
            if( pmij.lt.mijcut ) then
               flg = 1
               return
            endif
         enddo
       enddo

       end

      subroutine drcut_chk(n,p,drcut, flg)
c by Yoshitaro Takaesu -Aug/15/2011 @KEK
      implicit none

      integer i,j,n
      real*8 p(0:3,n),drcut
      integer flg

      real*8 dr,dRyij
      external dRyij
     
      flg = 0
     
      do i = 1,n-1
         do j = i+1,n
            dr = dRyij(p(0,i),p(0,j))
            if (dr.lt.drcut) then
               flg = 1
               return
            endif
         enddo
      enddo

      end


      subroutine rcut_chk(n,p,rcut, flg)
c by Yoshitaro Takaesu -Sep/3/2013 @KIAS
      implicit none
      integer i,j,n
      real*8 p(0:3,n),rcut
      integer flg
      real*8 r,Ryij
      external Ryij
     
      flg = 0
     
      do i = 1,n-1
         do j = i+1,n
            r = Ryij(p(0,i),p(0,j))
            if (r.lt.rcut) then
               flg = 1
               return
            endif
         enddo
      enddo

      end


      subroutine lptcut_chk(n,p,lptcut, flg)
c by Yoshitaro Takaesu -Aug/17/2011 @KEK
      implicit none

      integer i,j,n
      real*8 p(0:3,n),lptcut
      integer flg

      real*8 ptmax

      real*8 pt
      external pt
     
      flg = 0
      ptmax = 0
      
      do i = 1,n
         if (pt(p(0,i)).gt.ptmax) then
            ptmax = pt(p(0,i))
         endif
      enddo

      if (ptmax.lt.lptcut) flg = 1

      end

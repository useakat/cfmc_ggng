      subroutine ptcut_chk(n,p,ptcut, flg)
      implicit none

      integer i,n
      real*8 p(0:3,n),ptcut
      integer flg

      real*8 ppt
      external ppt
     
      flg = 0
     
      do i = 1,n
         if( ppt(p(0,i)).lt.ptcut ) flg = 1
      enddo

      end


      subroutine etacut_chk(n,p,etacut, flg)
      implicit none

      integer i,n
      real*8 p(0:3,n),etacut
      integer flg

      real*8 peta
      external peta
     
      flg = 0
     
      do i = 1,n
         if( abs(peta(p(0,i))).gt.etacut ) flg = 1
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

      return
      end


      subroutine mijcut_chk(n,p,mijcut, flg)
c by Yoshitaro Takaesu -Aug/04/2011 @KEK
      implicit none
      
      integer i,j
      
      integer n
      real*8 p(0:3,n),mijcut
      integer flg

      real*8 mij

      real*8 pmij
      external pmij

      flg = 0
      
      do i = 1,n-1
         do j = i+1,n
            mij = pmij(p(0,i),p(0,j))     
            if( mij.lt.mijcut ) then
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

      real*8 dr,pdeltar
      external pdeltar
     
      flg = 0
     
      do i = 1,n-1
         do j = i+1,n
            dr = pdeltar(p(0,i),p(0,j))
            if (dr.lt.drcut) then
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

      real*8 dr,ppt
      external ppt
     
      flg = 0
      ptmax = 0
      
      do i = 1,n
         if (ppt(p(0,i)).gt.ptmax) then
            ptmax = ppt(p(0,i))
         endif
      enddo

      if (ptmax.lt.lptcut) flg = 1

      end

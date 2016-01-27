      subroutine write_mom(n,p,lun)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK
      implicit none

      integer i
      integer n,lun
      real*8 p(0:3,n)
  
      write(lun,*) "momenta"
      do i = 1,n
         write(lun,*) p(0,i),p(1,i),p(2,i),p(3,i)
      enddo
      write(lun,*)

      return
      end


      subroutine write_wav(n,w,lun)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK
      implicit none

      integer i
      integer n,k,lun
      complex*16 w(6,n)
  
      write(lun,*) "wave functions"
      do i = 1,n
         write(lun,*) w(1,i),w(2,i),w(3,i),w(4,i),w(5,i),w(6,i)
      enddo
      write(lun,*)

      return
      end


      subroutine write_wav18(n,w,lun)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK
      implicit none

      integer i
      integer n,lun
      complex*16 w(18,n)
  
      write(lun,*) "wave functions"
      do i = 1,n
         write(lun,*) w(1,i),w(2,i),w(3,i),w(4,i)
         write(lun,*) w(5,i),w(6,i),w(7,i),w(8,i)
         write(lun,*) w(9,i),w(10,i),w(11,i),w(12,i)
         write(lun,*) w(13,i),w(14,i),w(15,i),w(16,i)
         write(lun,*) w(17,i),w(18,i)         
      enddo
      write(lun,*)

      return
      end

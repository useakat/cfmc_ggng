      subroutine save_histo(ii)
      implicit none
C     GLOBAL VARIABLES
      include 'cparam.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i
      character*4 iii
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (ich.lt.10) then
         open(31,file='bases_data/x1_'//ii(1:1)//'.top'
     &        ,status='replace')
         open(32,file='bases_data/x2_'//ii(1:1)//'.top'
     &        ,status='replace')
         open(33,file='bases_data/m2summ2_'//ii(1:1)//'.top'
     &        ,status='replace')

         do i = 1,ndim
            write(iii,*) i 
            if (i.lt.10) then
               open(33+i,file='bases_data/z'//iii(2:2)//'_'
     &              //ii(1:1)//'.top',status='replace')
            elseif (i.lt.100) then
               open(33+i,file='bases_data/z'//iii(2:3)//'_'
     &              //ii(1:1)//'.top',status='replace')
            endif
         enddo

         do i = 1,ndist
            write(iii,*) i
            if (i.lt.10) then
               open(40+i,file='bases_data/dist'//iii(2:2)//'_'
     &              //ii(1:1)//'.top',status='replace')
            elseif (i.lt.100) then
               open(40+i,file='bases_data/dist'//iii(2:3)//'_'
     &              //ii(1:1)//'.top',status='replace')
            endif
         enddo

      elseif (ich.lt.100) then
         open(31,file='bases_data/x1_'//ii(1:2)//'.top'
     &        ,status='replace')
         open(32,file='bases_data/x2_'//ii(1:2)//'.top'
     &        ,status='replace')
         open(33,file='bases_data/m2summ2_'//ii(1:2)//'.top'
     &        ,status='replace')

         do i = 1,ndim
            write(iii,*) i 
            if (i.lt.10) then
               open(33+i,file='bases_data/z'//iii(2:2)//'_'
     &              //ii(1:2)//'.top',status='replace')
            elseif (i.lt.100) then
               open(33+i,file='bases_data/z'//iii(2:3)//'_'
     &              //ii(1:2)//'.top',status='replace')
            endif
         enddo

         do i = 1,ndist
            write(iii,*) i
            if (i.lt.10) then
               open(40+i,file='bases_data/dist'//iii(2:2)//'_'
     &              //ii(1:2)//'.top',status='replace')
            elseif (i.lt.100) then
               open(40+i,file='bases_data/dist'//iii(2:3)//'_'
     &              //ii(1:2)//'.top',status='replace')
            endif
         enddo
      endif      

      do i = 31,33
         call xhsave(i,i)
         close(i)
      enddo
      do i = 34,33+ndim
         call xhsave(i,i)         
         close(i)
      enddo
      do i = 41,40+ndist
         call xhsave(i,i)
         close(i)
      enddo

      return
      end

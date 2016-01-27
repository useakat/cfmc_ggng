      program get_std
C     ****************************************************
C     
C     The driver for calculating standard deviation from 
C     a data file.
C
C     By Yoshitaro Takaesu @KEK Dec.24 2011
C     ****************************************************
      implicitnone
C     
C     LOCAL VARIABLES 
C     
      character*20 ffile
      character*2 iirow,nnrow
      integer nrow,irow
      real*8 mean,std

      call getarg(1,ffile)
      call getarg(2,nnrow)
      read(nnrow,*) nrow
      call getarg(3,iirow)
      read(iirow,*) irow

      call calc_std(ffile,nrow,irow,mean,std)

      write(*,*) mean,std,std/mean*100

      end


      subroutine calc_std(filename,nrow,irow, mean,sigma)
C     ****************************************************
C     
C     The subroutine for calculating standard deviation from 
C     a data file.
C
C     By Yoshitaro Takaesu @KEK Dec.24 2011
C     ****************************************************
      implicitnone
C     
C     CONSTANTS
C     
      integer maxline
      parameter(maxline=100000)
C     
C     LOCAL VARIABLES 
C     
      integer i
      character*20 filename
      integer nrow,line,irow
      real*8 dataa(nrow),total,a(maxline),sigma,mean,sigma2
C     ----------
C     BEGIN CODE
C     ----------
      open(1,file=filename,status='old')
      line = 0
      total = 0d0
      do i = 1,maxline
         read(1,*,end=99) dataa
         line = line +1
         a(line) = dataa(irow)
         total = total +a(line)
      enddo
 99   close(1)

      mean = total/dble(line)

      sigma2 = 0d0
      do i = 1,line
         sigma2 = sigma2 +( a(i) -mean )**2/dble(line-1)
      enddo

      sigma = dsqrt(sigma2)


      return
      end

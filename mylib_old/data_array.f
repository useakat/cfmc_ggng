      subroutine data_array(a,ndim,name,lname,lun)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 25 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer ndim,lname,lun
      integer a(ndim)
      character*99 name
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (lname.lt.10) then
         call data_array1(a,ndim,name,lname,lun)
      elseif (lname.lt.100) then
         call data_array2(a,ndim,name,lname,lun)
      endif

      return
      end


      subroutine data_array1(a,ndim,name,lname,lun)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 25 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer ndim,lname,lun
      integer a(ndim)
      character*9 name
      character*50 cfmt,title_fmt,end_fmt,mid_fmt
C     LOCAL VARIABLES 
      integer i,j,k,nrem,ltitle,lmid,lend
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      title_fmt = "6x,a4,1x,a1,1x,a1" 
      ltitle = 17
      mid_fmt = "5x,a1,5x,1x,1x" 
      lmid = 14
      end_fmt = "1x,a1"  
      lend = 5
      write(title_fmt(11:11),'(i1)') lname
      write(mid_fmt(10:10),'(i1)') lname

      if (ndim.le.10) then
         if (ndim-1.eq.0) then
            cfmt = "("//title_fmt(1:ltitle)//",i2,"
     &           //end_fmt(1:lend)//")"
            write(lun,cfmt(1:ltitle+lend+6)) "data",name(1:lname),
     &        "/",a,"/"
         else
            cfmt = "("//title_fmt(1:ltitle)//",i2,1(',',i2),"
     &           //end_fmt(1:lend)//")"
            write(cfmt(ltitle+6:ltitle+6),'(i1)') ndim-1
            write(lun,cfmt(1:ltitle+lend+16)) "data",name(1:lname),
     &        "/",a,"/"
         endif
      else
         i = ndim/10
         do j = 1,i
            if (j.eq.1) then
               cfmt = "("//title_fmt(1:ltitle)//",i2,9(',',i2))" 
               write(lun,cfmt(1:ltitle+15)) "data",name(1:lname),"/"
     &              ,(a(k),k=1,10)
            else
               cfmt = "("//mid_fmt(1:lmid)//",10(',',i2))"
               write(lun,cfmt(1:lmid+13)) "&",(a(k),k=10*(j-1)+1,10*j)
            endif
         enddo
         nrem = ndim -10*i
         if (nrem.eq.0) then
            cfmt = "("//mid_fmt(1:lmid)//","
     &           //end_fmt(1:lend)//")"
            write(lun,cfmt(1:lmid+lend+2)) "&",(a(j),j=10*i+1,ndim),"/"
         else
            cfmt = "("//mid_fmt(1:lmid)//",1(',',i2),"
     &           //end_fmt(1:lend)//")"
            write(cfmt(lmid+3:lmid+3),'(i1)') nrem
            write(lun,cfmt(1:lmid+lend+13)) "&",(a(j),j=10*i+1,ndim),"/"
         endif
      endif

      return
      end


      subroutine data_array2(a,ndim,name,lname,lun)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 25 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer ndim,lname,lun
      integer a(ndim)
      character*99 name
      character*50 cfmt,title_fmt,end_fmt,mid_fmt
C     LOCAL VARIABLES 
      integer i,j,k,nrem,ltitle,lmid,lend
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      title_fmt = "6x,a4,1x,a10,1x,a1" 
      ltitle = 18
      mid_fmt = "5x,a1,5x,10x,1x" 
      lmid = 15
      end_fmt = "1x,a1"  
      lend = 5
      write(title_fmt(11:12),'(i2)') lname
      write(mid_fmt(10:11),'(i2)') lname

      if (ndim.le.10) then
         cfmt = "("//title_fmt(1:ltitle)//",i2,1(',',i2),"
     &        //end_fmt(1:lend)//")"
         write(cfmt(ltitle+6:ltitle+6),'(i1)') ndim-1
         write(lun,cfmt(1:ltitle+lend+16)) "data",name,"/",a,"/"
      else
         i = ndim/10
         do j = 1,i
            if (j.eq.1) then
               cfmt = "("//title_fmt(1:ltitle)//",i2,9(',',i2))" 
               write(lun,cfmt(1:ltitle+15)) "data",name,"/"
     &              ,(a(k),k=1,10)
            else
               cfmt = "("//mid_fmt(1:lmid)//",10(',',i2))"
               write(lun,cfmt(1:lmid+13)) "&",(a(k),k=10*(j-1)+1,10*j)
            endif
         enddo
         nrem = ndim -10*i
         cfmt = "("//mid_fmt(1:lmid)//",1(',',i2),"
     &        //end_fmt(1:lend)//")"
         write(cfmt(lmid+3:lmid+3),'(i1)') nrem
         write(lun,cfmt(1:lmid+lend+13)) "&",(a(j),j=10*i+1,ndim),"/"
      endif

      return
      end

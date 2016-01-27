      program combine_dat
      implicit none
C     GLOBAL VARIABLES
      include 'cparam.inc'
      include 'len_run_name.inc'
C     ARGUMENTS 
      character*3 nndim,nndist,cchstart,cchp1
      character*(len_run_name) run_name
      integer ichstart,ichp1
C     LOCAL VARIABLES 
      integer max_line,nline,i,j,k,n
      parameter(max_line=1000)
      real*8 x(max_line),sumy(max_line),sumdy2(max_line),y,dy
C     ----------
C     BEGIN CODE
C     ----------
      call getarg(1,run_name)
      call getarg(2,nndim)
      call getarg(3,cchstart)
      call getarg(4,cchp1)
      read(nndim,*) ndim
      read(cchstart,*) ichstart
      read(cchp1,*) ichp1
c      call getarg(3,nndist)
c      read(nndist,*) ndist

      do i = 1,max_line
         sumy(i) = 0d0
         sumdy2(i) = 0d0
      enddo
      do i = ichstart,ichp1-1
         write(ii,*) i
         if (i.lt.10) then
            open(1,file="rslt_"//run_name//"/data/x1_"//ii(2:2)//
     &           ".dat"
     &           ,status="old")
         elseif (i.lt.100) then
            open(1,file="rslt_"//run_name//"/data/x1_"//ii(2:3)//
     &           ".dat"
     &           ,status="old")
         elseif (i.lt.1000) then
            open(1,file="rslt_"//run_name//"/data/x1_"//ii(2:4)//
     &           ".dat"
     &           ,status="old")
         endif
         nline = 0
         do j = 1,max_line
            nline = nline +1
            read(1,*,end=99) x(j),y,dy
            sumy(j) = sumy(j) +y
            sumdy2(j) = sumdy2(j) +dy**2
         enddo
 99      close(1)
      enddo      
      open(2,file="rslt_"//run_name//"/data/x1.dat"
     &        ,status="replace")
      do i = 1,nline
         write(2,*) x(i),sumy(i),dsqrt(sumdy2(i))
      enddo
      close(2)

      do i = 1,max_line
         sumy(i) = 0d0
         sumdy2(i) = 0d0
      enddo
      do i = ichstart,ichp1-1
         write(ii,*) i
         if (i.lt.10) then
            open(1,file="rslt_"//run_name//"/data/x2_"//ii(2:2)//
     &           ".dat",status="old")
         elseif (i.lt.100) then
            open(1,file="rslt_"//run_name//"/data/x2_"//ii(2:3)//
     &           ".dat",status="old")
         elseif (i.lt.1000) then
            open(1,file="rslt_"//run_name//"/data/x2_"//ii(2:4)//
     &           ".dat",status="old")
         endif
            nline = 0
         do j = 1,max_line
            nline = nline +1
            read(1,*,end=100) x(j),y,dy
            sumy(j) = sumy(j) +y
            sumdy2(j) = sumdy2(j) +dy**2
         enddo
 100     close(1)
      enddo      
      open(2,file="rslt_"//run_name//"/data/x2.dat"
     &        ,status="replace")
      do i = 1,nline
         write(2,*) x(i),sumy(i),dsqrt(sumdy2(i))
      enddo
      close(2)

      do i = 1,max_line
         sumy(i) = 0d0
         sumdy2(i) = 0d0
      enddo
      do i = ichstart,ichp1-1
         write(ii,*) i
         if (i.lt.10) then
            open(1,file="rslt_"//run_name//"/data/m2summ2_"
     &           //ii(2:2)//".dat",status="old")
         elseif (i.lt.100) then
            open(1,file="rslt_"//run_name//"/data/m2summ2_"
     &           //ii(2:3)//".dat",status="old")
         elseif (i.lt.1000) then
            open(1,file="rslt_"//run_name//"/data/m2summ2_"
     &           //ii(2:4)//".dat",status="old")
         endif
         nline = 0
         do j = 1,max_line
            nline = nline +1
            read(1,*,end=200) x(j),y,dy
            sumy(j) = sumy(j) +y
            sumdy2(j) = sumdy2(j) +dy**2
         enddo
 200     close(1)
      enddo      
      open(2,file="rslt_"//run_name//"/data/m2summ2.dat"
     &        ,status="replace")
      do i = 1,nline
         write(2,*) x(i),sumy(i),dsqrt(sumdy2(i))
      enddo
      close(2)

c      do k = 1,ndim
c         write(kk,*) k
c         do i = 1,max_line
c            sumy(i) = 0d0
c            sumdy2(i) = 0d0
c         enddo
c         do i = 1,nch
c            write(ii,*) i
c            if (k.lt.10) then
c               open(1,file="bases_plots/"//run_name//"/z"//kk(2:2)//"_"
c     &              //ii(2:2)//".dat",status="old")
c            elseif (k.lt.100) then
c               open(1,file="bases_plots/"//run_name//"/z"//kk(2:3)//"_"
c     &              //ii(2:2)//".dat",status="old")
c            endif
c            nline = 0
c            do j = 1,max_line
c               nline = nline +1
c               read(1,*,end=120) x(j),y,dy
c               sumy(j) = sumy(j) +y
c               sumdy2(j) = sumdy2(j) +dy**2
c            enddo
c 120        close(1)
c         enddo      
c         if (k.lt.10) then
c            open(2,file="bases_plots/"//run_name//"/z"//kk(2:2)//".dat"
c     &           ,status="replace")
c         elseif (k.lt.100) then
c            open(2,file="bases_plots/"//run_name//"/z"//kk(2:3)//".dat"
c     &           ,status="replace")
c         endif
c         do i = 1,nline
c            write(2,*) x(i),sumy(i),dsqrt(sumdy2(i))
c         enddo
c         close(2)
c      enddo

      do k = 1,ndist
         write(kk,*) k
         do i = 1,max_line
            sumy(i) = 0d0
            sumdy2(i) = 0d0
         enddo
         do i = ichstart,ichp1-1
            write(ii,*) i
            if (k.lt.10) then
               if (i.lt.10) then
                  open(1,file="rslt_"//run_name//"/data/dist"
     &                 //kk(2:2)//"_"//ii(2:2)//".dat",status="old")
               elseif (i.lt.100) then
                  open(1,file="rslt_"//run_name//"/data/dist"
     &                 //kk(2:2)//"_"//ii(2:3)//".dat",status="old")
               elseif (i.lt.1000) then
                  open(1,file="rslt_"//run_name//"/data/dist"
     &                 //kk(2:2)//"_"//ii(2:4)//".dat",status="old")
               endif
            elseif (k.lt.100) then
               if (i.lt.10) then
                  open(1,file="rslt_"//run_name//"/data/dist"
     &                 //kk(2:3)//"_"//ii(2:2)//".dat",status="old")
               elseif (i.lt.100) then
                  open(1,file="rslt_"//run_name//"/data/dist"
     &                 //kk(2:3)//"_"//ii(2:3)//".dat",status="old")
               elseif (i.lt.1000) then
                  open(1,file="rslt_"//run_name//"/data/dist"
     &                 //kk(2:3)//"_"//ii(2:4)//".dat",status="old")
               endif
            endif
            nline = 0
            do j = 1,max_line
               nline = nline +1
               read(1,*,end=130) x(j),y,dy
               sumy(j) = sumy(j) +y
               sumdy2(j) = sumdy2(j) +dy**2
            enddo
 130        close(1)
         enddo   
         if (k.lt.10) then
            open(2,file="rslt_"//run_name//"/data/dist"
     &           //kk(2:2)//".dat",status="replace")
         elseif (k.lt.100) then
            open(2,file="rslt_"//run_name//"/data/dist"
     &           //kk(2:3)//".dat",status="replace")
         endif
         do i = 1,nline
            write(2,*) x(i),sumy(i),dsqrt(sumdy2(i))
         enddo
         close(2)
      enddo

      return
      end

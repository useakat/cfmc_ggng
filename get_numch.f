      program get_numch
C     ****************************************************
C     The driver for reading the number of channels from
C     cparam.inc.
C
C     By Yoshitaro Takaesu @KEK Nov.4 2012
C     ****************************************************
      implicitnone
C     GLOBAL VARIABLES
      include 'cparam.inc'
C     ----------
C     BEGIN CODE
C     ----------
      ndim = 3*nfin -4 +ipdf*2

      open(1,file="nch.dat",status="replace")
      write(1,*) nch,isingle_ch
      close(1)

      open(1,file="ndist.dat",status="replace")
      write(1,*) ndist
      close(1)

      open(1,file="ndim.dat",status="replace")
      write(1,*) ndim
      close(1)

      open(2,file="ipdf.dat",status="replace")
      write(2,*) ipdf
      close(2)

      open(2,file="next.dat",status="replace")
      write(2,*) next,nini,nfin
      close(2)

      end

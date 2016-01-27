      program make_files
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN.23 2013
C     ****************************************************
      implicitnone
C     CONSTANTS     
C     ARGUMENTS 
      character*2 cnext
      integer next
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer nini,nfin,maxnch
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      call getarg(1,cnext)
      read(cnext,*) next 
      maxnch = 200
      nini = 2
      nfin = next -nini

      call make_cparam(nini,nfin,maxnch)
c      call make_bfunc

      end

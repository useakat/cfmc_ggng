      subroutine hmffrd
      implicit none
      include 'hmffdt.inc'
      include 'cparam.inc'
      
      character*80 cline0,cline1
      integer nwdmax
      parameter (nwdmax= 2)
      character*80 cword(nwdmax)
      integer lword(nwdmax) 
      integer nword
      character*1 charc/'c'/, cblk/' '/

      nkeys = 0

 100  continue
      read(lustdi,'(a)',end=200) cline0
      call tolow(cline0,cline1)

      if (cline1(1:1).ne.'c') then
      
         call comsep(cline1,nwdmax,nword,cword,lword)

         if (nword.eq.nwdmax) then
            nkeys = nkeys + 1
            ckey(nkeys) = cblk
            cval(nkeys) = cblk
            ckey(nkeys) = cword(1)(:lword(1))
            cval(nkeys) = cword(2)(:lword(2))
         endif

      endif
      goto 100

 200  continue

      return
      end

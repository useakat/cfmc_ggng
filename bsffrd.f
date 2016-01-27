      subroutine bsffrd
      implicit none
      include 'hmffdt.inc'
      include 'cparam.inc'
      integer ikey
      integer i,lnblnk

      do ikey=1,nkeys
         if(ckey(ikey).eq.'outf') then
            bsoutf = cval(ikey)
         else if(ckey(ikey).eq.'hstf') then
            bshstf = cval(ikey)
         else if(ckey(ikey).eq.'evtf') then
            spevtf = cval(ikey)
         else if(ckey(ikey).eq.'nvsp') then
            read(cval(ikey),*) nospev
         endif
      enddo

      if (ifprog.eq.1) then
         if(lnblnk(bsoutf).gt.0) then
            open(lubsf, file=bsoutf,status='unknown',form='unformatted')
         endif
      elseif (ifprog.eq.2) then
         if(lnblnk(bsoutf).gt.0) then
            open(lubsf, file=bsoutf,status='old',form='unformatted')
         else
            write(luerr,*) 'bsffrd: run spring without bases data'
            stop
         endif
         if(lnblnk(spevtf).gt.0) then
            open(luevtf, file=spevtf,status='unknown',
     &           form='formatted')
         endif
      endif

      if(lnblnk(bshstf).gt.0) then
         open(luhist, file=bshstf,status='unknown',
     &        form='formatted')
      endif

      return
      end

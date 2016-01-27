      subroutine impose_cuts(ierr)

      include 'cparam.inc'
      include 'pshk.inc'
      include 'bfunc.inc'

      ierr = 0
c      call ptorder_chk(nfin,Q(0,3),ierr2)
c      if( ierr2.ne.0 ) return
      call ptcut_chk(nfin,Q(0,3),ptcut,ierr)
      if( ierr.ne.0 ) return
      call etacut_chk(nfin,Q(0,3),ayecut,ierr)
      if( ierr.ne.0 ) return
      call ptijcut_chk(nfin,Q(0,3),ptijcut,ierr)
      if( ierr.ne.0 ) return      
      call drcut_chk(nfin,Q(0,3),rcut,ierr)
      if( ierr.ne.0 ) return      
      call lptcut_chk(nfin,Q(0,3),lptcut,ierr)
      if ( ierr.ne.0 ) return
c      call mijcut_chk(nfin,Q(0,3),8d0,ierr)
c      if ( ierr.ne.0 ) return

      return
      end

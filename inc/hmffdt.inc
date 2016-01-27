      integer nkeys
      common /hmffdt/ nkeys
      integer nkymax
      parameter (nkymax = 100)
      character*4 ckey
      character*80 cval
      common /hmffd0/ ckey(nkymax),cval(nkymax)
      save /hmffdt/, /hmffd0/

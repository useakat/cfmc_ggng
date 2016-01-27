      double precision puq(0:3),pux(0:3)
      double precision pini(0:3,nini)
      equivalence (puq(0),pini(0,1))
      equivalence (pux(0),pini(0,2))
      integer nhini(nini)
      integer idini(nini)

      double precision pfin(0:3,nfin)
      integer nhfin(nfin),idfin(nfin)

      double precision pj(0:3,nj)
      equivalence (pj(0,1),pfin(0,1))
      integer nhj(nj),idj(nj)
      equivalence (nhj(1),nhfin(1))
      equivalence (idj(1),idfin(1))

      double precision p(0:3,nexternal)
      integer nhel(nexternal),idpdg(nexternal)
      common /cfinal/ p,nhel,idpdg
      equivalence (nhini(1),nhel(1))
      equivalence (pini(0,1),p(0,1))
      equivalence (idini(1),idpdg(1))
      equivalence (nhfin(1),nhel(nini+1))
      equivalence (pfin(0,1),p(0,nini+1))
      equivalence (idfin(1),idpdg(nini+1))
      save /cfinal/

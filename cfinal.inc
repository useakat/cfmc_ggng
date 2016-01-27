      double precision puq(0:3),pux(0:3)
      integer idini(nini),ihini(nini)
      integer idfin(nfin),ihfin(nfin),icfin(nfin)
      double precision pini(0:3,nini),pfin(0:3,nfin)
      integer nhel(next)
      common/cfinal/idini,idfin,icfin,nhel
      equivalence (ihini(1),nhel(1))
      equivalence (ihfin(1),nhel(nini+1))
      save /cfinal/

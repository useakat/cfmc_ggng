      double precision function pdeltar(p1,p2)

      implicit none

      double precision p1(0:3),p2(0:3)

      double precision eta1,eta2
      double precision phi1,phi2

      double precision peta,pdelphi
      external peta,pdelphi

      eta1 = peta(p1)
      eta2 = peta(p2)

      phi1 = datan2(p1(2),p1(1))
      phi2 = datan2(p2(2),p2(1))

      pdeltar = dsqrt((eta1-eta2)**2+pdelphi(phi1,phi2)**2)

      return
      end

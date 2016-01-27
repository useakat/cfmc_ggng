      double precision function pdelphi(phi1, phi2)

      implicit none

      double precision phi1, phi2

      double precision x1,y1
      double precision x2,y2

      double precision sinphi,cosphi
      double precision phis,phic
      
      x1 = dcos(phi1)
      y1 = dsin(phi1)

      x2 = dcos(phi2)
      y2 = dsin(phi2)

      sinphi = x1*y2-y1*x2
      cosphi = x1*x2+y1*y2

      phis = dasin(sinphi)
      phic = dacos(cosphi)

      if (cosphi .ge. 0.d0) then
         pdelphi = phis
      else if (sinphi .ge. 0.d0) then
         pdelphi = phic
      else 
         pdelphi =  -phic
      endif

      return
      end

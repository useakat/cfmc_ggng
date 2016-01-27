      double precision ptcutt
      parameter (ptcutt = 20.d0)
      double precision etacutt
      parameter (etacutt = 2.5d0)
      double precision drcutt
      parameter (drcutt = 0.4d0)
        
      double precision whatmin, taulmin
      common /ccuts/ whatmin, taulmin
      save /ccuts/

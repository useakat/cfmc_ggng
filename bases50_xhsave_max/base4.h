      INTEGER NDMX,LENG
c      PARAMETER (NDMX = 50, LENG = 32768)
      PARAMETER (NDMX = 50, LENG = 65536)
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),
     .               ND,NG,NPG,MA(MXDIM)
      DOUBLE PRECISION XI,DX,DXD,DXP
      INTEGER ND,NG,NPG,MA
      SAVE /BASE4/
      subroutine lumfun(iset,ipar1,ipar2,Q,sqrts,shat,lf,erro)
C This subroutine calculate luminosity function.
      implicit none 

      include 'bases.inc'

      integer ipar1,ipar2,i,iset
      real*8 Q,ymax,ymin,s,shat,sqrts,lf,erro 

      integer ipar1_ex, ipar2_ex
      real*8 Q_ex,tau             
      common /export/ ipar1_ex, ipar2_ex, Q_ex
      common /para/ tau
 
      ipar1_ex=ipar1
      ipar2_ex=ipar2 
      Q_ex=Q 

      tau = (shat/sqrts)**2 
      call SetCT10(iset)

      call bsinit
      call userin
      call bases(bfunc, estim,sigma,ctime,it1,it2)
      
      lf = estim
      erro = sigma

      return 
      end             


      real*8 function bfunc(z)

      implicit none 

      real*8 z(1)

      integer ipar1, ipar2
      real*8 Q,tau
      common /export/ ipar1, ipar2, Q
      common /para/ tau

      real*8 y,x1,x2,CT10Pdf
      external CT10Pdf

      bfunc = 0d0
      y = z(1)
      x1 = dsqrt(tau)*dexp(y) 
      x2 = dsqrt(tau)*dexp(-y)

      if(x1.ge.1.d0) then 
        bfunc = 0.d0
      elseif(x2.ge.1.d0) then 
          bfunc = 0.d0 
      else
       if(ipar1.eq.ipar2) then
          bfunc = CT10Pdf(ipar1, x1,Q)*CT10Pdf(ipar2,x2,Q)
       else 
         bfunc = CT10Pdf(ipar1, x1,Q)*CT10Pdf(ipar2,x2,Q)
     & + CT10Pdf(ipar2, x1,Q)*CT10Pdf(ipar1,x2,Q)
       endif
      endif

      return 
      end 

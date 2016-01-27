      real*8 function dacosh(x)
C     ****************************************************
C     
C     By Yoshitaro Takaesu @KEK Jun.20 2012
C     
C     ArcCosh function
C
C     Input: 
C            real*8 x 
C     Output:
C
C     ****************************************************
      implicitnone

      real*8 x
C     ----------
C     BEGIN CODE
C     ----------
      if (x.gt.1) then
         dacosh = dlog( x +dsqrt( abs(x**2 -1) ) )
      else
         write(*,*) "Argument is less than 1. Cannot be evaluated."
      endif

      return
      end

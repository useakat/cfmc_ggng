	  SUBROUTINE GET_MOMENTA(NEXTERNAL,ENERGY,PMASS,P)
C---- auxiliary function to change convention between madgraph and rambo
c---- four momenta. 	  
	  IMPLICIT NONE
c	  INCLUDE "nexternal.inc"
C	  ARGUMENTS
	  INTEGER NEXTERNAL
	  REAL*8 ENERGY,PMASS(NEXTERNAL),P(0:3,NEXTERNAL),PRAMBO(4,10),WGT
C         LOCAL
         INTEGER I
         REAL*8 etot2,mom,m1,m2,e1,e2

	 include 'cparam.inc'

         ETOT2=energy**2
         m1=pmass(1)
         m2=pmass(2)
         mom=(Etot2**2 - 2*Etot2*m1**2 + m1**4 -
     -        2*Etot2*m2**2 - 2*m1**2*m2**2 + m2**4)/(4.*Etot2)
         mom=dsqrt(mom)
         e1=DSQRT(mom**2+m1**2)
         e2=DSQRT(mom**2+m2**2)
         write (*,*) e1+e2,mom

         if(nini.eq.2) then

            P(0,1)=e1
            P(1,1)=0d0
            P(2,1)=0d0
            P(3,1)=mom

            P(0,2)=e2
            P(1,2)=0d0
            P(2,2)=0d0
            P(3,2)=-mom
             
            call rambo(nexternal-2,energy,pmass(3),prambo,WGT)
            DO I=3, NEXTERNAL
               P(0,I)=PRAMBO(4,I-2)	
               P(1,I)=PRAMBO(1,I-2)
               P(2,I)=PRAMBO(2,I-2)
               P(3,I)=PRAMBO(3,I-2)	
            ENDDO
             
          elseif(nini.eq.1) then 
             
             P(0,1)=energy
             P(1,1)=0d0
             P(2,1)=0d0
             P(3,1)=0d0
             
             call rambo(nexternal-1,energy,pmass(2),prambo,WGT)
             DO I=2, NEXTERNAL
                P(0,I)=PRAMBO(4,I-1)	
                P(1,I)=PRAMBO(1,I-1)
                P(2,I)=PRAMBO(2,I-1)
                P(3,I)=PRAMBO(3,I-1)	
             ENDDO
          endif
          
	  RETURN
	  END

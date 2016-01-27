      subroutine uggggx(ga,gb,gc,gd,g,g2,Tmass,Twidth , ft)
c
c This subroutine computes an off-shell Tensor boson wavfunction from
c the 4 gluons-Kaluza-Klein tensor boson coupling, corresponding 
c to the color structure f^{a,b,e} f{c,d,e}. 
c
c To optain the complete amplitude, this subroutine must be called three
c times (once for each color structure) with the following permutations:
c     call uggggx(ga,gb,gc,gd,g,g2,Tmass,Twidth,ft)
c     call uggggx(ga,gc,gd,gb,g,g2,Tmass,Twidth,ft)
c     call uggggx(ga,gd,gb,gv,g,g2,Tmass,Twidth,ft)
c     corresponding to 
c	f^{a,b,e} f{c,d,e}
c	f^{a,c,e} f{d,b,e}
c	f^{a,d,e} f{b,c,e}
c
c input:
c       complex ga(6)          : first input gluon          ga
c       complex gb(6)          : second input gluon         gb
c       complex gc(6)          : third input gluon          gc
c       complex gd(6)          : fourth input gluon         gd
c
c       real    g              : coupling constant    gs
c       complex g2             : coupling constant    -k/2
c       real    Tmass          : mass  of OUTPUT KK Tensor T'
c       real    Twidth         : width of OUTPUT KK Tensor boson T'
c
c output:
c       complex ft(18)       : KK tensor current       j(T:g1,g2,g3,g4)     
c     
c- by Q.Li - Jan. 2008
c
      implicit none
      double complex ga(6),gb(6),gc(6),gd(6),tc(6,4),ft(18)
      double complex E1E3,E1E4,E2E3,E2E4,E1PT,E2PT,E3PT,E4PT

      double precision g,Tmass,Twidth

      integer a,b,i,j

      double precision pa(4), pb(4),pc(4),pd(4),pT(4)
      double precision MET(4,4)
      double complex cZero, d, g2
      double precision rZero, rTwo,PT2,r3
      parameter( rZero = 0.0d0, rTwo = 2.0d0, r3=3.0d0)
      parameter(cZero=(0.0d0,0.0d0))

      
      tc(5,1) = ga(5)+gb(5)+gc(5)+gd(5)
      tc(6,1) = ga(6)+gb(6)+gc(6)+gd(6)

      pa(1) = dreal(ga(5))
      pa(2) = dreal(ga(6))
      pa(3) = dimag(ga(6))
      pa(4) = dimag(ga(5))

      pb(1) = dreal(gb(5))
      pb(2) = dreal(gb(6))
      pb(3) = dimag(gb(6))
      pb(4) = dimag(gb(5))

      pc(1) = dreal(gc(5))
      pc(2) = dreal(gc(6))
      pc(3) = dimag(gc(6))
      pc(4) = dimag(gc(5))

      pd(1) = dreal(gd(5))
      pd(2) = dreal(gd(6))
      pd(3) = dimag(gd(6))
      pd(4) = dimag(gd(5))

     
      pT(1) = dreal(tc(5,1))
      pT(2) = dreal(tc(6,1))
      pT(3) = dimag(tc(6,1))
      pT(4) = dimag(tc(5,1))

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0

      PT2 = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
      
      E1PT = pT(1)*ga(1)-pT(2)*ga(2)-pT(3)*ga(3)-pT(4)*ga(4)
      E2PT = pT(1)*gb(1)-pT(2)*gb(2)-pT(3)*gb(3)-pT(4)*gb(4)
      E3PT = pT(1)*gc(1)-pT(2)*gc(2)-pT(3)*gc(3)-pT(4)*gc(4)
      E4PT = pT(1)*gd(1)-pT(2)*gd(2)-pT(3)*gd(3)-pT(4)*gd(4)
      
      E1E3 = ga(1)*gc(1)-ga(2)*gc(2)-ga(3)*gc(3)-ga(4)*gc(4)
      E1E4 = ga(1)*gd(1)-ga(2)*gd(2)-ga(3)*gd(3)-ga(4)*gd(4)
      E2E3 = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      E2E4 = gb(1)*gd(1)-gb(2)*gd(2)-gb(3)*gd(3)-gb(4)*gd(4)
   
      if ( Tmass.gt.rZero ) then
         d = - 1.0d0/dcmplx( PT2-Tmass**2, Tmass*Twidth )
      else
         d = - 1.0d0/dcmplx( PT2, rZero )
      end if
      
	    
      do a = 1,4
         do b=1,4
            
            tc(a,b) = -2*E2E4*ga(b)*gc(a) + 2*E1E4*gb(b)*gc(a)
     &- 2*E2E4*ga(a)*gc(b) + 2*E1E4*gb(a)*gc(b) + 2*E2E3*ga(b)*gd(a)
     & - 2*E1E3*gb(b)*gd(a) + 2*E2E3*ga(a)*gd(b) - 
     &  2*E1E3*gb(a)*gd(b) - 2*E1E4*E2E3*MET(a,b) + 2*E1E3*E2E4*MET(a,b) 
     &- (2*E1PT*E2E4*E3PT*rtwo*MET(a,b))/(r3*tmass**2) + 
     &  (2*E1E4*E2PT*E3PT*rtwo*MET(a,b))/(r3*tmass**2) 
     &+ (2*E1PT*E2E3*E4PT*rtwo*MET(a,b))/(r3*tmass**2)
     & - (2*E1E3*E2PT*E4PT*rtwo*MET(a,b))/(r3*tmass**2) - 
     &  (E1E4*E2E3*PT2*rtwo*MET(a,b))/(r3*tmass**2) 
     &+ (E1E3*E2E4*PT2*rtwo*MET(a,b))/(r3*tmass**2)
     & + (2*E2E4*E3PT*ga(b)*PT(a))/tmass**2 - 
     &  (2*E2E3*E4PT*ga(b)*PT(a))/tmass**2
     & - (2*E1E4*E3PT*gb(b)*PT(a))/tmass**2
     & + (2*E1E3*E4PT*gb(b)*PT(a))/tmass**2
     & + (2*E1PT*E2E4*gc(b)*PT(a))/tmass**2 - 
     &  (2*E1E4*E2PT*gc(b)*PT(a))/tmass**2
     & - (2*E1PT*E2E3*gd(b)*PT(a))/tmass**2
     & + (2*E1E3*E2PT*gd(b)*PT(a))/tmass**2
     & + (2*E2E4*E3PT*ga(a)*PT(b))/tmass**2 - 
     &  (2*E2E3*E4PT*ga(a)*PT(b))/tmass**2 
     &- (2*E1E4*E3PT*gb(a)*PT(b))/tmass**2
     & + (2*E1E3*E4PT*gb(a)*PT(b))/tmass**2 
     &+ (2*E1PT*E2E4*gc(a)*PT(b))/tmass**2 - 
     &  (2*E1E4*E2PT*gc(a)*PT(b))/tmass**2 
     &- (2*E1PT*E2E3*gd(a)*PT(b))/tmass**2 
     &+ (2*E1E3*E2PT*gd(a)*PT(b))/tmass**2
     & - (4*E1PT*E2E4*E3PT*PT(a)*PT(b))/tmass**4 + 
     &  (4*E1E4*E2PT*E3PT*PT(a)*PT(b))/tmass**4 
     &+ (4*E1PT*E2E3*E4PT*PT(a)*PT(b))/tmass**4
     & - (4*E1E3*E2PT*E4PT*PT(a)*PT(b))/tmass**4 - 
     &  (2*E1E4*E2E3*PT2*PT(a)*PT(b))/tmass**4
     & + (2*E1E3*E2E4*PT2*PT(a)*PT(b))/tmass**4
     & + (2*E1PT*E2E4*E3PT*rtwo*PT(a)*PT(b))/(r3*tmass**4) - 
     &  (2*E1E4*E2PT*E3PT*rtwo*PT(a)*PT(b))/(r3*tmass**4)
     & - (2*E1PT*E2E3*E4PT*rtwo*PT(a)*PT(b))/(r3*tmass**4)
     & + (2*E1E3*E2PT*E4PT*rtwo*PT(a)*PT(b))/(r3*tmass**4) + 
     &  (E1E4*E2E3*PT2*rtwo*PT(a)*PT(b))/(r3*tmass**4)
     & - (E1E3*E2E4*PT2*rtwo*PT(a)*PT(b))/(r3*tmass**4) 
     &+ (4*E1E4*E2E3*PT(a)*PT(b))/tmass**2 - 
     &  (4*E1E3*E2E4*PT(a)*PT(b))/tmass**2

            tc(a,b) = -tc(a,b)*d/2.0d0*g*g*g2
c     2.0 factor from propagator convention

         enddo
      enddo

      ft(1) = tc(1,1)
      ft(2) = tc(1,2)
      ft(3) = tc(1,3)
      ft(4) = tc(1,4)
      ft(5) = tc(2,1)
      ft(6) = tc(2,2)
      ft(7) = tc(2,3)
      ft(8) = tc(2,4)
      ft(9) = tc(3,1)
      ft(10) = tc(3,2)
      ft(11) = tc(3,3)
      ft(12) = tc(3,4)
      ft(13) = tc(4,1)
      ft(14) = tc(4,2)
      ft(15) = tc(4,3)
      ft(16) = tc(4,4)
      ft(17) = tc(5,1)
      ft(18) = tc(6,1)

      return
      end

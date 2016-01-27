      subroutine wrap_pshk(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rnd(3*nfin-2),Q(0:3,next),labQ(0:3),rrnd(3*nfin-2)

      Ierr = 0
      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk cannot be used with no PDF option. 
     &   (wrap_pshk)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshk(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshk(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshkpos(ipos,rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos,ipos
      real*8 rnd(3*nfin-2),Q(0:3,next),labQ(0:3),rrnd(3*nfin-2)

      Ierr = 0
c      delyflag = 1
      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk cannot be used with no PDF option. 
     &   (wrap_pshk)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshkpos(ipos,delyflag,ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshkpos(ipos,delyflag,ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshk2(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos,j
      real*8 rnd(30),Q(0:3,next),rrnd(30)

      ierr = 0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk cannot be used with no PDF option. 
     &   (wrap_pshk2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshk2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshk2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshkpos2(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos,j
      real*8 rnd(30),Q(0:3,next),rrnd(30)

      ierr = 0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk cannot be used with no PDF option. 
     &   (wrap_pshk2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshkpos2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshkpos2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshk2_2(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rrnd(3*nfin-2),rnd(3*nfin-2),Q(0:3,next),labQ(0:3),ipw12

      ierr = 0

c      if (ich.eq.1) ipw12 = 1d0
c      if (ich.eq.2) ipw12 = 2d0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk2_2 cannot be used with no PDF option. 
     &   (wrap_pshk2_2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshk2_2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshk2_2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2_2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshk2_2_0(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rrnd(3*nfin-2),rnd(3*nfin-2),Q(0:3,next),labQ(0:3),ipw12

      ierr = 0

c      if (ich.eq.1) ipw12 = 1d0
c      if (ich.eq.2) ipw12 = 2d0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk2_2 cannot be used with no PDF option. 
     &   (wrap_pshk2_2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshk2_2_0(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshk2_2_0(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2_2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshkpos2_2(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rrnd(3*nfin-2),rnd(3*nfin-2),Q(0:3,next),labQ(0:3),ipw12

      ierr = 0

c      if (ich.eq.1) ipw12 = 1d0
c      if (ich.eq.2) ipw12 = 2d0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk2_2 cannot be used with no PDF option. 
     &   (wrap_pshk2_2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshkpos2_2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshkpos2_2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2_2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshkpos2_2_0(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rrnd(3*nfin-2),rnd(3*nfin-2),Q(0:3,next),labQ(0:3),ipw12

      ierr = 0

c      if (ich.eq.1) ipw12 = 1d0
c      if (ich.eq.2) ipw12 = 2d0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk2_2 cannot be used with no PDF option. 
     &   (wrap_pshk2_2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshkpos2_2_0(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshkpos2_2_0(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2_2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshkpos2_3(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rrnd(3*nfin-2),rnd(3*nfin-2),Q(0:3,next),labQ(0:3),ipw12

      ierr = 0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk2_2 cannot be used with no PDF option. 
     &   (wrap_pshk2_2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshkpos2_3(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshkpos2_3(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2_2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end


      subroutine wrap_pshkpos2_4(rnd,sgncos,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      integer ierr,i,sgncos
      real*8 rrnd(3*nfin-2),rnd(3*nfin-2),Q(0:3,next),labQ(0:3),ipw12

      ierr = 0

c      if (ich.eq.1) ipw12 = 1d0
c      if (ich.eq.2) ipw12 = 2d0

      if (ipdf.eq.0) then
         write(luerr,*) "ERROR:  pshk2_2 cannot be used with no PDF option. 
     &   (wrap_pshk2_2)"
         stop
      elseif (ipdf.eq.1) then
         if (sgncos.eq.1) then
            call pshkpos2_2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         elseif (sgncos.eq.-1) then
            call pshkpos2_2(ss,nfin,rnd,x1,x2,wpsn,jacob,ierr)
c            call pshk2_2_d(ss,nfin,rrnd,ipw12,x1,x2,wpsn,jacob,ierr)
            if(ierr.ne.0) return
         else
            write(luerr,*) "ERROR: sgncos should be 1 or -1. 
     &      (wrap_pshk2_2)"
            stop
         endif
      endif
      call momconv(Q)

      return
      end

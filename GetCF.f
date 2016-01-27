      subroutine GetCF(rnd)
      implicitnone
      include 'cparam.inc'
      integer i,nperm,ipflag
      real*8 rnd

      if( flowtype.eq.1 ) then
         do i = 1,next
            floworder(i,1) = i
         enddo

         nperm = int( rnd*ncf ) + 1
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(floworder(2,1),next-1,ipflag)
            endif
         enddo
         call conjugate_cf(next,floworder(1,1),floworder(1,2))
      endif

      return
      end


      subroutine GetCF_i(nperm)
      implicitnone
      include 'cparam.inc'
      integer i,nperm,ipflag

      if( flowtype.eq.1 ) then
         do i = 1,next
            floworder(i,1) = i
         enddo

         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(floworder(2,1),next-1,ipflag)
            endif
         enddo
         call conjugate_cf(next,floworder(1,1),floworder(1,2))
      endif

      return
      end


      subroutine GetCF_basic(nperm)
      implicitnone
      include 'cparam.inc'
      integer i,nperm,pos

      if ( flowtype.eq.1 ) then
         call gen_basiccf(nperm,next,floworder(1,1))
         call conjugate_cf(next,floworder(1,1),floworder(1,2))
      endif

      return
      end


      subroutine GetCF_opt
      implicitnone
      include 'cparam.inc'
      include 'opt_info.inc'
      integer i,ichopt

      ichopt = per_opt(ich) ! conversion of the system ch number to opt-ch number
      cffact = opt_ncf(ichopt)
      i = int(cffact*rand()) +1
      ichcf = opt_cf(i,ichopt)
      call GetCF_basic(ichcf)

      return
      end


      subroutine GetCF_opt_co1(rnd,icflow)
      implicitnone
      include 'cparam.inc'
      include 'opt_info.inc'
      integer i,j,fact,a(next),iflag,cf(next,next-1)
      integer icfflag,bcf_count,nperm,ibcf(next-1),icflow
      integer ichopt
      real*8 rnd
      external fact

      ichopt = per_opt(ich)
      cffact = opt_ncf(ichopt)
      i = int(cffact*rnd) +1
      ichcf = opt_cf(i,ichopt)
      call GetCF_basic(ichcf)
      
      do i = 1,next
         a(i) = i
      enddo
      do i =1,next-1
         call gen_basiccf(i,next,cf(1,i))
      enddo
      bcf_count = 0

      nperm = fact(next-1)
      do i = 1,nperm
         if (i.ne.1) then
            call ipnext(a,next,iflag)
         endif
         do j = 1,next-1
            call check_cf(next,cf(1,j),a,icfflag)
            if (icfflag.eq.1) then
               bcf_count = bcf_count +1
               ibcf(bcf_count) = i
               exit
            endif
         enddo
      enddo

      icflow = ibcf(ichcf)

      return
      end

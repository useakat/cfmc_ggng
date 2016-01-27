      subroutine GenWgt_per(P)
C
C     This subroutine is for phsk_v6 or before, where the priority is put on 
C     non-splitting leg generations, and a splitting-leg momentum is calculated
C     using 4-momentum conservation.
C
      implicitnone
      include 'cparam.inc'
      include 'opt_info.inc'
      integer i,j,nsch
      real*8 ipw1,ipw2,ipw3,ipwm,P(0:3,next)
      real*8 tjacob,sjacob,sptjacob
      real*8 ps4bd22jacob
      external tjacob,sjacob,sptjacob,ps4bd22jacob

      ipw1 = 1.5d0
      ipw2 = 1.5d0
      ipw3 = 2d0 ! for the pt-jacobian of a splitting leg
      ipwm = 1d0 ! for the m-jacobian of a splitting leg

      do i = 1,nch_per  ! the number of all channels
         AA(i) = 1d0
         if ((perch(1,2,i).gt.0).and.(perch(2,2,i).gt.0)) then
            AA(i) = AA(i)*tjacob(perch(1,2,i),perch(2,2,i),ipw1,ipw1,P)
         elseif (perch(1,2,i).gt.0) then
            AA(i) = AA(i)*tjacob(perch(1,2,i),0,ipw1,0d0,P)
         elseif (perch(2,2,i).gt.0) then
            AA(i) = AA(i)*tjacob(perch(2,2,i),0,ipw1,0d0,P)
         endif
         nsch = 0
         do j = 3,max_natm
            ! searching for s-ch splitting
            if (perch(j,1,i).gt.0) then
               if (perch(j,2,i).gt.0) then
                  ! For found s-ch splitting
                  nsch = nsch +1
                  ! if the pshk leg number for the one of the splitting particles is next
                  if (chmom(perch(j,1,i)-2,i).eq.next) then
                     ! calculate s-sh jacobian without ptij jacobian
                     AA(i) = AA(i)*sjacob(perch(j,1,i),perch(j,2,i),ipwm,P)
                  else
                     ! calculate s-sh jacobian with ptij jacobian
                     AA(i) = AA(i)*sptjacob(perch(j,1,i),perch(j,2,i),ipw3,ipwm,P)
                  endif
               else
                  ! For the found non-splitting leg
                  ! if there is a splitting-leg, or the pshk leg number for this leg is 
                  ! not next
                  if ((nsch.gt.0).or.(j.lt.max_natm)) then
                     AA(i) = AA(i)*tjacob(perch(j,1,i),0,ipw1,0d0,P)
                  endif
               endif
            endif
         enddo
         AA(i) = AA(i)*jfact(perch_type(i))
      enddo

      return
      end


      subroutine GenWgt_per2(P)
C
C     This subroutine is for phsk_v7, where the priority is put on 
C     splitting leg generations, and a non-splitting-leg momentum is 
C     calculated using 4-momentum conservation.
C
      implicitnone
      include 'cparam.inc'
      include 'opt_info.inc'
      integer i,j,nsch
      real*8 ipw1,ipw2,ipw3,ipwm,P(0:3,next)
      real*8 tjacob,sjacob,sptjacob
      real*8 ps4bd22jacob
      external tjacob,sjacob,sptjacob,ps4bd22jacob

      ipw1 = 1.5d0
      ipw2 = 1.5d0
      ipw3 = 1d0
      ipwm = 1d0

      do i = 1,nch_per  ! the number of all channels
         AA(i) = 1d0
         if ((perch(1,2,i).gt.0).and.(perch(2,2,i).gt.0)) then
            AA(i) = AA(i)*tjacob(perch(1,2,i),perch(2,2,i),ipw1,ipw1,P)
         elseif (perch(1,2,i).gt.0) then
            AA(i) = AA(i)*tjacob(perch(1,2,i),0,ipw1,0d0,P)
         elseif (perch(2,2,i).gt.0) then
            AA(i) = AA(i)*tjacob(perch(2,2,i),0,ipw1,0d0,P)
         endif
         nsch = 0
         do j = 3,max_natm
            ! searching for s-ch splitting
            if (perch(j,1,i).gt.0) then
               if (perch(j,2,i).gt.0) then
                  ! For found s-ch splitting
                  nsch = nsch +1
                  ! if the larger particle number of the splitting particles is next
                  if (chmom(perch(j,1,i)-2,i).eq.next) then
                     ! calculate s-sh jacobian without ptij jacobian
                     AA(i) = AA(i)*sjacob(perch(j,1,i),perch(j,2,i),ipwm,P)
                  else
                     ! calculate s-sh jacobian with ptij jacobian
                     AA(i) = AA(i)*sptjacob(perch(j,1,i),perch(j,2,i),ipw3,ipwm,P)
                  endif
               else
                  if ((nsch.gt.0).or.(j.lt.max_natm)) then
                     AA(i) = AA(i)*tjacob(perch(j,1,i),0,ipw1,0d0,P)
                  endif
               endif
            endif
         enddo
         AA(i) = AA(i)*jfact(perch_type(i))
      enddo

      return
      end

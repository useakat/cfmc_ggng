      subroutine analyse_perch(natm,atm,opt_flag)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 24 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'chinfo.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer natm,opt_flag
      integer atm(2,natm)
C     LOCAL VARIABLES 
      integer i,j,nsch,ipos1,ipos2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nsch = 0

      do i = 1,natm
         if (i.gt.2) then
            if (atm(2,i).ne.0) nsch = nsch +1
         endif
      enddo

      if (nsch.eq.0) then
         nch_per = nch_per +1
         ntch_per = ntch_per +1
         perch_type(nch_per) = 0
         per_opt(nch_per) = 0
         if (opt_flag.eq.1) then
            nch_opt = nch_opt +1
            ntch_opt = ntch_opt +1
            optch_type(nch_opt) = 0
            opt_per(nch_opt) = nch_per
            per_opt(nch_per) = nch_opt
         endif
         do i = 1,natm
            do j = 1,2
               perch(i,j,nch_per) = atm(j,i)
            enddo
         enddo
         if ((atm(2,1)+atm(2,2)).eq.2*ngluons-1) then
            ipos_tch(ntch_per) = 2
         else
            if ((atm(2,1).eq.ngluons).or.(atm(2,2).eq.ngluons)) then
               ipos_tch(ntch_per) = 1
            else
               ipos_tch(ntch_per) = 0
            endif
         endif
         call get_chmom(natm,atm,nsch)
      elseif (opt_flag.eq.1) then
         nch_opt = nch_opt +1
         nch_per = nch_per +1
         opt_per(nch_opt) = nch_per
         per_opt(nch_per) = nch_opt
         if (nsch.eq.1) then
            nsch_opt = nsch_opt +1
            optch_type(nch_opt) = 1
            nsch_per = nsch_per +1
            perch_type(nch_per) = 1
         elseif (nsch.eq.2) then
            n2sch_opt = n2sch_opt +1
            n2sch_per = n2sch_per +1
            if (ngluons.eq.6) then
               optch_type(nch_opt) = 20
               perch_type(nch_per) = 20
            else   
               optch_type(nch_opt) = 2
               perch_type(nch_per) = 2
            endif
         elseif (nsch.eq.3) then
            n3sch_opt = n3sch_opt +1
            optch_type(nch_opt) = 3
            n3sch_per = n3sch_per +1
            perch_type(nch_per) = 3
         elseif (nsch.eq.4) then
            n4sch_opt = n4sch_opt +1
            optch_type(nch_opt) = 4
            n4sch_per = n4sch_per +1
            perch_type(nch_per) = 4
         endif
         do i = 1,natm
            do j = 1,2
               perch(i,j,nch_per) = atm(j,i)
            enddo
         enddo
         call get_chmom(natm,atm,nsch)
      endif

      return
      end

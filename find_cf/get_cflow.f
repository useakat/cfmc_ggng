      subroutine get_cflow(next,natm,atm,a,maxk,ncf,cflow)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C     ****************************************************
      implicit none

C     GLOBAL VARIABLES

C     CONSTANTS

C     ARGUMENTS 
      integer natm,next
      integer a(natm),atm(2,natm),maxk(natm),ncf
      integer cflow(next,1000)
C     LOCAL VARIABLES 
      integer i,j,cfpos,k1,k2,k3,k4,k5,k6,k7,k8
      integer aatm(2,natm)
C     EXTERNAL FUNCTIONS

C     ----------
C     BEGIN CODE
C     ----------
      ncf = 0

      if (natm.eq.3) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do j = 1,2
                     aatm(j,1) = atm(mod(j*k1,3),1)
                     aatm(j,2) = atm(mod(j*k2,3),2)
                     aatm(j,3) = atm(mod(j*k3,3),3)
                  enddo
                  ncf = ncf +1
                  call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
               enddo 
            enddo
         enddo
      elseif (natm.eq.4) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do j = 1,2
                        aatm(j,1) = atm(mod(j*k1,3),1)
                        aatm(j,2) = atm(mod(j*k2,3),2)
                        aatm(j,3) = atm(mod(j*k3,3),3)
                        aatm(j,4) = atm(mod(j*k4,3),4)
                     enddo
                     ncf = ncf +1
                     call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                  enddo 
               enddo
            enddo
         enddo
      elseif (natm.eq.5) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do j = 1,2
                           aatm(j,1) = atm(mod(j*k1,3),1)
                           aatm(j,2) = atm(mod(j*k2,3),2)
                           aatm(j,3) = atm(mod(j*k3,3),3)
                           aatm(j,4) = atm(mod(j*k4,3),4)
                           aatm(j,5) = atm(mod(j*k5,3),5)
                        enddo
                        ncf = ncf +1
                        call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                     enddo
                  enddo 
               enddo
            enddo
         enddo      
      elseif (natm.eq.6) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do k6 = 1,maxk(6)
                           do j = 1,2
                              aatm(j,1) = atm(mod(j*k1,3),1)
                              aatm(j,2) = atm(mod(j*k2,3),2)
                              aatm(j,3) = atm(mod(j*k3,3),3)
                              aatm(j,4) = atm(mod(j*k4,3),4)
                              aatm(j,5) = atm(mod(j*k5,3),5)
                              aatm(j,6) = atm(mod(j*k6,3),6)
                           enddo
                           ncf = ncf +1
                          call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                        enddo
                     enddo
                  enddo 
               enddo
            enddo
         enddo      
      elseif (natm.eq.7) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do k6 = 1,maxk(6)
                           do k7 = 1,maxk(7)
                              do j = 1,2
                                 aatm(j,1) = atm(mod(j*k1,3),1)
                                 aatm(j,2) = atm(mod(j*k2,3),2)
                                 aatm(j,3) = atm(mod(j*k3,3),3)
                                 aatm(j,4) = atm(mod(j*k4,3),4)
                                 aatm(j,5) = atm(mod(j*k5,3),5)
                                 aatm(j,6) = atm(mod(j*k6,3),6)
                                 aatm(j,7) = atm(mod(j*k7,3),7)
                              enddo
                              ncf = ncf +1
                          call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                       enddo
                    enddo
                 enddo
              enddo 
           enddo
        enddo
      enddo      
      elseif (natm.eq.8) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do k6 = 1,maxk(6)
                           do k7 = 1,maxk(7)
                              do k8 = 1,maxk(8)
                                 do j = 1,2
                                    aatm(j,1) = atm(mod(j*k1,3),1)
                                    aatm(j,2) = atm(mod(j*k2,3),2)
                                    aatm(j,3) = atm(mod(j*k3,3),3)
                                    aatm(j,4) = atm(mod(j*k4,3),4)
                                    aatm(j,5) = atm(mod(j*k5,3),5)
                                    aatm(j,6) = atm(mod(j*k6,3),6)
                                    aatm(j,7) = atm(mod(j*k7,3),7)
                                    aatm(j,8) = atm(mod(j*k8,3),8)
                                 enddo
                                 ncf = ncf +1
                          call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo 
               enddo
            enddo
         enddo
      endif

      return
      end
      

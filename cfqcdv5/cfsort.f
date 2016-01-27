      subroutine cfsort(a,np,ng,splitter, b,nabelian)

      implicitnone

      integer np,ng,i,j,k,splitter
      integer a(np),b(ng),nabelian,abelianflag,qline(2,4)

      abelianflag = 0
      nabelian = 0
      do j = 1,np
         if (a(j).lt.splitter) then
            if (abelianflag.ne.1) then
                b(j) = a(j)
            else
               b(j-1) = a(j)
               nabelian = nabelian +1
            endif
         else
            qline(1,1) = a(j)
            qline(2,1) = j
            abelianflag = 1
         endif
      enddo
      
      return
      end

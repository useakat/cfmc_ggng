C   21/01/86 602101130 MEMBER NAME  TOUP     (S)           FORTRAN
************************************************************************
*     ------------------------------
      subroutine tolow(ctxin,ctxout)
*     ------------------------------
*
*     (Purpose)
*        Convert character string with upper case characters
*        into lower case character string.
*
*     (Input)
*        CTXIN ;  Input character string.
*
*     (Output)
*        CTXOUT ; Converted character string.
*
*     (History)
*        2008.08.03 : First version by J.Kanzaki. (converted from toup.f)
*
************************************************************************

      character*(*) ctxin,ctxout
      character*1 char

      character*26 charup  /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      character*26 charlow /'abcdefghijklmnopqrstuvwxyz'/
      integer nchar
      parameter (nchar = 26)

      ltext=len(ctxin)

      do i=1,ltext
         char=ctxin(i:i)
         do j=1,nchar
            if (char.eq.charup(j:j)) then
               char = charlow(j:j)
               goto 100
            endif
         enddo
 100     continue
         ctxout(i:i)=char
      enddo

      return
      end

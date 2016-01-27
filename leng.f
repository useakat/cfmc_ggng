C   15/06/85 506171243 MEMBER NAME  LENG     (S)           FORTRAN
      FUNCTION LENG(CSTR)

      CHARACTER*(*) CSTR
      CHARACTER*1 CBLK/' '/


      LSTR=LEN(CSTR)

      IF(LSTR.GT.0) THEN
         DO 10 I=LSTR,1,-1
         IF(CSTR(I:I).NE.CBLK) THEN
            LENG=I
            RETURN
         END IF
   10    CONTINUE
      END IF

      LENG=0

      RETURN
      END

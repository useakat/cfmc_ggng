C   13/02/85 601202147 MEMBER NAME  COMSEP   (S)           FORTRAN
************************************************************************
*    --------------------------------------------------
      SUBROUTINE COMSEP(CLINE,MAXWRD,NWORD,CWORD,LWORD)
*    --------------------------------------------------
*
*     (Purpose)
*        Separate words in a sentense - STRING - and store them to arrar
*
*     (Input)
*        STRING : Input string which has words.
*        MAXWRD : Upper limit of word# in STRING - you have to define th
*                 of arrary WORD and LWORD more than MAXWRD in calling p
*
*     (Outputs)
*        NWORD  : Total# of words found in STRING.
*        WORD   : Array which has words found.
*        LWORD  : Array which has # of characters in each word
*
*     (Memo)
*        This routine interpret following characters as delimiters.
*            space( ), comma(,), slash(/)
*
*        All spaces at a head of a word will be ignored.
*
*     (Creation Date)
*        1985.02.13 : First version by  K.Amako.
*
************************************************************************

      CHARACTER*(*) CLINE,CWORD
      CHARACTER*1   CHAR,CBLK,CCOM,CSLSH,CAPOS

      DIMENSION CWORD(1),LWORD(1)

      DATA CBLK/' '/,CCOM/','/,CSLSH/'/'/,CAPOS/''''/


***** Initialization
      NWORD=0

      DO 10 I=1,MAXWRD
      CWORD(I)=CBLK
      LWORD(I)=0
   10 CONTINUE

      MAXLEN=LENG(CLINE)
      IF(MAXLEN.LE.0) RETURN

      IFWORD=0
      IFAPOS=0
      IFCONT=0

      DO 50 ICHAR=1,MAXLEN

      CHAR=CLINE(ICHAR:ICHAR)

      IF(CHAR.EQ.CBLK .OR. CHAR.EQ.CCOM) THEN
         IF(IFWORD.EQ.0) GO TO 50
         IF(IFAPOS.EQ.0 .OR. IFCONT.EQ.1) THEN
            LWORD(NWORD)=LW
            IFAPOS=0
            IFWORD=0
         ELSE
            IFCONT=0
            LW=LW+1
            CWORD(NWORD)(LW:LW)=CHAR
         END IF
      ELSE IF(CHAR.EQ.CAPOS) THEN
         IF(IFWORD.EQ.0) THEN
            IFAPOS=1
            IFWORD=1
            LW=0
            NWORD=NWORD+1
            IF(NWORD.GT.MAXWRD) THEN
               NWORD=NWORD-1
               RETURN
            END IF
         ELSE
            IF(IFCONT.EQ.1) THEN
               LW=LW+1
               CWORD(NWORD)(LW:LW)=CHAR
            ELSE
               IFCONT=1
            END IF
         END IF
      ELSE
         IFCONT=0
         IF(IFWORD.EQ.0) THEN
            IFWORD=1
            LW=0
            NWORD=NWORD+1
            IF(NWORD.GT.MAXWRD) THEN
               NWORD=NWORD-1
               RETURN
            END IF
         END IF
         LW=LW+1
         CWORD(NWORD)(LW:LW)=CHAR
      END IF

   50 CONTINUE

      IF(IFWORD.EQ.1) LWORD(NWORD)=LW

      RETURN
      END

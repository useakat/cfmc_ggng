      COMPLEX*16 FUNCTION FLOW1(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. 1.4.1, 2012-06-02
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 1
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: g g > g g g g g QCD=5 QED=0 WEIGHTED=5 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=70)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=7)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=49)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/-1,-1,1,1,1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL VXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL VXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL VXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL VXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL VXXXXX(P(0,IP(7)),ZERO,NHEL(IP(7)),1*IC(IP(7)),W(1,7))
      CALL VVV1_1(W(1,4),W(1,3),GC_4,ZERO, ZERO, W(1,8))
      CALL VVV1_1(W(1,2),W(1,1),GC_4,ZERO, ZERO, W(1,9))
      CALL VVV1_1(W(1,6),W(1,5),GC_4,ZERO, ZERO, W(1,10))
      CALL VVV1_1(W(1,7),W(1,6),GC_4,ZERO, ZERO, W(1,11))
      CALL VVV1_1(W(1,5),W(1,4),GC_4,ZERO, ZERO, W(1,12))
      CALL VVV1_1(W(1,3),W(1,2),GC_4,ZERO, ZERO, W(1,13))
      CALL VVV1_1(W(1,7),W(1,1),GC_4,ZERO, ZERO, W(1,14))
      CALL VVV1_1(W(1,5),W(1,8),GC_4,ZERO, ZERO, W(1,15))
      CALL VVV1_1(W(1,12),W(1,3),GC_4,ZERO, ZERO, W(1,16))
      CALL VVVV1_1(W(1,5),W(1,4),W(1,3),GC_6,ZERO, ZERO, W(1,17))
      CALL VVVV4_1(W(1,5),W(1,4),W(1,3),GC_6,ZERO, ZERO, W(1,18))
      CALL SUMV4(-2.*ONE,W(1,15),-2.*ONE,W(1,16),-2.*ONE,W(1,17),2.
     $ *ONE,W(1,18),W(1,19))
      CALL VVV1_1(W(1,7),W(1,9),GC_4,ZERO, ZERO, W(1,20))
      CALL VVV1_1(W(1,2),W(1,14),GC_4,ZERO, ZERO, W(1,21))
      CALL VVVV1_1(W(1,7),W(1,2),W(1,1),GC_6,ZERO, ZERO, W(1,22))
      CALL VVVV3_1(W(1,7),W(1,2),W(1,1),GC_6,ZERO, ZERO, W(1,23))
      CALL SUMV4(2.*ONE,W(1,20),2.*ONE,W(1,21),2.*ONE,W(1,22),2.
     $ *ONE,W(1,23),W(1,24))
      CALL VVV1_1(W(1,7),W(1,10),GC_4,ZERO, ZERO, W(1,25))
      CALL VVV1_1(W(1,11),W(1,5),GC_4,ZERO, ZERO, W(1,26))
      CALL VVVV1_1(W(1,7),W(1,6),W(1,5),GC_6,ZERO, ZERO, W(1,27))
      CALL VVVV4_1(W(1,7),W(1,6),W(1,5),GC_6,ZERO, ZERO, W(1,28))
      CALL SUMV4(-2.*ONE,W(1,25),-2.*ONE,W(1,26),-2.*ONE,W(1,27),2.
     $ *ONE,W(1,28),W(1,29))
      CALL VVV1_1(W(1,6),W(1,12),GC_4,ZERO, ZERO, W(1,30))
      CALL VVV1_1(W(1,10),W(1,4),GC_4,ZERO, ZERO, W(1,31))
      CALL VVVV1_1(W(1,6),W(1,5),W(1,4),GC_6,ZERO, ZERO, W(1,32))
      CALL VVVV4_1(W(1,6),W(1,5),W(1,4),GC_6,ZERO, ZERO, W(1,33))
      CALL SUMV4(-2.*ONE,W(1,30),-2.*ONE,W(1,31),-2.*ONE,W(1,32),2.
     $ *ONE,W(1,33),W(1,34))
      CALL VVV1_1(W(1,3),W(1,9),GC_4,ZERO, ZERO, W(1,35))
      CALL VVV1_1(W(1,13),W(1,1),GC_4,ZERO, ZERO, W(1,36))
      CALL VVVV1_1(W(1,3),W(1,2),W(1,1),GC_6,ZERO, ZERO, W(1,37))
      CALL VVVV4_1(W(1,3),W(1,2),W(1,1),GC_6,ZERO, ZERO, W(1,38))
      CALL SUMV4(-2.*ONE,W(1,35),-2.*ONE,W(1,36),-2.*ONE,W(1,37),2.
     $ *ONE,W(1,38),W(1,39))
      CALL VVV1_1(W(1,4),W(1,13),GC_4,ZERO, ZERO, W(1,40))
      CALL VVV1_1(W(1,8),W(1,2),GC_4,ZERO, ZERO, W(1,41))
      CALL VVVV1_1(W(1,4),W(1,3),W(1,2),GC_6,ZERO, ZERO, W(1,42))
      CALL VVVV4_1(W(1,4),W(1,3),W(1,2),GC_6,ZERO, ZERO, W(1,43))
      CALL SUMV4(-2.*ONE,W(1,40),-2.*ONE,W(1,41),-2.*ONE,W(1,42),2.
     $ *ONE,W(1,43),W(1,44))
      CALL VVV1_1(W(1,6),W(1,14),GC_4,ZERO, ZERO, W(1,45))
      CALL VVV1_1(W(1,11),W(1,1),GC_4,ZERO, ZERO, W(1,46))
      CALL VVVV3_1(W(1,7),W(1,6),W(1,1),GC_6,ZERO, ZERO, W(1,47))
      CALL VVVV4_1(W(1,7),W(1,6),W(1,1),GC_6,ZERO, ZERO, W(1,48))
      CALL SUMV4(-2.*ONE,W(1,45),2.*ONE,W(1,46),-2.*ONE,W(1,47),
     $ -2.*ONE,W(1,48),W(1,49))
C     Amplitude(s) for diagram number 1
      CALL VVV1_0(W(1,6),W(1,19),W(1,24),GC_4,AMP(1))
C     Amplitude(s) for diagram number 2
      CALL VVVV1_0(W(1,6),W(1,5),W(1,8),W(1,24),GC_6,AMP(2))
      CALL VVVV4_0(W(1,6),W(1,5),W(1,8),W(1,24),GC_6,AMP(3))
C     Amplitude(s) for diagram number 3
      CALL VVVV1_0(W(1,7),W(1,6),W(1,19),W(1,9),GC_6,AMP(4))
      CALL VVVV4_0(W(1,7),W(1,6),W(1,19),W(1,9),GC_6,AMP(5))
C     Amplitude(s) for diagram number 4
      CALL VVVV1_0(W(1,7),W(1,10),W(1,8),W(1,9),GC_6,AMP(6))
      CALL VVVV4_0(W(1,7),W(1,10),W(1,8),W(1,9),GC_6,AMP(7))
C     Amplitude(s) for diagram number 5
      CALL VVV1_0(W(1,29),W(1,8),W(1,9),GC_4,AMP(8))
C     Amplitude(s) for diagram number 6
      CALL VVV1_0(W(1,10),W(1,8),W(1,24),GC_4,AMP(9))
C     Amplitude(s) for diagram number 7
      CALL VVVV1_0(W(1,11),W(1,5),W(1,8),W(1,9),GC_6,AMP(10))
      CALL VVVV4_0(W(1,11),W(1,5),W(1,8),W(1,9),GC_6,AMP(11))
C     Amplitude(s) for diagram number 8
      CALL VVV1_0(W(1,11),W(1,19),W(1,9),GC_4,AMP(12))
C     Amplitude(s) for diagram number 9
      CALL VVV1_0(W(1,7),W(1,34),W(1,39),GC_4,AMP(13))
C     Amplitude(s) for diagram number 10
      CALL VVVV1_0(W(1,7),W(1,6),W(1,12),W(1,39),GC_6,AMP(14))
      CALL VVVV4_0(W(1,7),W(1,6),W(1,12),W(1,39),GC_6,AMP(15))
C     Amplitude(s) for diagram number 11
      CALL VVV1_0(W(1,34),W(1,3),W(1,24),GC_4,AMP(16))
C     Amplitude(s) for diagram number 12
      CALL VVVV1_0(W(1,6),W(1,12),W(1,3),W(1,24),GC_6,AMP(17))
      CALL VVVV4_0(W(1,6),W(1,12),W(1,3),W(1,24),GC_6,AMP(18))
C     Amplitude(s) for diagram number 13
      CALL VVVV1_0(W(1,7),W(1,34),W(1,3),W(1,9),GC_6,AMP(19))
      CALL VVVV4_0(W(1,7),W(1,34),W(1,3),W(1,9),GC_6,AMP(20))
C     Amplitude(s) for diagram number 14
      CALL VVVV1_0(W(1,11),W(1,12),W(1,3),W(1,9),GC_6,AMP(21))
      CALL VVVV4_0(W(1,11),W(1,12),W(1,3),W(1,9),GC_6,AMP(22))
C     Amplitude(s) for diagram number 15
      CALL VVV1_0(W(1,11),W(1,12),W(1,39),GC_4,AMP(23))
C     Amplitude(s) for diagram number 16
      CALL VVV1_0(W(1,29),W(1,4),W(1,39),GC_4,AMP(24))
C     Amplitude(s) for diagram number 17
      CALL VVVV1_0(W(1,7),W(1,10),W(1,4),W(1,39),GC_6,AMP(25))
      CALL VVVV4_0(W(1,7),W(1,10),W(1,4),W(1,39),GC_6,AMP(26))
C     Amplitude(s) for diagram number 18
      CALL VVVV1_0(W(1,10),W(1,4),W(1,3),W(1,24),GC_6,AMP(27))
      CALL VVVV4_0(W(1,10),W(1,4),W(1,3),W(1,24),GC_6,AMP(28))
C     Amplitude(s) for diagram number 19
      CALL VVVV1_0(W(1,29),W(1,4),W(1,3),W(1,9),GC_6,AMP(29))
      CALL VVVV4_0(W(1,29),W(1,4),W(1,3),W(1,9),GC_6,AMP(30))
C     Amplitude(s) for diagram number 20
      CALL VVVV1_0(W(1,11),W(1,5),W(1,4),W(1,39),GC_6,AMP(31))
      CALL VVVV4_0(W(1,11),W(1,5),W(1,4),W(1,39),GC_6,AMP(32))
C     Amplitude(s) for diagram number 21
      CALL VVV1_0(W(1,5),W(1,44),W(1,49),GC_4,AMP(33))
C     Amplitude(s) for diagram number 22
      CALL VVVV1_0(W(1,5),W(1,4),W(1,13),W(1,49),GC_6,AMP(34))
      CALL VVVV4_0(W(1,5),W(1,4),W(1,13),W(1,49),GC_6,AMP(35))
C     Amplitude(s) for diagram number 23
      CALL VVVV1_0(W(1,6),W(1,5),W(1,44),W(1,14),GC_6,AMP(36))
      CALL VVVV4_0(W(1,6),W(1,5),W(1,44),W(1,14),GC_6,AMP(37))
C     Amplitude(s) for diagram number 24
      CALL VVVV1_0(W(1,6),W(1,12),W(1,13),W(1,14),GC_6,AMP(38))
      CALL VVVV4_0(W(1,6),W(1,12),W(1,13),W(1,14),GC_6,AMP(39))
C     Amplitude(s) for diagram number 25
      CALL VVV1_0(W(1,34),W(1,13),W(1,14),GC_4,AMP(40))
C     Amplitude(s) for diagram number 26
      CALL VVV1_0(W(1,12),W(1,13),W(1,49),GC_4,AMP(41))
C     Amplitude(s) for diagram number 27
      CALL VVVV1_0(W(1,10),W(1,4),W(1,13),W(1,14),GC_6,AMP(42))
      CALL VVVV4_0(W(1,10),W(1,4),W(1,13),W(1,14),GC_6,AMP(43))
C     Amplitude(s) for diagram number 28
      CALL VVV1_0(W(1,10),W(1,44),W(1,14),GC_4,AMP(44))
C     Amplitude(s) for diagram number 29
      CALL VVV1_0(W(1,19),W(1,2),W(1,49),GC_4,AMP(45))
C     Amplitude(s) for diagram number 30
      CALL VVVV1_0(W(1,5),W(1,8),W(1,2),W(1,49),GC_6,AMP(46))
      CALL VVVV4_0(W(1,5),W(1,8),W(1,2),W(1,49),GC_6,AMP(47))
C     Amplitude(s) for diagram number 31
      CALL VVVV1_0(W(1,6),W(1,19),W(1,2),W(1,14),GC_6,AMP(48))
      CALL VVVV4_0(W(1,6),W(1,19),W(1,2),W(1,14),GC_6,AMP(49))
C     Amplitude(s) for diagram number 32
      CALL VVVV1_0(W(1,10),W(1,8),W(1,2),W(1,14),GC_6,AMP(50))
      CALL VVVV4_0(W(1,10),W(1,8),W(1,2),W(1,14),GC_6,AMP(51))
C     Amplitude(s) for diagram number 33
      CALL VVVV1_0(W(1,12),W(1,3),W(1,2),W(1,49),GC_6,AMP(52))
      CALL VVVV4_0(W(1,12),W(1,3),W(1,2),W(1,49),GC_6,AMP(53))
C     Amplitude(s) for diagram number 34
      CALL VVVV1_0(W(1,34),W(1,3),W(1,2),W(1,14),GC_6,AMP(54))
      CALL VVVV4_0(W(1,34),W(1,3),W(1,2),W(1,14),GC_6,AMP(55))
C     Amplitude(s) for diagram number 35
      CALL VVVV1_0(W(1,7),W(1,34),W(1,13),W(1,1),GC_6,AMP(56))
      CALL VVVV4_0(W(1,7),W(1,34),W(1,13),W(1,1),GC_6,AMP(57))
C     Amplitude(s) for diagram number 36
      CALL VVVV1_0(W(1,11),W(1,12),W(1,13),W(1,1),GC_6,AMP(58))
      CALL VVVV4_0(W(1,11),W(1,12),W(1,13),W(1,1),GC_6,AMP(59))
C     Amplitude(s) for diagram number 37
      CALL VVV1_0(W(1,29),W(1,44),W(1,1),GC_4,AMP(60))
C     Amplitude(s) for diagram number 38
      CALL VVVV1_0(W(1,29),W(1,4),W(1,13),W(1,1),GC_6,AMP(61))
      CALL VVVV4_0(W(1,29),W(1,4),W(1,13),W(1,1),GC_6,AMP(62))
C     Amplitude(s) for diagram number 39
      CALL VVVV1_0(W(1,7),W(1,10),W(1,44),W(1,1),GC_6,AMP(63))
      CALL VVVV4_0(W(1,7),W(1,10),W(1,44),W(1,1),GC_6,AMP(64))
C     Amplitude(s) for diagram number 40
      CALL VVVV1_0(W(1,11),W(1,5),W(1,44),W(1,1),GC_6,AMP(65))
      CALL VVVV4_0(W(1,11),W(1,5),W(1,44),W(1,1),GC_6,AMP(66))
C     Amplitude(s) for diagram number 41
      CALL VVVV1_0(W(1,29),W(1,8),W(1,2),W(1,1),GC_6,AMP(67))
      CALL VVVV4_0(W(1,29),W(1,8),W(1,2),W(1,1),GC_6,AMP(68))
C     Amplitude(s) for diagram number 42
      CALL VVVV1_0(W(1,11),W(1,19),W(1,2),W(1,1),GC_6,AMP(69))
      CALL VVVV4_0(W(1,11),W(1,19),W(1,2),W(1,1),GC_6,AMP(70))
      JAMP(1)=+1./2.*IMAG1*AMP(1)-IMAG1*AMP(2)+IMAG1*AMP(3)-IMAG1
     $ *AMP(4)+IMAG1*AMP(5)+2*IMAG1*AMP(6)-2*IMAG1*AMP(7)-IMAG1*AMP(8)
     $ -IMAG1*AMP(9)+2*IMAG1*AMP(10)-2*IMAG1*AMP(11)-IMAG1*AMP(12)
     $ +1./2.*IMAG1*AMP(13)-IMAG1*AMP(14)+IMAG1*AMP(15)+1./2.*IMAG1
     $ *AMP(16)-IMAG1*AMP(17)+IMAG1*AMP(18)-IMAG1*AMP(19)+IMAG1
     $ *AMP(20)+2*IMAG1*AMP(21)-2*IMAG1*AMP(22)-IMAG1*AMP(23)
     $ +1./2.*IMAG1*AMP(24)-IMAG1*AMP(25)+IMAG1*AMP(26)-IMAG1*AMP(27)
     $ +IMAG1*AMP(28)-IMAG1*AMP(29)+IMAG1*AMP(30)-IMAG1*AMP(31)
     $ +IMAG1*AMP(32)+1./2.*IMAG1*AMP(33)-IMAG1*AMP(34)+IMAG1*AMP(35)
     $ +IMAG1*AMP(36)-IMAG1*AMP(37)-2*IMAG1*AMP(38)+2*IMAG1*AMP(39)
     $ +IMAG1*AMP(40)-IMAG1*AMP(41)-2*IMAG1*AMP(42)+2*IMAG1*AMP(43)
     $ +IMAG1*AMP(44)+1./2.*IMAG1*AMP(45)-IMAG1*AMP(46)+IMAG1*AMP(47)
     $ +IMAG1*AMP(48)-IMAG1*AMP(49)-2*IMAG1*AMP(50)+2*IMAG1*AMP(51)
     $ -IMAG1*AMP(52)+IMAG1*AMP(53)+IMAG1*AMP(54)-IMAG1*AMP(55)
     $ -IMAG1*AMP(56)+IMAG1*AMP(57)+2*IMAG1*AMP(58)-2*IMAG1*AMP(59)
     $ +1./2.*IMAG1*AMP(60)-IMAG1*AMP(61)+IMAG1*AMP(62)-IMAG1*AMP(63)
     $ +IMAG1*AMP(64)-IMAG1*AMP(65)+IMAG1*AMP(66)-IMAG1*AMP(67)
     $ +IMAG1*AMP(68)-IMAG1*AMP(69)+IMAG1*AMP(70)

      FLOW1=JAMP(1)
      RETURN
      END

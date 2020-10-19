**********************************************************************
    
!    Copyright 2017, I.J. Thompson
!     
!    This file is part of FRESCO.
!
!    FRESCO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!     
!    FRESCO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!     
!    You should have received a copy of the GNU General Public License
!    along with FRESCO. If not, see <http://www.gnu.org/licenses/>.
!     
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

*****FRXX6**************************************************************
*****TWONN*************************************************************
      SUBROUTINE TWONN(KN1,KN2,HDN,MR,RNN,ISC,IPC,NLN,NN,
     &  QNF, EGS,KIND,IN, IAK,IAMIN,IAMAX,IB,JEX,BAND,MXP,MXX,ICR,ICP,
     &  COPY,KNZR,LMIN,LMAX,SMIN,SMAX,J12MIN,J12MAX,T,D0,RIN,DM,CRM,
     &  NK,TNT,COEF,NFL,cxwf,FORML,MAXNLN,MSP,FORMR,RMI,FMSCAL,
     &  NN2,NNU,BE,FK,NRAD,EP2,COM,NK1)
	use io
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 J12,J12MIN,J12MAX,COEF(NK1),J,JA,JB,D0(MSP,2)
      INTEGER QNF(19,MSP),T,TNT(4,NK1),BAND(2,MXP,MXX),
     &      COPY(2,MXP,MXX,2)
      REAL*8 FORML(MAXNLN,MSP,2),FORMR(NRAD),WFC,FFR4
      REAL*8 FNN(NLN,NN2),REL(NN2),WREL(NN2)
     &  , BE(MSP,3),XG(6),WG(6),JEX(6,MXP,MXX),JCORE,JCOMB,JCOM
      REAL*8 FN4(NLN,NN2)
	real*8,allocatable:: QERN(:,:,:,:)
      CHARACTER*8 COM
      LOGICAL FRAC,FAIL3,PRES,GAUS,cxwf
      DATA NWW,XG(4),XG(5),XG(6),WG(4),WG(5),WG(6)/3,
     1   .2386191861D0,.6612093865D0,.9324695142D0,
     2   .4679139346D0,.3607615730D0,.1713244924D0/
      FRAC(X) = ABS(X-NINT(X)).GT.1D-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
      Z = 0.0
      NW=2*NWW
      DO 1 N=1,NWW
      JJ=NW-N+1
      XG(N)=-XG(JJ)
  1   WG(N)=WG(JJ)
      IF(NFL.NE.0) REWIND ABS(NFL)
      KMAX = 0
      IF(NFL.LE.0) THEN
      DO 3 IK=1,NK
       IF(TNT(1,IK).LT.0.OR.TNT(3,IK).LT.0) GO TO 3
      DO 2 KNA=TNT(1,IK),MSP
      DO 2 KNB=TNT(2,IK),MSP
      IF(QNF(1,KNA).NE.TNT(1,IK) .OR. QNF(1,KNB).NE.TNT(2,IK))GOTO 2
      KMAX = MAX(KMAX, QNF(9,KNA)+QNF(9,KNB) )
2     CONTINUE
3     CONTINUE
      KMAX = KMAX + SMAX + J12MAX + LMAX + 1.1
      ENDIF
      KMAX1 = KMAX + 1
	allocate (QERN(NN2,NLN,KMAX1,NK))
      GAUS = ISC .LE. 0
      RM = DM * CRM / (DM + CRM)
      RED = RM/DM
      DN = DM*0.5
C
C    TRANSFORM 'NK' PAIRS OF S.P. WFNTNS IN [R1,R2],J12 COORDINATES
C                             TO [RR,R] COORDINATES (C.M. & RELATIVE)
C    USE  R1 = RR + R/2
C      &  R2 = RR - R/2      SO
      A = 1.0
      B = 0.5
      P = 1.0
      Q =-0.5

      if(TNT(4,1).eq.0) then
C				NOW USER JACOBIAN COORDINATES::!!!!!!!
C    USE  R1 = RR + R/2
C      &  R2 = P*RR + Q*R      SO
      A = 1.0
      B = 0.5
      P = CRM/(CRM+DN)
      Q =-0.5 * (DM+CRM)/(DN+CRM)
      endif
C
      RINTP = HDN * MR
           IC3 =(NN +NW-1 ) / NW
            CI = (RNN-RMI)*NW/NN
            C1 = 0.5 * CI
            C2 = C1 + RMI
               N = 0
            DO 14 K=1,IC3
            DO 13 JJ=1,NW
               N = N + 1
                  IF(N.GT.NN) GO TO 13
             IF(.NOT.GAUS) THEN
               REL(N) = C1 * (JJ-3.5)/3.0 + C2
               WREL(N) = SQRT(C1 / 3.0)
             ELSE
               REL(N) = C1 * XG(JJ) + C2
               WREL(N) = SQRT(C1 * WG(JJ))
             ENDIF
               REL(N+NN) = REL(N)
               WREL(N+NN) = 0.0
13          CONTINUE
14          C2 = C2 + CI

         IF(NFL.LT.0) WRITE(ABS(NFL),4) GAUS,IAMIN,IAMAX,J12MIN,J12MAX,
     &                 SMIN,SMAX,LMIN,LMAX,NLN,NN2,RINTP,RNN,RMI,RED
4         FORMAT(L4,2I4,4F5.1,4I4,4F8.4)
         IF(NFL.GT.0) READ(NFL,4) GAUS,IAMIN,IAMAX,J12MIN,J12MAX,
     &                 SMIN,SMAX,LMIN,LMAX
C
      IF(NFL.LE.0)
     &CALL QNN(QERN,NLN,NN,KMAX1,NK,TNT,A,B,P,Q,NN2,FORML,
     &         NNU,MAXNLN,RINTP,RIN,MSP,QNF,IPC,XG,WG,NW,REL)
      EPS = 1D-14
      WN = 0.0
      RMS = 0.0
      RMC = 0.0
      WM = 0.0
      KN = KN1
       JCOM = JEX(IN,ICP,MAX(IB,1))
      NJ12=NINT(J12MAX-J12MIN)*2
      NS=NINT(SMAX-SMIN)*2
      DO 100 IA=IAMIN,IAMAX
         JCORE = JEX(IN,ICR,IA)
         IF(IPC.GE.5) WRITE(KO,*) IN,ICP,IB,JCOM,ICR,IA,JCORE
!      DO 100 J12=J12MIN,J12MAX,0.5
      DO 100 IJ12=0,NJ12
      J12=J12MIN+IJ12*0.5
      DO 100 LL=0,NINT(SMAX+J12+LMAX)
      DO 100 L=LMIN,LMAX
!      DO 100 S=SMIN,SMAX,0.5
      DO 100 IS=0,NS
      S=SMIN+IS*0.5
!      DO 100 J=ABS(L-S),L+S
      NJ=NINT(L+S-ABS(L-S))
      DO 100 IJ=0,NJ
      J=ABS(L-S)+IJ
c         IF(IPC.GE.5) WRITE(KO,*) 'J12,LL,L,S,J',J12,LL,L,S,J
         IF(FAIL3(LL+Z,J,J12)) GO TO 100
         IF(IAK.GE.1  .AND. (COPY(IN,ICR,IA,1).NE.0 .OR.
     & (-1)**(LL+L).NE.SIGN(1,BAND(IN,ICP,IB)*BAND(IN,ICR,IA))))GOTO 100
      IF(KIND.EQ.9.AND.FAIL3(J12,JCORE,JCOM) ) GOTO 100
      IF(MOD(NINT(L + S + T),2).eq.0) go to 100
      NMAX = MIN(NN,KN2-KN+1)
      IF(KNZR.GT.0) NMAX = MIN(1,NMAX)
      DO 17 N=1,NMAX
      QNF(3,KN+N-1) = IA
      QNF(9,KN+N-1) = LL
      QNF(10,KN+N-1) = NINT(2.*J)
      QNF(11,KN+N-1) = NINT(2.*J12)
      QNF(12,KN+N-1) = N
       IF(KNZR.NE.0) QNF(12,KN+N-1) = 0
      QNF(13,KN+N-1) = L
      QNF(14,KN+N-1) = NINT(2*S)
      QNF(1,KN+N-1) = 0
      BE(KN+N-1,1) = EGS
      DO 17 I=2,3
      BE(KN+N-1,I) = 0.0
17    D0(KN+N-1,I-1)=0.0
         DO 18 N=1,NN2
         DO 18 I=1,NLN
18       FNN(I,N)=0.
         IF(KNZR.EQ.0 .OR.  KN.GT.KN2) GO TO 191
         DO 19 JJ=1,2
         DO 19 N=1,NMAX
         DO 19 I=1,NLN
19       FORML(I,KN+N-1,JJ)=0.
191   PRES = .FALSE.
      IF(NFL.LE.0) THEN
      R0 = SQRT((2.*LL+1.)*(2.*L+1.)*(2.*S+1.)*(2.*J+1.))
      DO 802 IK=1,NK
      IF(TNT(1,IK).eq.-1) THEN
       CALL EXTERN1(TNT(1,IK),TNT(2,IK),COEF(IK),LL,L,S,J,J12,JCORE,JCOM
     &             ,PRES,IPC,FNN,NLN,NN,NN2,RINTP,REL,EGS)
         BE(KN:KN+NMAX-1,1) = EGS
      ELSE IF(TNT(1,IK).le.-2) THEN
       CALL EXTERN2(TNT(1,IK),TNT(2,IK),COEF(IK),LL,L,S,J,J12,JCORE,JCOM
     &     ,KIND,PRES,IPC,FNN,NLN,NN,NN2,RINTP,REL,WREL,max(NN,NLN),EGS)
         BE(KN:KN+NMAX-1,1) = EGS
      ELSE
      DO 801 KNA=TNT(1,IK),MSP
         IF(QNF(1,KNA).NE.TNT(1,IK)) GO TO 801
      DO 80 KNB=TNT(2,IK),MSP
         IF(QNF(1,KNB).NE.TNT(2,IK))GOTO 80
         KINDA = QNF(7,KNA)
         KINDB = QNF(7,KNB)
         IF(QNF(8,KNA).EQ.0 .OR. QNF(8,KNB).EQ.0) GO TO 80
         LA = QNF(9,KNA)
         LB = QNF(9,KNB)
         JA = QNF(11,KNA) * 0.5
         JB = QNF(11,KNB) * 0.5
         SA = QNF(10,KNA) * 0.5
         SB = QNF(10,KNB) * 0.5
      IF(TNT(3,IK).GE.0 .AND.KINDA.NE.1.AND.KINDB.NE.1) THEN
          IF(FAIL3(S,SA,SB).OR.MOD(LL+L+LA+LB,2).NE.0) GO TO  80
      R2 = COEF(IK)
      R3 = (1 - (-1)**NINT(L + S + T))/SQRT(2.)
      IF(ABS(R2*R3).LT.EPS .AND. IPC.LT.8) GO TO 80
      IF(LA.EQ.LB .AND. 1E-5.GT.ABS(SA-SB)+ABS(JA-JB)+
     &   ABS(QNF(8,KNA)-QNF(8,KNB)) ) R3 = R3 / SQRT(2.)
      R1 = R0 * SQRT((2.*JA+1.)*(2.*JB+1.)*(2.*LA+1.)*(2.*LB+1.))
      DO 70 NA=0,LA
      DO 70 NB=0,LB
         LAMNA = LA - NA
         LBMNB = LB - NB
      R4 = EXP( (FACT(2*LA+1 +1)-FACT(2*NA +1)-FACT(2*LAMNA+1 +1)
     &        + FACT(2*LB+1  +1)-FACT(2*NB +1)-FACT(2*LBMNB+1 +1)) * .5)
      R11 = R1 * SQRT((2.*LAMNA+1.)*(2.*LBMNB+1.))
C
      DO 70 K=0,KMAX1
         K1 = K + 1
      R8 = 0.0
      DO 30 LAM=ABS(NINT(S-J12)),NINT(S+J12)
         IF(IPC.LT.8.AND.FAIL3(LAM+Z,LA+Z,LB+Z)) GO TO 30
      R51 = (2.*LAM+1.) * (-1)**(L+LL-LAM)
      R52 = RACAH(LL+Z,L+Z,J12,S,LAM+Z,J)
      R53 = WIG9J(LA +Z,   LB +Z,   LAM+Z,
     &            SA,      SB,      S,
     &            JA,      JB,      J12)
      R5 = R51 * R52 * R53
      DO 20 LAMA=ABS(LL-K),LL+K
      DO 20 LAMB=ABS(L-K),L+K
         IF(IPC.LT.8.AND.FAIL3(LAMA+Z, LAMB+Z, LAM+Z)) GO TO 20
      R61 = (2.*LAMA+1.) * (2.*LAMB+1.) * (-1)**(LA+LB+LL+LAMB)
      R62 = WIG3J(LAMNA+Z,NB+Z,LAMA+Z, Z,Z,Z)
      R63 = WIG3J(LBMNB+Z,NA+Z,LAMB+Z, Z,Z,Z)
      IF(ABS(R62*R63).LT.EPS .AND. IPC.LT.8) GO TO 20
      R64 = WIG9J(LAMNA+Z, NB+Z,    LAMA+Z,
     &             NA+Z,    LBMNB+Z, LAMB+Z,
     &             LA+Z,    LB+Z,    LAM+Z)
      R6 = R61 * R62 * R63 * R64
      R7  = (2*K+1.) * WIG3J(LAMA+Z,LL+Z,K+Z,Z,Z,Z)
     &               * WIG3J(LAMB+Z,L +Z,K+Z,Z,Z,Z)
     &               * RACAH(LAMA+Z,LL+Z,LAMB+Z,L+Z,K+Z,LAM+Z)
      R57 = R5 * R6 * R7
      R8 = R8 + R57
      IF(IPC.GE.5) WRITE(KO,*) LL,L,S,J,IK,KNA,KNB,LA,LB,NA,NB,K
     & ,LAMNA,LBMNB,LAM,LAMA,LAMB,R5,R6,R7,R57,R8
      IF(IPC.GE.6) WRITE(KO,*) R51,R52,R53,R61,R62,R63,R64,R7
20    CONTINUE
30    CONTINUE
      R9 = R11 * R2 * R3 * R4 * R8
      IF(ABS(R9).LT.EPS .AND. IPC.LT.8) GO TO 70
        PRES = .TRUE.
      TC = R9 * A**LAMNA * B**NA * P**NB * Q**LBMNB
      IF(IPC.GE.4) WRITE(KO,32) LL,L,S,J,IK,KNA,KNB,LA,LB,NA,NB,K
     &       ,R2,R3,R4*R11*R8,R9,TC
32    FORMAT('0',2I5,2F5.1,8I4,3F10.5,2F12.5)
C     IF(K.GT.KMAX.AND.IPC.LT.8) STOP 'KMAX'
      IF(K.GT.KMAX.AND.IPC.LT.8) CALL ABEND(32)
      IF(ABS(R9).LT.EPS) GO TO 70
        NRF = LBMNB + NA
        DO 50 I=1,NLN
        RR = (I-1)*RINTP
        TT = TC
        IF(I.GT.1) TT = TC * RR**(LAMNA + NB)
        DO 50 JJ=1,NN2
  50    FNN(I,JJ) = FNN(I,JJ) + TT*QERN(JJ,I,K1,IK)*REL(JJ)**NRF
70    CONTINUE
      ELSE
C            INSERT PRODUCT OF CLUSTER WAVE FUNCTIONS
C
      IF(LA.NE.LL.OR.ABS(JA-J12)+ABS(SA-J).GT.1E-5.OR.LB.NE.L) GO TO 80
      IF(KINDB.EQ.0) THEN
         IF(ABS(S-SB) + ABS(SA-JB) .GT.1E-5) GO TO 80
      ELSE IF(KINDB.EQ.1) THEN
           JCOMB = JEX(QNF(6,KNB),QNF(4,KNB),QNF(5,KNB))
         IF(ABS(S-JB) + ABS(SA-JCOMB) .GT. 1E-5) GOTO 80
      ENDIF
      IF(IPC.GE.3) WRITE(KO,32) LL,L,S,J,IK,KNA,KNB,LA,LB
      DO 75 JJ=1,NN
            R = REL(JJ)
        TT = COEF(IK) * FFR4(R*RIN,FORML(1,KNB,1),NLN) * R**LB
              R1 = 1.0
           DO 75 I=1,NLN
              IF(I.GT.1) R1 = ((I-1)*RINTP)**LA
           FNN(I,JJ)    = FNN(I,JJ)    +TT*R1*FORML(I,KNA,1)
75         FNN(I,JJ+NN) = FNN(I,JJ+NN) +TT*R1*FORML(I,KNA,2)
      PRES = .TRUE.
      ENDIF
80    CONTINUE
801   CONTINUE
      ENDIF
802   CONTINUE
      ENDIF
      	FK = SQRT(ABS(EGS) * FMSCAL * RM)
         IF(NFL.LT.0) WRITE(ABS(NFL),4) PRES,LL,L,S,J,J12
         IF(NFL.GT.0) READ(NFL,4)       PRES
      IF(.NOT.PRES) GO TO 100
        IF(NFL.LT.0) WRITE(ABS(NFL),804) FNN
C             READ IN FNN FROM FILE NFL > 0
        IF(NFL.GT.0) READ(NFL,804) FNN
 804    FORMAT(1P,6E12.4)
           WFC = 0.0
            CI = 0.0
            DO 81 JJ=1,NN2
               R = REL(JJ)
                TT = WREL(MOD(JJ-1,NN)+1)**2
            DO 81 I=1,NLN
                 RR = (I-1)*RINTP
            WF        = FNN(I,JJ) * RR * R
            FN4(I,JJ) = WF
            IF(KNZR.NE.0 .AND. KN.LE.KN2) THEN
             DO 805 KNZ=KNZR,MSP
                IF(QNF(1,KNZ).NE.KNZR) GO TO 805
             IF(L.NE.QNF(9,KNZ)) GO TO 805
             WFC = FFR4(R*RIN,FORML(1,KNZ,1),NLN) * R**(QNF(9,KNZ)+1)
             FORML(I,KN,1+(JJ-1)/NN) = FORML(I,KN,1+(JJ-1)/NN)
     &                  + WFC * WF * TT
805          CONTINUE
            ENDIF
             IF(JJ.GT.NN) GOTO 81
               CI = CI + WF**2 * TT*RINTP
               RMS= RMS+WF**2 * ((RED*RR)**2+(0.5*R)**2) * TT*RINTP
               RMC= RMC + WF**2 *  ((1-RED)*RR)**2       * TT*RINTP
81           CONTINUE
        WM = WM + CI
        IF(CI.LT.WM * EP2*0.01 - EPS) GO TO 100
         IF(IPC.GE.2.AND.KNZR.LE.0) THEN
            WRITE(KO,82) LL,L,S,J,CI
82          FORMAT(/' TWO-PARTICLE FORM FACTOR for LL,l,s,j =',2I4,2F4.1
     &             ,'  : Sq. Norm =',F10.6)
            CALL DISPLY(FN4,NLN,NN,NLN,SCALE)
            IF(ISC.NE.0) THEN
               DO 825 JJ=1,NN
825            WRITE(KO,826) REL(JJ),RINTP,(FN4(I,JJ),I=1,NLN)
826            FORMAT(' Wfn at r =',F6.2,' for RR from 0 in steps of',
     &                F6.3,' fm. :'/(1X,12F10.5))
               ENDIF
            IF(SCALE.EQ.Z) GO TO 100
            IF(.NOT.IPC.GE.5) GO TO 84
            WRITE(KO,83)
83          FORMAT('0Form of Interaction Potential * Wave Function :'/)
            CALL DISPLY(FN4(1,NN+1),NLN,NN,NLN,SCALE)
        ENDIF
84    DO 90 I=1,NLN
         RR = (I-1)*RINTP
         R1 = 1.0
         IF(I.GT.1) R1 = RR**(-LL)
         DO 90 JJ=1,NMAX
            KM = KN+JJ-1
       IF(KNZR.EQ.0) THEN
         FORML(I,KM,1) = FNN(I,JJ) * WREL(JJ) * R1 * REL(JJ)
         FORML(I,KM,2) = FNN(I,JJ+NN)*WREL(JJ)* R1 * REL(JJ)
       ELSE
         FN4(I,JJ) = FORML(I,KM,1)
         FN4(I,JJ+NN)=FORML(I,KM,2)
         FORML(I,KM,1) = FORML(I,KM,1) * R1 / MAX(RR,EPS)
         FORML(I,KM,2) = FORML(I,KM,2) * R1 / MAX(RR,EPS)
       ENDIF
         BE(KM,2)= BE(KM,2)+ FN4(I,JJ)**2
         BE(KM,3)= BE(KM,3)+(FN4(I,JJ)*RR)**2
         D0(KM,1) = D0(KM,1) + FN4(I,JJ+NN)*RR/R1
         D0(KM,2) = D0(KM,2) + FN4(I,JJ+NN)*SINH(FK*RR)/(R1*FK)
90       CONTINUE
      TT = 1.0 * RINTP * SQRT(4*PI)
      DO 91 N=1,LL
91    TT = TT / (2.*N+1.)
      DO 96 N=1,NMAX
         KM = KN + N - 1
      D0(KM,1) = D0(KM,1) * TT
      D0(KM,2) = D0(KM,2) * TT
      BE(KM,2) = SQRT(BE(KM,2) * RINTP)
      BE(KM,3) = SQRT(BE(KM,3) * RINTP) / (BE(KM,2)+1e-20)
         IF(KNZR.EQ.0) WN = WN + (BE(KM,2) * WREL(N))**2
         IF(KNZR.NE.0) WN = WN + BE(KM,2)**2
      QNF(1,KM) = KN1
      QNF(8,KM) = 1
      R1 = 1E-10 * MAX(BE(KM,2),1E-20+Z)
      DO 92 I=3,NINT(NLN - REL(NN)*MAX(ABS(B),ABS(Q))/RINTP)
92    IF(FN4(I,N)*FN4(I-1,N).LT.-R1) QNF(8,KM) = QNF(8,KM)+1
      DO 93 I=1,2
93    IF(LL.NE.0 .OR. KNZR.NE.0)
     &            FORML(1,KM,I) = 2.*FORML(2,KM,I) - FORML(3,KM,I)
      DO 95 I=1,NRAD
         R = (I-1)*HDN
      IF(KNZR.EQ.0) THEN
       FORMR(I) = FF(R,RIN,FNN(1,N),NLN) * REL(N)  * WREL(N)**2
      ELSE
       FORMR(I) = FFR4(R*RIN,FORML(1,KM,1),NLN) * MAX(R,EPS)**LL
      ENDIF
95    CONTINUE
      IF(IPC.GE.2 .AND. KNZR.NE.0)
     &     WRITE(KO,955) KM,((I-1)*RINTP,FORMR((I-1)*MR+1),I=1,NLN/2,1)
955   FORMAT('0FORM FACTOR',I3,' IS :' /(4(4X,F8.3,2X,E11.4,4X)) )
      IF(IPC.GE.3 .AND. KNZR.NE.0)
     &WRITE(KO,955) -KM,(I,(I-1)*RINTP,FORML(I,KN,2)*((I-1)*RINTP)**LL
     &  /FORMR((I-1)*MR+1)    ,    I=2,NLN/3)
96    WRITE(8,REC=KM) FORMR
      KN = KN + NMAX
100   CONTINUE
      KN2 = KN - 1
         IF(IPC.GE.1) WRITE(KO,98) WM,WN,COM,SQRT(RMS/MAX(WM,EPS))
     &                                  ,SQRT(RMC/MAX(WM,EPS))
98    FORMAT('0Total square norm of all form factors calculated =',
     &F10.5,',   and of all those stored =',F10.5/
     & '  Overall rms radius of the two nucleons relative to ',
     & 'the c.m. of the ',A8,' =',F9.4,' fm.',
     & ', and of the core =',F9.4,' fm.')
	deallocate (QERN)
      RETURN
      END
      SUBROUTINE QNN(QERN,NLN,NLO,MAXL1,NK,TNT,A,B,P,Q,NL2,FORML,
     &               NNT,MAXNLN,RINTP,RIN,MSP,QNF,IPC,XG,WG,NW,REL)
	use io
	use drier
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8  QERN(NL2,NLN,MAXL1,NK),XG(6),WG(6),REL(NLO)
      INTEGER TNT(4,NK),QNF(19,MSP)
      REAL*8 VR,WF1,WF2,FFR4
      REAL*8 FORML(MAXNLN,MSP,2)
      REAL*8 PLEG(NNT,MAXL1),UK(NNT),WT(NNT),RN1(NNT),RN2(NNT)

C
      AB2 = A*B*2
      PQ2 = P*Q*2
         THMAX = 3.14159
           IC3 = NNT/NW
           NNU = IC3 * NW
            CI = THMAX/IC3
            C1 = 0.5D0 * CI
            C2 = C1
               K = 0
            DO 14  J=1,IC3
            DO 13 NN=1,NW
               K = K + 1
               TH = C1 * XG(NN) + C2
               UK(K) =  COS(TH)
               WT(K) = C1 * WG(NN) * SIN(TH)
               PLEG(K,1) = 1
               if(MAXL1>1) PLEG(K,2) = UK(K)
13             CONTINUE
14           C2 = C2 + CI
               DO 15 L=2,MAXL1-1
               DO 15 K=1,NNU
               WD = UK(K)
15         PLEG(K,L+1) = ((2*L-1)*WD*PLEG(K,L-1+1) -(L-1)*PLEG(K,L-2+1))
     &                       / DBLE(L)
C
      IF(IPC.GE.4) WRITE(KO,*) NLO,(REL(K),K=1,NLO)
      DO 80 I=1,NLN
         RR = (I-1) * RINTP
         DO 29 J=1,NL2
         DO 29 IN=1,NK
         DO 29 L1=1,MAXL1
29       QERN(J,I,L1,IN) = 0.0
C
         DO 50 J=1,NLO
            R = REL(J)
           IF(I+1.GE.NLN- R * MAX(ABS(B),ABS(Q))*RIN)  GO TO 80
            ABR = (A*RR)**2 + (B*R)**2
            PQR = (P*RR)**2 + (Q*R)**2
            RRR  = R*RR
            DO 31 K=1,NNU
               RN22 =     PQR + PQ2*RRR*UK(K)
               RN12 =     ABR + AB2*RRR*UK(K)
               RN2(K) = SQRT(ABS(RN22)) * RIN
31             RN1(K) = SQRT(ABS(RN12)) * RIN
               DO 50 IN=1,NK
               IF(TNT(1,IN).LT.0.OR.TNT(3,IN).LT.0) GO TO 50
                   IF2 = TNT(2,IN)
                   IF1 = TNT(1,IN)
                IF(QNF(8,IF1).EQ.0 .OR. QNF(8,IF2).EQ.0) GO TO 50
                   WR = 0.5
                     IF(TNT(3,IN).GT.0)
     &                  WR = .5 * FFR4(R*RIN,FORML(1,TNT(3,IN),1),NLN)
     &                              * R**QNF(9,TNT(3,IN))
               DO 40 K=1,NNU
                   WF2= FFR4(RN2(K),FORML(1,IF2,1),NLN)
                   WF1= FFR4(RN1(K),FORML(1,IF1,1),NLN)
                      VR=FFR4(RN1(K),FORML(1,IF1,2),NLN) * WF2
     &                  +FFR4(RN2(K),FORML(1,IF2,2),NLN) * WF1
               WD= VR * WT(K) * WR
               WF= (WF1*WF2) * WT(K) * WR
            DO 36 L1=1,MAXL1
36          QERN(J+NLO,I,L1,IN)=QERN(J+NLO,I,L1,IN) + WD * PLEG(K,L1)
               DO 37 L1=1,MAXL1
37          QERN(J,I,L1,IN)    = QERN(J,I,L1,IN)    + WF * PLEG(K,L1)
40         CONTINUE
50      CONTINUE
80    CONTINUE
      RETURN
      END
      SUBROUTINE EXTERN1(ISC,NFL,COEF,LL,L,S,J,J12,JCORE,JCOM,
     &                  PRES,IPC,FNN,NLN,NN,NN2,RINTP,REL,EGS)
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MR=66, MT=56, MMAX=66, MMA=18,  MR1=MR+1,MT1=MT+1)
      REAL*8 S,J,J12,JCORE,JCOM,FNN(NLN,NN2),J3,REL(NN)
      LOGICAL PRES
      REAL*8 THET(0:MT1),R(0:MR1),CO(0:MT1),SI(0:MT1),TR(MMA)
      INTEGER TNQ(MMA,5)
      CHARACTER*10 TA(5)
      REAL*8 XLA(MR+4),XMU(MT+4),RHO(1),TH(1)
      PARAMETER(LWRK=4*MMAX+MMAX+4 , LIWRK=MMAX+MMAX)
      REAL*8 WRK(LWRK)
      INTEGER IWRK(LIWRK)
      real*8,allocatable:: V(:)
       Z = 0.0
C
      REWIND NFL
      READ(NFL,500) TA,H2SM,CX,CY
500   FORMAT(5A10,3E10.3)
      READ(NFL,510) DR,AC,NR,NT,JC1,JC2,JC3
510   FORMAT(2E10.3,5I10)
      READ(NFL,520) EGS
	EGS = abs(EGS)
520   FORMAT(E10.3)
      READ(NFL,530) MA
530   FORMAT(8I10)
         CALL CHECK(NR,MR,28)
         CALL CHECK(NT,MT,28)
         CALL CHECK(MA,MMA,28)
         CALL CHECK(MAX(NR,NT),MMAX,28)
      DO 20 IA=1,MA
20    READ(NFL,530) (TNQ(IA,IQ),IQ=1,5)
540   FORMAT(8E10.3)
C
      CALL GRIPOL(NR,DR,AC,NT,JC1,JC2,JC3, DT,R,THET,CO,SI)
C
      SUM = 0.0
      DO 25 IA=1,MA
      TR(IA) = 0.0
         IF(LL.NE.TNQ(IA,1)/2 .OR.
     &      L .NE.TNQ(IA,3)/2 .OR.
     &      ABS(S*2-TNQ(IA,4)) + ABS(J*2-TNQ(IA,5)) .GT. 0.01) GO TO 25
      J3 = TNQ(IA,2) * 0.5
      TR(IA) = (-1)**NINT(J12-L-J) * SQRT((2.*J12+1.)*(2.*J3+1.))
     &         * RACAH(J,LL+Z,JCOM,JCORE,J12,J3)
      IF(IPC.GE.4) WRITE(KO,30) IA,(TNQ(IA,I)*.5,I=1,5),TR(IA)
30    FORMAT(' From input ch.',I3,' ( nos.',5F5.1,') ',/
     &       '  Transformation coefficient =',F12.5)
C     IF(IPC.GE.5) WRITE(KO,*) LL,L,S,J,J12,JCORE,JCOM,J3
25    SUM = SUM + ABS(TR(IA))
      IF(SUM.LT.1E-20) RETURN
      C = COEF *  (CX * CY) ** (-0.25)
      NTR = NT * NR
      PRES = .TRUE.
	 allocate (V(NTR*2))
C
      DO 110 IA=1,MA
      READ(NFL,540) (V(I),I=1,NTR)
      IF(ABS(TR(IA)).LT.1E-20) GO TO 110
      IFAIL=0
      CALL E01DAF(NR,NT,R(1),THET(1),V,IPX,IPY,XLA,XMU,V,V(NTR+1),IFAIL)
      DO 105 I=2,NLN
         RR = (I-1)*RINTP
         Y = RR / SQRT(CY)
      DO 100 IX=1,NN
         X = REL(IX) / SQRT(CX)
C     IF(IPC.GE.5) WRITE(KO,*) 'I,Y,IX,X =',I,Y,IX,X
         RHO(1) = SQRT(X*X + Y*Y)
         TH(1)  = ATAN(Y/X)
        VAL = 0.0
        if(RHO(1).gt.R(NR)) go to 90
            IFAIL = 001
      CALL E02DFF(1,1,IPX,IPY,RHO,TH,XLA,XMU,V,VAL,
     X           	  WRK,LWRK,IWRK,LIWRK,IFAIL)
 90    FNN(I,IX) = VAL * ( C * TR(IA) / RR) / REL(IX)
c       write(98,*) 'RHO, XLA(4), XLA(IPX-4),VAL=',
c     X       real(RHO(1)), real(XLA(4)), real(XLA(IPX-4)),real(VAL)
100   CONTINUE
105   CONTINUE
110   CONTINUE
	 deallocate (V)
      RETURN
      END
      SUBROUTINE GRIPOL(NR,DR,AC,NT,JC1,JC2,JC3,DT,R,TET,CO,SI)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 R(0:NR)
      REAL*8 TET(0:NT)
      REAL*8 CO(0:NT),SI(0:NT)
      DATA PIS2/1.570796327/
      NULL=0
C****
C**** R   DIMENSION 1+NR+1
C**** TET DIMENSION 1+NT+1
C**** CO  DIMENSION 1+NT+1
C**** SI  DIMENSION 1+NT+1
C****
      DRC=DR
      R(NULL)=0.
      DO 10 IR=1,NR
      R(IR)=R(IR-1)+DRC
      DRC=DRC*AC
   10 CONTINUE
      R(NR+1)=R(NR)+DRC
C****
      K1=8
      TET(NULL)=0.
      DO 15 IT=1,NT
      TET(IT)=TET(IT-1)+K1
      IF(IT.EQ.JC1.OR.IT.EQ.JC2.OR.IT.EQ.JC3)K1=K1/2
   15 CONTINUE
      TET(NT+1)=TET(NT)+K1
      DT=PIS2/TET(NT+1)
      DO 20 IT=1,NT
      TET(IT)=TET(IT)*DT
      CO(IT)=COS(TET(IT))
      SI(IT)=SIN(TET(IT))
   20 CONTINUE
      TET(NT+1)=TET(NT+1)*DT
      CO(NULL)=1.
      SI(NULL)=0.
      CO(NT+1)=0.
      SI(NT+1)=1.
      RETURN
      END
      SUBROUTINE EXTERN2(ISC,NFL,COEF,LL,L,S,J,J12,JCORE,JCOM,
     &           KIND,PRES,IPC,FNN,NLN,NN,NN2,RINTP,REL,WREL,MM,EGS)
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 S,J,J12,JCORE,JCOM,FNN(NLN,NN2),REL(NN),JC,JTOT,
     &        	WREL(NN)
      LOGICAL PRES,SKIP
      REAL*8,allocatable:: RV(:),WK(:),XLA(:),XMU(:),WRK(:),V(:)
      INTEGER,allocatable:: IWRK(:)
       Z = 0.0
       FPI = 4d0*atan(1d0)*4d0
C
        CALL CHECK(NLN,MM,27)
        CALL CHECK(NN,MM,27)
	if(IPC>4) write(6,*) 'Read 2NT wf from file ',NFL
      REWIND NFL
      READ(NFL,'(i4,f10.5,i4,f10.5)',err=2) NRXY,RSTEP,I,EGS
2     if(IPC>4) write(6,'(i4,f10.5,i4,f10.5)') NRXY,RSTEP,I,EGS
	EGS = abs(EGS)
      NRXY=NRXY+1
      MRXY=max(NRXY,NLN)
      MRXY2=MRXY**2
        LWRK=4*MM+MRXY+4 ; LIWRK=MM+MRXY
	allocate (RV(MRXY),WK((MRXY+6)**2),XLA(MRXY+4),XMU(MRXY+4),
     X		  WRK(4*MM+MRXY+4),IWRK(MM+MRXY))
         CALL CHECK(NRXY,MRXY,26)
         CALL CHECK(NLN,MRXY,26)
	 allocate (V(2*MRXY**2))
       WN = 0.0
       VN = 0.0
      DO 200 IVERT=1,3
	SKIP = IVERT.ge.2 .and. IVERT.ne.abs(ISC)
      READ(NFL,*) 
      READ(NFL,*) 
       SUM = 0.0
      IF(IPC.GE.5) WRITE(KO,*) 'Reqd: IVERT,LL,L,S,J,J12,JCORE,JCOM',
     X       			IVERT,LL,L,S,J,J12,JCORE,JCOM
      DO 110 IA=1,100
      READ(NFL,915) LL1,LL2,LLL,ISS,JNN,IC,JC,JTOT
       if(LL1.lt.0) go to 120
915       format(3x,6i3,2f4.1)
       IF(IPC.GE.5) WRITE(KO,*) 'Have:',LL1,LL2,LLL,ISS,JNN,IC,JC,JTOT

C
      TR = 0.0
         IF(LL.NE.LL2 .OR. L .NE.LL1 .OR.
     &      ABS(S-ISS) + ABS(J12-JNN) .GT. 0.01) GO TO 40
         IF(KIND.eq.9.and.abs(JCORE-JC)+abs(JCOM-JTOT).gt.0.01) go to 40
      TR = (-1)**NINT(J12-LLL-S) * SQRT((2.*J+1.)*(2.*LLL+1.))
!!!  &         * RACAH(L+Z,LL+Z,JNN+Z,S,LLL+Z,J)   ! WRONG
     &         * RACAH(LL+Z,L+Z,JNN+Z,S,LLL+Z,J)   ! Order corrected 18 May 12 for L/=LL cases!
     &         * (-1)**LLL    ! factor added, for efaddy coordinates, 3 Nov 08
      IF(IPC.GE.4)WRITE(KO,30) IA,LL1,LL2,LLL,ISS,JNN,IC,JC,JTOT,TR,COEF
30    FORMAT(' From input ch.',I3,' ( nos.',6i3,2f4.1,') ',
     X       '    Transformation coefficient =',2F12.5)
      SUM = SUM + ABS(TR)
C
40      RV(1) = 0.0
        DO 50 IX=2,NRXY
       RV(IX) = (IX-1)*RSTEP
       read(NFL,*) IIX
50       READ(NFL,*) (V(MRXY2+IY+(IX-1)*NRXY),IY=2,NRXY)
	if(SKIP) then
	  read(NFL,*) 
	  go to 110
	  endif
       DO 52 IX=1,NRXY
       IY=1
         V(MRXY2+IY+(IX-1)*NRXY) = 0.0
52       V(MRXY2+IX+(IY-1)*NRXY) = 0.0
       XN = 0.0
        DO 53 IX=1,NRXY
        DO 53 IY=1,NRXY
       if(IVERT.eq.1) XN = XN + V(MRXY2+IY+(IX-1)*NRXY)**2 * RSTEP**2
53     if(IVERT.gt.1) XN = XN + V(MRXY2+IY+(IX-1)*NRXY) * RSTEP**2 *FPI

        read(NFL,*) PNORM
      IF(ABS(TR).LT.1E-20) GO TO 110
      IF(IPC.GE.5) WRITE(KO,*) 'Have:',LL1,LL2,LLL,ISS,JNN,IC,JC,JTOT
      IF(IPC.GE.5) WRITE(KO,*) 'Input partial integral:',real(XN)
      IFAIL=0
      CALL E01DAF(NRXY,NRXY,RV,RV,V(MRXY2+1),IPX,IPY,XLA,XMU,V,WK,IFAIL)

      NLNI = 0
      do 55 I=1,NLN
       RV(I) = (I-1)*RINTP
55    if(RV(I).lt.(NRXY-1)*RSTEP) NLNI = I
      NNI = 0
      do 56 IX=1,NN
56    if(REL(IX).lt.(NRXY-1)*RSTEP) NNI = IX

      call flush(6)
      IFAIL=0
      CALL E02DFF(NNI,NLNI,IPX,IPY,REL,RV,XLA,XMU,V,V(MRXY2+1),
     X           	  WRK,LWRK,IWRK,LIWRK,IFAIL)
      DO 105 IY=2,NLNI
      DO 100 IX=1,NNI
       VAL = V(MRXY2+NLNI*(IX-1)+IY) / (RV(IY)*REL(IX))
       if(IVERT.eq.1) FNN(IY,IX) = VAL * COEF * TR 
       if(IVERT.gt.1) FNN(IY,IX+NN) = VAL * COEF * TR
       if(IPC.ge.5) then
       if(IVERT.eq.1) WN = WN + FNN(IY,IX)**2 * WREL(IX)**2
     X       		* (RV(IY)*REL(IX))**2
       if(IVERT.gt.1) VN = VN + FNN(IY,IX)*FNN(IY,IX+NN) * WREL(IX)**2
     X       		* (RV(IY)*REL(IX))**2
       endif
100   CONTINUE
105   CONTINUE
110   CONTINUE
120        read(NFL,*) TNORM
c      IF(IPC.GE.5) WRITE(KO,*) 'Input summed integrals:',real(TNORM)
200   CONTINUE
	 deallocate (V,RV,WK,XLA,XMU,WRK,IWRK)
      IF(SUM.LT.1E-20) RETURN
      PRES = .TRUE.
      IF(IPC.lt.5) return
       WN = WN*RINTP
       VN = VN*RINTP
        WRITE(KO,*) 'Stored wf norm :',real(WN)
        WRITE(KO,*) 'Vertex potential me:',real(VN)
      RETURN
      END

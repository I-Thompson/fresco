!***********************************************************************
!     
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

*****FRXX10*************************************************************
      SUBROUTINE KERNEM(FFR,KCOEF,NLL,NONO,IREM,ICV,NM,LDMIN,LDMAX,
     &               C1,C2,IC1,IC2,REPEAT,  LVAL,ICTO,ICFROM,REV,PART,
     &               NK,NG,FPT,GPT,CHNO,CP, MCG,MAXMV,MKNL,NL0,QNF,
     &               KNL,NUMLT,LTMIN,LTRANS,IC7,MXLN1,MXLNP1,MMX1,
     X               MXMV1,MXMVP1)
	use parameters
	use io
	use kcom
	use drier
	use trace
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KCOEF(MAXQRN,NUMLT),MCG(3,MAXMV,NG,MKNL)
      INTEGER LVAL(MAXCH),D,C,PART(MAXCH),GPT(2,NG),CP,C1,C2,
     &    QNF(19,MSP),CHNO(MFNL,6),FPT(6,NK),NM(MAXQRN,MFNL)
      LOGICAL REV,PRES,IFAIL3,THERE,FFR,REPEAT,C1FR,NONO,
     X        LTRANS(MAXQRN)
      DATA ONE/1D0/
      IFAIL3(I,J,K) = I.GT.J+K .OR. I.LT.ABS(J-K)
C
!        write(0,*) ' Called KERNEM'

      Z = 0.0
      EPS = 1E-14
      PI = 4.0*ATAN(ONE)
      PISQ8 = 8.0*PI**2
      REPEAT = .TRUE.
      C1FR = ICFROM.EQ.IC1
      IF(C1FR) THEN
         C = C1
         D = C2
      ELSE
         C = C2
         D = C1
      ENDIF
      IF(ICTO.NE.PART(D).OR.ICFROM.NE.PART(C)) STOP 7
C       SO TO 'D' FROM 'C' CHANNEL
C  AND WITH NAME(1,IC1) LIKE D & NAME(1,IC2) LIKE P IN (D,P) REACTION
      LD = LVAL(D)
      LC = LVAL(C)
      LDMIN = MIN(LD,LDMIN)
      LDMAX = MAX(LD,LDMAX)
      THERE = .FALSE.
      DO 8 IN=1,NG
      DO 8 ILT=1,NUMLT
 8     IF(ABS(KCOEF(IN,ILT)).GT.EPS.AND..NOT.LTRANS(IN)) THERE = .TRUE.
      IF(.NOT.THERE) GO TO 100

      NL = NL + 1
      DRY = DRY .OR. NL.GT.MFNL
      IF(NL.GT.MFNL)  WRITE(KO,9) NL,MFNL
9     FORMAT(//' ****** NOT ENOUGH ROOM FOR',I4,' NON-LOCAL FORM',
     &' FACTORS  IN MFNL ARRAY OF',I4,' *****')
      KNL = NL-NL0
      CHNO(NL,1) = D
      IF(.NOT.REV) CHNO(NL,1) = -D
      CHNO(NL,2) = C
      CHNO(NL,3) = 2
      IF(NONO) CHNO(NL,3) = IC7
      IF(.NOT.FFR) CHNO(NL,3) = CHNO(NL,3) + 10
      CHNO(NL,4) = NLL
      CHNO(NL,6) = 1
      PRES = .FALSE.
C
      DO 40 IN=1,NG
       NM(IN,KNL) = 0
       IF(LTRANS(IN)) GO TO 40
C
C   HERE, TO START :    GPT(1  IS NO. OF B.S. IN CHANNEL 'C1' ("DEUT")
C                   AND GPT(2  IS NO. OF B.S. IN CHANNEL 'C2' ("PROT")
         IF(C1FR) THEN
           KNP = GPT(1,IN)
           KN = GPT(2,IN)
         ELSE
           KNP = GPT(2,IN)
           KN  = GPT(1,IN)
         ENDIF
C
C   NOW, MORE PROPERLY, KNP IS NO. OF B.S. IN CHANNEL 'C' (FROM)
C                   AND KN  IS NO. OF B.S. IN CHANNEL 'D' (TO)
         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
         MXLN1 = MAX(MXLN1,LN+1)
         MXLNP1= MAX(MXLNP1,LNP+1)
         THERE = .FALSE.
         DO 16 ILT=1,NUMLT
 16       IF(ABS(KCOEF(IN,ILT)).GT.EPS) THERE = .TRUE.
         IF(.NOT.THERE) GO TO 40
      IF(LISTCC.GT.1) WRITE(KO,18) NL,D,C,IN,KNP,KN,NUMLT,LTMIN
 18   FORMAT('0NL. M-interaction #',I3,' to',I3,' from',I3,' :',5I3)
C
      IM = 0
      DO 30 MVP=-LNP,LNP
      DO 30 MV =-LN,LN
        MM = MVP - MV
        IF(ABS(MM).GT.LD) GO TO 30
        IF(MM.LT.0) GO TO 30
          R0 = 2.0
          IF(MM.EQ.0 .AND. MV.EQ.0) R0 = 1.0
          IF(MM.EQ.0 .AND. MV.LT.0) GO TO 30
        IM = IM + 1
        MCG(1,IM,IN,KNL) = MV
        MCG(2,IM,IN,KNL) = MVP
        MMX1 = MAX(ABS(MM)+1,MMX1)
        MXMV1 = MAX(ABS(MV)+1,MXMV1)
        MXMVP1 = MAX(ABS(MVP)+1,MXMVP1)
        MCG(3,IM,IN,KNL) = 0.0
        R1 = PISQ8 * YLMC(LN,MV) * YLMC(LNP,MVP)
     X             * YLMC(LC,0) *  YLMC(LD,MM)
       if(C1FR) R1 = R1 * (-1)**iabs(MM)
         DO 25 ILT=1,NUMLT
            LTOTAL = LTMIN + ILT - 1
            IF(LTOTAL.LT.0.OR.ABS(MM+MV).GT.LTOTAL) GO TO 25
            IF(IFAIL3(LD,LN,LTOTAL).OR.IFAIL3(LC,LNP,LTOTAL)) GO TO 25
          R2 = KCOEF(IN,ILT)
          R3 = WIG3J(LD+Z,LN+Z, LTOTAL+Z,MM+Z,MV+Z, Z-MM-MV)
          R4 = WIG3J(LC+Z,LNP+Z,LTOTAL+Z,Z   ,MVP+Z,Z-MVP)
          R = R0 * R1 * R2 * R3 * R4
          IF(LISTCC.GT.2) WRITE(KO,23) D,C,IN,MVP,MV,LTOTAL,
     X          R0,R1,R2,R3,R4,R
23        FORMAT('0KERNEM:',6I3,5F8.4,F9.6)
          PRES = PRES .OR. ABS(R) .GT. EPS
         MCG(3,IM,IN,KNL) = MCG(3,IM,IN,KNL) + R
25       CONTINUE
          IF(LISTCC.GT.1) WRITE(KO,23) D,C,IN,MVP,MV,IM,MCG(3,IM,IN,KNL)
30    CONTINUE
      NM(IN,KNL) = IM
40    CONTINUE
      IF(.NOT.PRES) NL = NL - 1
100   RETURN
      END
      SUBROUTINE QERNEM(NLN,NLM,NLO,NK,NG,XA,XB,XP,XQ,FPT,CUTOFF,LTRANS,
     &             MCG,MAXMV,NL0,KNL,NM,MMX1,C1FR,
     &             IC7,ICV,IREM,NNT,RINTO,EPC,HNL,NLT,LRANG1,LDMIN,CP,
     &             VFOLD,VCORE,THM,RINS,RINC,NLL,QNF,NLOC,
     &             FORML,VSP,MINT,RIN,NLC,CENTRE,NONO,HF,HT,
     &             FFREAL,JACOBN,LISTCC,CHNO,RERR,NCH,LVAL,NLPLOT,
     X             MXMV1,MXMVP1,MXLN1,MXLNP1,FORMC,VSC,cxwf)
	use parameters
	use io
	use drier
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MCG(3,MAXMV,NG,KNL),ALT(NNT,MMX1,LRANG1),
     X       ALN(NNT,MXMV1,MXLN1),ALNP(NNT,MXMVP1,MXLNP1)
      REAL*8 JACOBN,COFFIN(4),XG(6),WG(6),C2,CI,TH,THM(MAXNLN)
      INTEGER FPT(6,NK),QNF(19,MSP),CUTOFF,NM(MAXQRN,KNL),CHNO(MFNL,6),
     X        CP,LVAL(NCH),KLJ(NLM),KLNEW(NLM)
      COMPLEX*16 VCOR(MAXNNU),VSUBI,DV,FFC,SUMT,FNC(NLL,NLO,KNL),
     X          VFOLD(MAXNLN),VCORE(MAXNLN),SUM(NNT,KNL),WC(MAXNNU),
     X          FORMC(MAXNLC,MSP),VSC(MAXNLC,MSP),WF1C,WF2C,VC,FFC4
      LOGICAL FFREAL,C1FR,NONO,LTRANS(NK),REP,cxwf
      REAL*8 FNL(NLL,NLO,KNL),NLOC(NLM),VR,WF2,WF1,FFR,FFR4,RCOR2,S,T,
     X       VSP(MAXNLN,MSP),FORML(MAXNLN,MSP),
     X       RUM(NNT,KNL),WD(MAXNNU),CG(MAXNNU)
      REAL*8 UK(MAXNNU),GK(MAXNNU),COST(MAXNNU),COSF(MAXNNU),
     &       WT(MAXNNU),R1(MAXNNU),R2(MAXNNU),RCOR(MAXNNU)
      DATA NWW,XG(4),XG(5),XG(6),WG(4),WG(5),WG(6)/3,
     1   .2386191861D0,.6612093865D0,.9324695142D0,
     2   .4679139346D0,.3607615730D0,.1713244924D0/, PI /3.14159D0/
C
!	write(0,*) ' Called QERNEM'
      CALL CCTIME(IPTIME)
      EP = EPC * 0.01
      IF(EP.lt.1e-5) EP = 1E-5
      HIRS = 1.0/RINS
      HIRC = 1.0/RINC
      RFRTM = 1.0/((NLL-1)*RINTO)**2
      NW=2*NWW
      DO 1 N=1,NWW
      NN=NW-N+1
      XG(N)=-XG(NN)
  1   WG(N)=WG(NN)
      MMX = MMX1 - 1
      MXLN = MXLN1 - 1
      MXLNP= MXLNP1 -1
      MXMV = MXMV1 - 1
      MXMVP= MXMVP1 -1
C
!	   if(LISTCC>1) write(KO,601)  'VFOLD:',VFOLD
!	   if(LISTCC>1) write(KO,601)  'VCORE:',VCORE
!601	   format(1x,a6,' potential:',/(1x,1p,6g14.5))
C
           IC3 = NNT/NW
           NNU = IC3 * NW
            CI = 1D0/IC3
            C1 = 0.5D0 * CI
C
C   In QERNEM, we use a pure initial/final set of parameters,
C      as we are not concerned with any projectile/target asymmetries.
C      (These exist in the coupling order, but the transfer kernel
C
C       <LT,LN; LTOTAL | Vbind + Vcore - Vfold | LF,LNP; LTOTAL>
C
C       being calculated here IS symmetric)
C
C    FROM channel LF, radius RF=ri, M-projection = 0,
C         internal b.s radius R1, l,M = LNP,MVP, state KNF
C
C    TO   channel LT, radius RT=rf, M-projection = MM,
C         internal b.s radius R2, l,M = LN,MV,   state KNT
C
C   NOTE THAT A,B,P,Q ADHERE TO THIS CONVENTION:
C       R1 = P * RF + Q * RT = P * ri + Q * rf     (ri & rf in Buttle's
C       R2 = A * RF + B * RT = A * ri + B * rf       DAISY)
C   AS, IN QERNEL, A,B, P,Q ARE DEFINED FOR:
C       RDINT(proj) = XP * RF + XQ * RT
C       RN(targ)    = XA * RF + XB * RT
C
      IF(C1FR) THEN
C                  from: R1=RDINT, to: R2=RN
        A = XA
        B = XB
        P = XP
        Q = XQ
       ELSE
C                  from: R1=RN, to: R2=RDINT
        A = XP
        B = XQ
        P = XA
        Q = XB
       ENDIF
      PA = P - A
      QB = Q - B
      HFHT = HF/HT
      PQ = P*HFHT + Q
      AB = A*HFHT + B
      PAQB2 = PA*QB*2.
      PAQB  = PA*HFHT+QB
      AB2 = A*B*2
      PQ2 = P*Q*2
      IF(LISTCC.GE.3) WRITE(KO,20) ' A,B,P,Q,JACOBN',A,B,P,Q,JACOBN
 20    FORMAT(1X,'QERNEM: ', A20 / (1X,10F8.4))
C 21    FORMAT(1X, 9F8.4)
         THMAX = PI
         DO 5 IL=1,KNL
         DO 5 I=1,NLL
         DO 5 JJ=1,NLO
          IF(FFREAL)      FNL(I,JJ,IL) = 0.0
          IF(.NOT.FFREAL) FNC(I,JJ,IL) = 0.0
5       CONTINUE
C------------------------------------------START R-TO LOOP
        REP = .FALSE.
        XTI = 0.0
	KLIM= 0
      DO 100 I=2,NLL
         RT = (I-1) * RINTO
         IF(RT.LT.(CUTOFF-1)*HT) GO TO 100
         DO 8 J=1,NLM
8        KLNEW(J) = NNU
         IF(I.GT.3.AND.REP) THEN
           IF(KLIM.EQ.NNU) TH = MIN(THMAX*1.333,PI)
           IF(KLIM.LT.NNU) TH = (ACOS(UK(KLIM)+1.) + THMAX) * 0.5
C        WRITE(KO,*) ' AT I =',I,', KLJ(*)=',(KLJ(J),J=1,NLM)
           DO 10 J=1,NLM
             T = ACOS(UK(KLJ(J)) + 1.) / TH
             S = COS(T * THMAX) - 1.0
             DO 9 K=1,NNU
9            IF(UK(K).LT.S) GO TO 10
             K = NNU
10           KLNEW(J) = K
             THMAX = TH
          ENDIF
       if(LISTCC.GE.5) then
	 WRITE(KO,105) I,KLIM,THMAX,XTI
105      format(' AT I =',i4,', KLIM,THMAX,XTI =',i4,2g20.8)
         WRITE(KO,11) (KLNEW(J),J=1,NLM,4)
11       FORMAT(' KLNEW(*)=',/,(1x,20i4))
	endif
         KLIM = NNU/2
         DO 12 J=1,NLM
 12      KLJ(J) = NNU/4
         THM(I) = THMAX
         REP = .TRUE.
         XTL = XTI
         XTI = 0.0
            C2 = C1
               K = 0
            DO 14  J=1,IC3
            DO 13 NN=1,NW
               K = K + 1
               X = C1 * XG(NN) + C2
               TH = (3D0*X*X + 1D0)*X*THMAX*0.25D0
               R1(K) =  COS(TH)
               UK(K) = R1(K)- 1D0
               GK(K) = C1 * WG(NN)
               WT(K) =  SIN(TH) * (9D0 *X*X+1D0 )*THMAX*0.25D0
13             CONTINUE
14           C2 = C2 + CI
       LDMAX = LDMIN + LRANG1 -1
       CALL PLMV(R1,NNU,LDMAX,MMX,NNT,MMX1,LRANG1,ALT,KLTT,R2)
C      IF(LISTCC.GE.3)
C
      VSUBI = 0.0
      DO 25 K=1,NNU
25    VCOR(K) = 0.0
      IF(IREM.NE.0.AND.IC7.EQ.0) THEN
         IF(.NOT.FFREAL) VSUBI = FFC(RT*HIRS,VFOLD,NLN)
         IF(FFREAL)      VSUBI = FFR(RT*HIRS,VFOLD,NLN)
       ENDIF
C------------------------------------------START R-FROM LOOP
         G = 4./3.
         G = 1.
         DO 85 J=1,NLM
            DNL = (J - NLM/2 - 1) * HNL  +  CENTRE
            RT0 = RT*HFHT
            RF  = RT0    + DNL
            IF(RF.LT.(CUTOFF-1)*HF .OR. RF.LE.0.01) GO TO 85
            PQR = P * DNL + PQ * RT
            ABR = A * DNL + AB * RT
            PQRS = PQR*PQR
            ABRS = ABR*ABR
            PAQBR = (PA * DNL + PAQB*RT)**2
            RFRT = RF*RT
C           NNJ = KLNEW(J)   (this shortcut gives inaccurate results!)
            NNJ = NNU
               DO 30 IL=1,KNL
               DO 30 K=1,NNJ
                IF(FFREAL)      RUM(K,IL) = 0.0
                IF(.NOT.FFREAL) SUM(K,IL) = 0.0
  30            CONTINUE
            DO 31 K=1,NNJ
               R1SQ  =  PQRS + PQ2*RFRT*UK(K)
               R2SQ  =  ABRS + AB2*RFRT*UK(K)
               R1(K) =  SQRT(ABS(R1SQ))
               R2(K) =  SQRT(ABS(R2SQ))
C              COST(K) = (A * RF + B * RT * (UK(K)+1.0))/R2(K)
C 31           COSF(K) = (P * RF + Q * RT * (UK(K)+1.0))/R1(K)
               COST(K) = (ABR + B * RT * UK(K))/R2(K)
  31           COSF(K) = (PQR + Q * RT * UK(K))/R1(K)
C
C Find Plm coefficients for 0<l<MXLN, 0 <M< MXLN(P),  COS(T/F),K=1,NNJ
C      results in ALN(K,M+1,KL),   RCOR temporary storage.
         CALL PLMV(COST,NNJ,MXLN, MXMV, NNT,MXMV1, MXLN1, ALN, KL, RCOR)
         CALL PLMV(COSF,NNJ,MXLNP,MXMVP,NNT,MXMVP1,MXLNP1,ALNP,KLP,RCOR)
C
           IF(IREM.NE.0) THEN
            IF(IC7.EQ.1.AND..NOT.FFREAL) VSUBI= FFC(RF*HIRS,VFOLD,NLN)
            IF(IC7.EQ.1.AND.     FFREAL) VSUBI= FFR(RF*HIRS,VFOLD,NLN)
                  DO 315 K=1,NNJ
                  RCOR2 = PAQBR + PAQB2*RFRT * UK(K)
 315              RCOR(K) = SQRT(ABS(RCOR2)) * HIRC
               DO 316 K=1,NNJ
 316           VCOR(K)=FFC(RCOR(K),VCORE,NLN) - VSUBI
            ENDIF
C
               DO 50 INF=1,NK
                   IN = FPT(3,INF)
                   IF(LTRANS(IN)) GO TO 50
C  Sum all the INF with the coefficient MCG(3,IM,IN,IL) (i.e. using IN)
                     IF(C1FR) THEN
                       KNF = FPT(1,INF)
                       KNT = FPT(2,INF)
                     ELSE
                       KNF = FPT(2,INF)
                       KNT = FPT(1,INF)
                     ENDIF
                   LNP= QNF(9,KNF)
                   LN = QNF(9,KNT)
                DO 40 K=1,NNJ
                if(.not.cxwf) then
                 WF1= FFR4(R1(K)*RIN,FORML(1,KNF),MINT)
                 WF2= FFR4(R2(K)*RIN,FORML(1,KNT),MINT)
                   IF(NONO) THEN
                      VR = WF1 * WF2
                   ELSE
                 IF(ICV.EQ.0) VR=FFR4(R1(K)*RIN,VSP(1,KNF),MINT) * WF2
                 IF(ICV.EQ.1) VR=FFR4(R2(K)*RIN,VSP(1,KNT),MINT) * WF1
                   ENDIF
                else  ! cxwf
                 WF1C= FFC4(R1(K)*RIN,FORMC(1,KNF),MINT)
                 WF2C= FFC4(R2(K)*RIN,FORMC(1,KNT),MINT)
                   IF(NONO) THEN
                      VC = WF1C * WF2C
                   ELSE
                 IF(ICV.EQ.0) VC=FFC4(R1(K)*RIN,VSC(1,KNF),MINT) * WF2C
                 IF(ICV.EQ.1) VC=FFC4(R2(K)*RIN,VSC(1,KNT),MINT) * WF1C
                   ENDIF
                 endif ! cxwf
         IF(.NOT.FFREAL) THEN  ! necessarily, if cxwf=T
C                            Complex form factors
                if(.not.cxwf) then
                  DV = VR
                  IF(IREM.NE.0) DV =  DV + VCOR(K) * WF2 * WF1
                else
                  DV = VC
                  IF(IREM.NE.0) DV =  DV + VCOR(K) * WF2C * WF1C
                endif
               WC(K) = DV * WT(K)

                   W = ABS(WC(K))
           ELSE
C                        Real form factors
                IF(IREM.NE.0) VR = VR + DBLE(VCOR(K)) *WF2*WF1
               WD(K) = VR * WT(K)
                   W = ABS(WD(K))
           ENDIF
           IF(W.GT.EP*XTL) KLJ(J) = MAX(KLJ(J),K)
             XTI = MAX(XTI,W)
             KLIM = MAX(KLIM,KLJ(J))
 40        CONTINUE
C
           DO 49 IL=1,KNL
              IF(CHNO(IL+NL0,6).EQ.0) GO TO 49
              LT = LVAL(IABS(CHNO(IL+NL0,1)))
              KLT =MOD(LT,LRANG1)+1
               DO 48 IM=1,NM(IN,IL)
                  MV  = NINT(MCG(1,IM,IN,IL))
                  MVP = NINT(MCG(2,IM,IN,IL))
                  MM = MVP - MV
               DO 43 K=1,NNJ
 43            CG(K) =  GK(K) * MCG(3,IM,IN,IL)
     X                * ALT(K,IABS(MM)+1,KLT)
     X                * ALN(K,IABS(MV)+1,LN+1)
     X                * ALNP(K,IABS(MVP)+1,LNP+1)
               IF(FFREAL) THEN
                 DO 45 K=1,NNJ
 45              RUM(K,IL) = RUM(K,IL) + WD(K) * CG(K)
               ELSE
                 DO 46 K=1,NNJ
 46              SUM(K,IL) = SUM(K,IL) + WC(K) * CG(K)
               ENDIF
 48            CONTINUE
 49            CONTINUE
C
 50    CONTINUE
       C2 = RFRT * RFRTM
       CALL SPLINT( DNL/(HNL*NLT) + NLC    ,NLO,IJ,NP,COFFIN)
       G = 2. - G
       S  = (G / NLT) * JACOBN  * RFRT
       DO 80 IL=1,KNL
              IF(CHNO(IL+NL0,6).EQ.0) GO TO 80
               IF(FFREAL) THEN
                  RUMT = 0.0
                  DO 55 K=1,NNJ
 55               RUMT = RUMT + RUM(K,IL)
                  NLOC(J)=NLOC(J) + RUMT**2 * C2
               ELSE
                  SUMT = 0.0
                  DO 58 K=1,NNJ
 58               SUMT = SUMT + SUM(K,IL)
                  NLOC(J)=NLOC(J) + ABS(SUMT)**2 * C2
               ENDIF
        DO 70 M=1,NP
            T  =  COFFIN(M) * S
            JJ = IJ + M - 1
         IF(FFREAL) THEN
            FNL(I,JJ,IL) = FNL(I,JJ,IL) + RUMT * T
         ELSE
            FNC(I,JJ,IL) = FNC(I,JJ,IL) + SUMT * T
         ENDIF
70      CONTINUE
80    CONTINUE
85    CONTINUE
100   CONTINUE
      IF(DRY) RETURN
C
      DO 150 IL=1,KNL
      IF(CHNO(IL+NL0,6).EQ.0) WRITE(KO,*) 'Q ERROR: ',IL,' NOT REQD!'
C!UC  IF(MACH.EQ.7) CHNO(IL+NL0,5) = GETPOS(12)
      IF(FFREAL) THEN
             WRITE(12) ((FNL(I,J,IL),I=1,NLL),J=1,NLO)
      ELSE
             WRITE(12) ((FNC(I,J,IL),I=1,NLL),J=1,NLO)
      ENDIF
C
      IF(NLPLOT.LE.0) GO TO 150
      WRITE(KO,110) CHNO(IL+NL0,1),CHNO(IL+NL0,2),CP
 110  FORMAT(' NL m-interaction v-c(',I4,') c(',I4,') of CP=',I3,' is'/)
      IF(FFREAL) CALL DISPLY(FNL(1,1,IL),NLL,NLO,NLL,SCALE)
      IF(.NOT.FFREAL) THEN
         CALL DISPLR(FNC(1,1,IL),NLL,NLO,NLL,SCALR)
           SCALI = SCALR * RERR * 1E3
         CALL DISPLI(FNC(1,1,IL),NLL,NLO,NLL,SCALI)
         SCALE = MAX(SCALR,SCALI)
         ENDIF
      IF(SCALE.lt.1e-10) GO TO 150
      NLPLOT = NLPLOT - 1
C       IF(FFREAL) THEN
C         DO 95 I=2,NLL
C         FNLD(I,1) = 0
C         DO 95 J=1,NLO
C95       FNLD(I,1) = FNLD(I,1) + FNL(I,J,IL) * HNL *MLT
C         WRITE(KO,98) (FNLD(I,1),I=2,NLL)
C98       FORMAT(1X,18F7.2)
C       ENDIF
150   CONTINUE
C
      CALL CCTIME(IPTIME)
       PT = IPTIME * 0.01
      IF(LISTCC.GE.4) WRITE(KO,992) PT,(THM(I)*180/PI,I=2,NLL)
992   FORMAT(/' Theta - maxima (after',F6.2,' secs) are'/(1X,20F6.1))
      RETURN
      END
      SUBROUTINE NLSTAT(ICP,NLOC,NLM,EP,HNL,CENTRE)
	use parameters
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 NLOC(NLM),XM,VR
C
      XM = 0.0
      DO 210 J=1,NLM
210   XM = MAX(NLOC(J),XM)
      VR = 100.*EP * 0.1
      IF(ABS(XM).LT.1.E-30) GO TO 219
      JMIN = 10000
      JMAX = -10000
      DO 215 J=1,NLM
            DNL = (J - NLM/2-1) * HNL + CENTRE
      NLOC(J) = SQRT( NLOC(J) / XM) * 100.
      IF(NLOC(J).LE.VR) GO TO 215
         JMAX = MAX(JMAX,J)
         JMIN = MIN(JMIN,J)
215   CONTINUE
      IF(NLOC(1).GT.VR)
     X JMIN = 1 - MAX(INT(LOG(NLOC(1)/VR)/LOG(NLOC(2)/NLOC(1))),0)
      IF(NLOC(NLM).GT.VR)
     &JMAX =NLM+MAX(INT(LOG(NLOC(NLM)/VR)/LOG(NLOC(NLM-1)/NLOC(NLM))),0)
      W  = (JMAX - JMIN) * HNL
      CI = ((JMAX + JMIN) * 0.5 - NLM/2-1) * HNL + CENTRE
219   WRITE(KO,220) ICP,W,CI,JMIN,JMAX,(NLOC(J),J=1,NLM)
220   FORMAT('0For Coupling ',I3,':',
     x     ' Recommended non-local width is greater than',F8.2,' fm., ',
     &     'and Centration ',F6.2,'  (from Jmin,max =',2I6,')',/,
     &     '0Relative non-local usages (rms) are'/(1X,15F8.3))
      RETURN
      END
      SUBROUTINE PLMV(X,NX,N,M,NXD,NAD,NC,PL,KL,S)
C
C Find P(l,m) coefficients for 0 < l <N, 0 < m < M,  for all X(k),k=1,NX
C  Results in PL(k,m+1,KL), for KL= l mod NC +1.   S = temporary storage
C
      use io
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PL(NXD,NAD,NC),L,X(NX),S(NX)
      IC(I) = MOD(I-1,NC) + 1
      N1 = N+1
      M1 = M+1
      IF(N1.GT.NC.AND.NC.LT.3.OR.NX.GT.NXD.OR.M1.GT.NAD) THEN
         WRITE(KO,*) 'NX,N,M,NXD,NAD,NC =',NX,N,M,NXD,NAD,NC
         STOP 'PLMV'
         ENDIF
      DO 10 I=1,NC
      DO 10 J=1,M1
      DO 10 K=1,NX
10    PL(K,J,I) = 0.
      DO 11 K=1,NX
      PL(K,1,1) = 1.
11    S(K) = SQRT(ABS(1.-X(K)*X(K)))
      IF(N1.GE.2) THEN
        DO 12 K=1,NX
12      PL(K,1,2) = X(K)
      ENDIF
      IF(N1.GE.2.AND.M1.GE.2) THEN
        DO 14 K=1,NX
14      PL(K,2,2) = S(K)
      ENDIF
      DO 20 I=3,N1
      L = I-1
      RL = 1D0/L
CDIR$ IVDEP
      DO 15 K=1,NX
15    PL(K,1,IC(I))=((2.*L-1.)*X(K)*PL(K,1,IC(I-1))
     X                     - (L-1.)*PL(K,1,IC(I-2)))*RL
      JM = MIN(I,M1)
      DO 20 J=2,JM
CDIR$ IVDEP
      DO 20 K=1,NX
      PL(K,J,IC(I)) = (2.*L-1.)*S(K)*PL(K,J-1,IC(I-1)) + PL(K,J,IC(I-2))
20    CONTINUE
      KL = IC(N1)
      RETURN
      END
****ASYMPTOPIA************************************************************
	SUBROUTINE ASYMPTOPIA(NCH,ipmax,CI,STREN,LVAL,MASS,PART,
     X    ryev,LAMBDA,NF,CLIST,NFLIST,NCLIST,CF)
	use parameters
      	implicit none
	integer nch,LVAL(nch),IF,ipmax,NF,LAMBDA(NF),
     x   NCLIST(MAXCH,MAXCH),
     X   NFLIST(MAXCH,MAXCH,MCLIST),NC,C,C2,i,IC,PART(MAXCH,3)
	REAL*8 cfp(ipmax),STREN(NF),T,CF(MAXCH,MAXCH,ipmax),
     X	 	MASS(4,MXP+1),ryev
	complex*16 CI,CLIST(MAXCH,MAXCH,MCLIST),S
C
        do 344 c=1,NCH
	 IC = PART(c,1)
        do 344 c2=1,NCH
          S = CI**(-LVAL(C2)+LVAL(C))
	  cfp(1:ipmax) = 0d0
	  do NC=1,NCLIST(c,c2)
	    IF = NFLIST(c,c2,NC)
	    t = S*CLIST(c,c2,NC)
	    i = LAMBDA(IF)+1
  	    if(i>1) cfp(i) = cfp(i)+t*STREN(IF)
	  enddo
	  if(c==c2) cfp(1) = cfp(1) + MASS(3,IC) * MASS(4,IC) * COULCN
344	CF(c,c2,1:ipmax) = cfp(1:ipmax)/ryev

        RETURN
	END SUBROUTINE ASYMPTOPIA
      SUBROUTINE WRTMAT (A,N,M,L,NCOL,IOUT)
      IMPLICIT INTEGER (A-Z)
      REAL*8 A
      DIMENSION A(L,*)
*
*     WRTMAT : PRINTS AN N*M MATRIX STORED IN AN L*L ARRAY
*
      NTIM = M / NCOL
      NF = 0
      IF ( NTIM .GT. 0 ) THEN
         DO 20 I = 1, NTIM
            NI = NF + 1
            NF = NF + NCOL
            DO 10 J = 1, N
               WRITE (IOUT,1000) (A(J,K),K = NI,NF)
   10       CONTINUE
            WRITE (IOUT,1010)
   20    CONTINUE
      ENDIF
      NI = NF + 1
      IF ( NI .LE. M ) THEN
         DO 40 J = 1, N
            WRITE (IOUT,1000) (A(J,K),K = NI,M)
   40    CONTINUE
      ENDIF
      RETURN
*
c1000 FORMAT (7G16.7)
c1000 FORMAT (7F16.8)
c1000 FORMAT (8F13.8)
 1000 FORMAT (10G13.4)
 1010 FORMAT (/)
      END

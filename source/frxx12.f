      SUBROUTINE A02AAF(XXR,XXI,YR,YI)
      DOUBLE PRECISION  XXI, XXR, YI, YR
      DOUBLE PRECISION  H, HALF, ONE, XI, XR, ZERO
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
      INTRINSIC         ABS, SQRT
      DATA              ZERO/0.0D0/, HALF/0.5D0/, ONE/1.0D0/
      XR = ABS(XXR)
      XI = XXI
      IF (XR.GT.ONE) H = SQRT(XR*HALF+A02ABF(XR*HALF,XI*HALF))
      IF (XR.LE.ONE) H = SQRT(XR+A02ABF(XR,XI))*SQRT(HALF)
      IF (XI.NE.ZERO) XI = XI/(H+H)
      IF (XXR.LT.ZERO) GO TO 20
      YR = H
      YI = XI
      RETURN
   20 IF (XI.LT.ZERO) GO TO 40
      YR = XI
      YI = H
      RETURN
   40 YR = -XI
      YI = -H
      RETURN
      END
      DOUBLE PRECISION FUNCTION A02ABF(XXR,XXI)
      DOUBLE PRECISION                 XXI, XXR
      DOUBLE PRECISION                 H, ONE, XI, XR, ZERO
      INTRINSIC                        ABS, SQRT
      DATA                             ZERO/0.0D0/, ONE/1.0D0/
      XR = ABS(XXR)
      XI = ABS(XXI)
      IF (XI.LE.XR) GO TO 20
      H = XR
      XR = XI
      XI = H
   20 IF (XI.NE.ZERO) GO TO 40
      A02ABF = XR
      RETURN
   40 H = XR*SQRT(ONE+(XI/XR)**2)
      A02ABF = H
      RETURN
      END
      SUBROUTINE A02ACF(XXR,XXI,YYR,YYI,ZR,ZI)
      DOUBLE PRECISION  XXI, XXR, YYI, YYR, ZI, ZR
      DOUBLE PRECISION  A, H, ONE
      INTRINSIC         ABS
      DATA              ONE/1.0D0/
      IF (ABS(YYR).LE.ABS(YYI)) GO TO 20
      H = YYI/YYR
      A = ONE/(H*YYI+YYR)
      ZR = (XXR+H*XXI)*A
      ZI = (XXI-H*XXR)*A
      RETURN
   20 H = YYR/YYI
      A = ONE/(H*YYR+YYI)
      ZR = (H*XXR+XXI)*A
      ZI = (H*XXI-XXR)*A
      RETURN
      END
      SUBROUTINE C02AFF(A,N,SCALE,Z,WORK,IFAIL)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C02AFF')
      INTEGER           IFAIL, N
      LOGICAL           SCALE
      DOUBLE PRECISION  A(2,0:N), WORK(4*(N+1)), Z(2,N)
      DOUBLE PRECISION  BIG
      INTEGER           I, IER, NDEG, NREC
      LOGICAL           SC
      CHARACTER*80      REC(2)
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
      EXTERNAL          C02AFZ
      INTRINSIC         SQRT
      IER = IFAIL
      SC = SCALE
      NREC = 0
      NDEG = N
      IF (N.LT.1) THEN
         IER = 1
         WRITE (REC,FMT=99999) N
         NREC = 2
      ELSE IF ((A(1,0).EQ.ZERO) .AND. (A(2,0).EQ.ZERO)) THEN
         IER = 1
         WRITE (REC,FMT=99998)
         NREC = 1
      ELSE
         BIG = 1.0D0/(SQRT(2.0D0)*X02AMF())
         DO 10 I = 1, N
            Z(1,I) = -BIG
            Z(2,I) = -BIG
   10    CONTINUE
         CALL C02AFZ(A,NDEG,SC,Z,WORK(1),WORK(2*N+3),IER)
         IF (IER.EQ.0) THEN
            IFAIL = 0
            GO TO 20
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
   20 RETURN
99999 FORMAT (' ** On entry, N.lt.1:',/'    N = ',I16)
99998 FORMAT (' ** On entry, A(1,0).eq.0 and A(2,0).eq.0')
      END
      SUBROUTINE C02AFW(ARI,AII,BRI,BII,CR,CI,FAIL)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  AII, ARI, BII, BRI, CI, CR
      LOGICAL           FAIL
      DOUBLE PRECISION  AI, AR, BI, BR, BIG, DIV, FLMAX, FLMIN, NUMI,
     *                  NUMR, TEMP
      LOGICAL           FIRST
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
      INTRINSIC         ABS, MAX, SIGN
      SAVE              BIG, FIRST, FLMAX
      DATA              FIRST/.TRUE./
      AR = ARI
      AI = AII
      BR = BRI
      BI = BII
      IF (AR.EQ.ZERO .AND. AI.EQ.ZERO) THEN
         CR = ZERO
         CI = ZERO
         IF (BR.EQ.ZERO .AND. BI.EQ.ZERO) THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
         IF (FIRST) THEN
            FIRST = .FALSE.
            FLMIN = X02AMF()
            FLMAX = 1/FLMIN
            BIG = FLMAX/2
         END IF
         TEMP = MAX(ABS(AR),ABS(AI),ABS(BR),ABS(BI))
         IF (TEMP.GE.BIG) THEN
            AR = AR/2
            AI = AI/2
            BR = BR/2
            BI = BI/2
         END IF
         IF (BR.EQ.ZERO .AND. BI.EQ.ZERO) THEN
            CR = SIGN(FLMAX,AR)
            CI = SIGN(FLMAX,AI)
            FAIL = .TRUE.
         ELSE
            IF (ABS(BR).GE.ABS(BI)) THEN
               TEMP = BI/BR
               DIV = BR + TEMP*BI
               NUMR = AR + TEMP*AI
               NUMI = AI - TEMP*AR
            ELSE
               TEMP = BR/BI
               DIV = BI + TEMP*BR
               NUMR = AI + TEMP*AR
               NUMI = TEMP*AI - AR
            END IF
            IF (ABS(DIV).GE.ONE) THEN
               CR = NUMR/DIV
               CI = NUMI/DIV
               FAIL = .FALSE.
            ELSE
               TEMP = ABS(DIV)*FLMAX
               IF ((ABS(NUMR).LE.TEMP) .AND. (ABS(NUMI).LE.TEMP)) THEN
                  CR = NUMR/DIV
                  CI = NUMI/DIV
                  FAIL = .FALSE.
               ELSE
                  IF (DIV.GE.ZERO) THEN
                     CR = SIGN(FLMAX,NUMR)
                     CI = SIGN(FLMAX,NUMI)
                  ELSE
                     CR = SIGN(FLMAX,-NUMR)
                     CI = SIGN(FLMAX,-NUMI)
                  END IF
                  FAIL = .TRUE.
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
      SUBROUTINE C02AFX(AR,AI,BR,BI,CR,CI,ZSM,ZLG)
      DOUBLE PRECISION  HALF, ONE, ZERO, TWO
      PARAMETER         (HALF=0.5D0,ONE=1.0D0,ZERO=0.0D0,TWO=2.0D0)
      DOUBLE PRECISION  AI, AR, BI, BR, CI, CR
      DOUBLE PRECISION  ZLG(2), ZSM(2)
      DOUBLE PRECISION  DEPS, FINITY, SQRTFY, SQRTTY, TINY
      INTEGER           EMAX, EMIN, EXPDEP, LRGEXP
      LOGICAL           OVFLOW, UNFLOW
      DOUBLE PRECISION  SC, XC1, XC2
      INTEGER           EXPBSQ, SCLEXP
      DOUBLE PRECISION  A(2), B(2), C(2), CT(2), D(2)
      DOUBLE PRECISION  C02AGR
      INTEGER           C02AGX
      EXTERNAL          C02AGR, C02AGX
      EXTERNAL          A02AAF, A02ACF
      INTRINSIC         ABS, MAX, MIN
      COMMON            /AC02AF/OVFLOW, UNFLOW
      COMMON            /BC02AF/DEPS, FINITY, SQRTFY, SQRTTY, TINY,
     *                  EMAX, EMIN, EXPDEP, LRGEXP
      DOUBLE PRECISION  APPABS
      APPABS(XC1,XC2) = MAX(ABS(XC1),ABS(XC2))
      A(1) = AR
      A(2) = AI
      B(1) = -BR
      B(2) = -BI
      C(1) = CR
      C(2) = CI
      IF (APPABS(A(1),A(2)).NE.ZERO) THEN
         IF (APPABS(C(1),C(2)).NE.ZERO) THEN
            SCLEXP = (C02AGX(APPABS(A(1),A(2)))+C02AGX(APPABS(C(1),C(2))
     *               ))/2
            IF (APPABS(B(1),B(2)).NE.ZERO) THEN
               EXPBSQ = 2*(C02AGX(APPABS(B(1),B(2)))-SCLEXP)
            ELSE
               EXPBSQ = -2*EXPDEP
            END IF
            IF (EXPBSQ.LE.EXPDEP) THEN
               SCLEXP = MIN(SCLEXP+1,EMAX)
               SCLEXP = MAX(SCLEXP,EMIN)
               SC = C02AGR(ONE,SCLEXP)
               IF (EXPBSQ.LT.-EXPDEP) THEN
                  B(1) = ZERO
                  B(2) = ZERO
               ELSE
                  B(1) = (B(1)/SC)*HALF
                  B(2) = (B(2)/SC)*HALF
               END IF
               A(1) = A(1)/SC
               A(2) = A(2)/SC
               C(1) = C(1)/SC
               C(2) = C(2)/SC
               CT(1) = B(1)*B(1) - B(2)*B(2) - A(1)*C(1) + A(2)*C(2)
               CT(2) = TWO*B(2)*B(1) - A(2)*C(1) - A(1)*C(2)
               CALL A02AAF(CT(1),CT(2),D(1),D(2))
               IF (D(1)*B(1)+D(2)*B(2).LE.ZERO) THEN
                  D(1) = -D(1)
                  D(2) = -D(2)
               END IF
               B(1) = B(1) + D(1)
               B(2) = B(2) + D(2)
            END IF
            CALL A02ACF(B(1),B(2),A(1),A(2),ZLG(1),ZLG(2))
            CALL A02ACF(C(1),C(2),B(1),B(2),ZSM(1),ZSM(2))
         ELSE
            ZSM(1) = ZERO
            ZSM(2) = ZERO
            CALL A02ACF(B(1),B(2),A(1),A(2),ZLG(1),ZLG(2))
         END IF
      ELSE
         OVFLOW = .TRUE.
         ZLG(1) = FINITY
         ZLG(2) = ZERO
         IF (APPABS(B(1),B(2)).EQ.ZERO .AND. APPABS(C(1),C(2)).NE.ZERO)
     *       THEN
            ZSM(1) = -ZLG(1)
            ZSM(2) = -ZLG(2)
         ELSE
            IF (APPABS(B(1),B(2)).EQ.ZERO) THEN
               ZSM(1) = ZLG(1)
               ZSM(2) = ZLG(2)
            ELSE
               CALL A02ACF(C(1),C(2),B(1),B(2),ZSM(1),ZSM(2))
            END IF
         END IF
      END IF
      RETURN
      END
      SUBROUTINE C02AFY(DX,DY,NDEG,A,P,PPRIME,PDPRIM,ERROR,DEFLAT)
      DOUBLE PRECISION  SXTEEN, THREE, TWO, TO3RDS, ZERO
      PARAMETER         (SXTEEN=16.0D0,THREE=3.0D0,TWO=2.0D0,
     *                  TO3RDS=TWO/THREE,ZERO=0.0D0)
      DOUBLE PRECISION  DX, DY, ERROR
      INTEGER           NDEG
      DOUBLE PRECISION  A(2,0:NDEG), DEFLAT(2,0:NDEG), P(2), PDPRIM(2),
     *                  PPRIME(2)
      DOUBLE PRECISION  DEPS, FINITY, SQRTFY, SQRTTY, TINY
      INTEGER           EMAX, EMIN, EXPDEP, LRGEXP
      LOGICAL           OVFLOW, UNFLOW
      DOUBLE PRECISION  ABSX, DT, DVI, DVR, WI, WR
      INTEGER           I
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
      INTRINSIC         DBLE
      COMMON            /AC02AF/OVFLOW, UNFLOW
      COMMON            /BC02AF/DEPS, FINITY, SQRTFY, SQRTTY, TINY,
     *                  EMAX, EMIN, EXPDEP, LRGEXP
      IF (NDEG.GE.1) THEN
         ABSX = A02ABF(DX,DY)
         WR = ZERO
         WI = ZERO
         DVR = A(1,0)
         DVI = A(2,0)
         DEFLAT(1,0) = DVR
         DEFLAT(2,0) = DVI
         DEFLAT(1,1) = A(1,1) + (DX*DEFLAT(1,0)-DY*DEFLAT(2,0))
         DEFLAT(2,1) = A(2,1) + (DX*DEFLAT(2,0)+DY*DEFLAT(1,0))
         DO 20 I = 2, NDEG
            DT = DVR + (DX*WR-DY*WI)
            WI = DVI + (DX*WI+DY*WR)
            WR = DT
            DT = DEFLAT(1,I-1) + (DX*DVR-DY*DVI)
            DVI = DEFLAT(2,I-1) + (DX*DVI+DY*DVR)
            DVR = DT
            DEFLAT(1,I) = A(1,I) + (DX*DEFLAT(1,I-1)-DY*DEFLAT(2,I-1))
            DEFLAT(2,I) = A(2,I) + (DX*DEFLAT(2,I-1)+DY*DEFLAT(1,I-1))
   20    CONTINUE
         P(1) = DEFLAT(1,NDEG)
         P(2) = DEFLAT(2,NDEG)
         IF ( .NOT. OVFLOW) THEN
            ERROR = TO3RDS*A02ABF(A(1,0),A(2,0))
            DO 40 I = 1, NDEG - 1
               ERROR = A02ABF(DEFLAT(1,I),DEFLAT(2,I)) + ABSX*ERROR
   40       CONTINUE
            ERROR = SXTEEN*DEPS*(A02ABF(DEFLAT(1,NDEG),DEFLAT(2,NDEG))
     *              +THREE*ABSX*ERROR)
            IF (OVFLOW) THEN
               ERROR = DBLE(NDEG)*ERROR
               OVFLOW = .FALSE.
            END IF
            PPRIME(1) = DVR
            PPRIME(2) = DVI
            PDPRIM(1) = WR
            PDPRIM(2) = WR
         END IF
      ELSE IF (NDEG.EQ.0) THEN
         P(1) = A(1,0)
         P(2) = A(2,0)
         ERROR = ZERO
         PPRIME(1) = ZERO
         PPRIME(2) = ZERO
         PDPRIM(1) = ZERO
         PDPRIM(2) = ZERO
      END IF
      RETURN
      END
      SUBROUTINE C02AFZ(A,NDEG,SCALE,Z,DU,DEFLAT,IER)
      DOUBLE PRECISION  GAMA, THETA
      PARAMETER         (GAMA=1.0D0,THETA=2.0D0)
      DOUBLE PRECISION  HALF, ONE, SMALL, BIGONE, SMLONE, RCONST,
     *                  ONEPQT, ZERO, TWO
      PARAMETER         (HALF=0.5D0,ONE=1.0D0,SMALL=1.0D-3,
     *                  BIGONE=1.0001D0,SMLONE=0.99999D0,RCONST=1.445D0,
     *                  ONEPQT=1.25D0,ZERO=0.0D0,TWO=2.0D0)
      INTEGER           IER, NDEG
      LOGICAL           SCALE
      DOUBLE PRECISION  A(2,0:NDEG), DEFLAT(2,0:NDEG), DU(2,0:NDEG),
     *                  Z(2,NDEG)
      DOUBLE PRECISION  DEPS, FINITY, SQRTFY, SQRTTY, TINY
      INTEGER           EMAX, EMIN, EXPDEP, LRGEXP
      LOGICAL           OVFLOW, UNFLOW
      DOUBLE PRECISION  ABDIR, ABDIRO, ABSCL, DX, DZ0I, DZ0R, DZNI,
     *                  DZNR, E, F0, FEJER, FN, G, LOWERB, MXCOEF, R,
     *                  RATIO, RTN, S, T, UPPERB, X2N, X2N1, XC1, XC2,
     *                  XN, XN1, XN2, XN2N
      INTEGER           I, IERS, IHALF, ISPIR, ITER, K, MXCFEX, N, NERR,
     *                  SCBYEX
      LOGICAL           CAUCHY, CONTIN, OVF, SAVO, SAVU, SPIRAL, STARTD,
     *                  UNF
      DOUBLE PRECISION  C(2), CDIR(2), CDIRO(2), CF(2), CF1(2), CF2(2),
     *                  CL(2), CR(2), CSPIR(2), CTEMP(2)
      CHARACTER*80      REC(2)
      DOUBLE PRECISION  A02ABF, C02AGY, X02AJF, X02AKF, X02ALF
      INTEGER           C02AGX, X02BJF, X02BKF, X02BLF
      LOGICAL           C02AGS
      EXTERNAL          A02ABF, C02AGY, X02AJF, X02AKF, X02ALF, C02AGX,
     *                  X02BJF, X02BKF, X02BLF, C02AGS
      EXTERNAL          A02ACF, C02AFW, C02AFX, C02AFY, X04AAF, X04BAF
      INTRINSIC         ABS, EXP, LOG, MAX, MIN, DBLE, SQRT
      COMMON            /AC02AF/OVFLOW, UNFLOW
      COMMON            /BC02AF/DEPS, FINITY, SQRTFY, SQRTTY, TINY,
     *                  EMAX, EMIN, EXPDEP, LRGEXP
      DOUBLE PRECISION  APXABS
      APXABS(XC1,XC2) = ABS(XC1) + ABS(XC2)
      TINY = X02AKF()
      SQRTTY = SQRT(TINY)
      FINITY = X02ALF()
      SQRTFY = SQRT(FINITY)
      EXPDEP = X02BJF() + 1
      EMIN = X02BKF() - 1
      EMAX = X02BLF() - 1
      LRGEXP = EMAX + 1 - EXPDEP
      DEPS = X02AJF()
      IERS = IER
      IER = 0
      ITER = 0
      IHALF = 0
      ISPIR = 0
      N = NDEG
      SAVO = OVFLOW
      SAVU = UNFLOW
      OVF = .FALSE.
      UNF = .FALSE.
      MXCOEF = ZERO
      DO 20 I = 0, N
         DU(1,I) = A(1,I)
         DU(2,I) = A(2,I)
         MXCOEF = MAX(MXCOEF,APXABS(A(1,I),A(2,I)))
   20 CONTINUE
      IF (MXCOEF.EQ.ZERO) THEN
         DO 40 I = 1, N
            Z(1,I) = FINITY
            Z(2,I) = ZERO
   40    CONTINUE
         N = 0
         OVF = .TRUE.
      ELSE
         MXCFEX = C02AGX(MXCOEF)
         IF (MXCFEX.GT.LRGEXP) THEN
            SCBYEX = 0
            SCALE = .FALSE.
         ELSE
            SCBYEX = LRGEXP - MXCFEX
         END IF
      END IF
      CAUCHY = .FALSE.
   60 IF (N.GT.2) THEN
         IF (SCALE) THEN
            IF (SCBYEX.NE.0) THEN
               DO 80 I = 0, N
                  DU(1,I) = C02AGY(DU(1,I),SCBYEX)
                  DU(2,I) = C02AGY(DU(2,I),SCBYEX)
   80          CONTINUE
               SCBYEX = 0
            END IF
         END IF
         UNF = UNFLOW .OR. UNF
         DO 100 I = 0, N - 1
            IF (C02AGS(DU(1,I)) .AND. C02AGS(DU(2,I))) THEN
               Z(1,N-I) = FINITY
               Z(2,N-I) = ZERO
            ELSE
               GO TO 120
            END IF
  100    CONTINUE
  120    IF (I.NE.0) THEN
            DO 140 K = I, N
               DU(1,K-I) = DU(1,K)
               DU(2,K-I) = DU(2,K)
  140       CONTINUE
            N = N - I
            IF (SCBYEX.EQ.-EXPDEP) THEN
               IER = 3
               IF (IERS.NE.1) THEN
                  CALL X04AAF(0,NERR)
                  WRITE (REC,FMT=99999)
                  CALL X04BAF(NERR,REC(1))
                  CALL X04BAF(NERR,REC(2))
               END IF
               GO TO 300
            END IF
            OVF = .TRUE.
            GO TO 60
         END IF
         DO 160 I = N, 1, -1
            IF (C02AGS(DU(1,I)) .AND. C02AGS(DU(2,I))) THEN
               Z(1,I) = ZERO
               Z(2,I) = ZERO
            ELSE
               GO TO 180
            END IF
  160    CONTINUE
  180    IF (I.NE.N) THEN
            N = I
            GO TO 60
         END IF
         OVFLOW = .FALSE.
         UNFLOW = .FALSE.
         IF ( .NOT. CAUCHY) THEN
            XN = DBLE(N)
            XN1 = DBLE(N-1)
            XN2 = DBLE(N-2)
            X2N = TWO/XN
            X2N1 = X2N/XN1
            XN2N = XN2/XN
            RTN = SQRT(XN)
            G = EXP((LOG(A02ABF(DU(1,N),DU(2,N)))-LOG(A02ABF(DU(1,0),
     *          DU(2,0))))/XN+SMALL)
            OVFLOW = .FALSE.
            CALL C02AFW(DU(1,N-1),DU(2,N-1),DU(1,N),DU(2,N),CR(1),CR(2),
     *                  OVFLOW)
            IF (OVFLOW) THEN
               Z(1,N) = ZERO
               Z(2,N) = ZERO
               N = N - 1
               GO TO 60
            END IF
            CTEMP(1) = X2N1*DU(1,N-2)
            CTEMP(2) = X2N1*DU(2,N-2)
            CF2(1) = X2N*DU(1,N-1)
            CF2(2) = X2N*DU(2,N-1)
            CALL C02AFX(CTEMP(1),CTEMP(2),CF2(1),CF2(2),DU(1,N),DU(2,N),
     *                  C,CF1)
            CR(1) = XN2N*CR(1)
            CR(2) = XN2N*CR(2)
            CTEMP(1) = (C(1)*CR(1)-C(2)*CR(2)) + XN1
            CTEMP(2) = C(2)*CR(1) + C(1)*CR(2)
            CALL A02ACF(C(1),C(2),CTEMP(1),CTEMP(2),CDIRO(1),CDIRO(2))
            ABDIRO = APXABS(CDIRO(1),CDIRO(2))
            G = MIN(G,BIGONE*MIN(APXABS(C(1),C(2)),RTN*ABDIRO))
            R = G
            S = BIGONE*G
            UNFLOW = .FALSE.
            DO 200 I = 0, N
               DEFLAT(1,I) = A02ABF(DU(1,I),DU(2,I))
  200       CONTINUE
  220       IF (R.LT.S) THEN
               T = DEFLAT(1,0)
               S = ZERO
               OVFLOW = .FALSE.
               DO 240 I = 1, N - 1
                  S = R*S + T
                  T = R*T + DEFLAT(1,I)
  240          CONTINUE
               S = R*S + T
               T = (R*T-DEFLAT(1,N))/S
               S = R
               R = R - T
               GO TO 220
            END IF
            IF (OVFLOW) THEN
               SCBYEX = -EXPDEP
               GO TO 60
            END IF
            CAUCHY = .TRUE.
            UPPERB = MIN(RCONST*XN*R,G)
            LOWERB = SMLONE*S
            UNF = UNFLOW .OR. UNF
         END IF
         FEJER = UPPERB
         G = UPPERB
         CDIR(1) = CDIRO(1)
         CDIR(2) = CDIRO(2)
         ABDIR = ABDIRO
         RATIO = ABDIR/G
         DZNR = ZERO
         DZNI = ZERO
         FN = A02ABF(DU(1,N),DU(2,N))
         F0 = FN
         SPIRAL = .FALSE.
         STARTD = .FALSE.
         CONTIN = .TRUE.
  260    IF (CONTIN) THEN
            ITER = ITER + 1
            IF (RATIO.GT.THETA) THEN
               IF (STARTD) THEN
                  IHALF = IHALF + 1
                  ABSCL = HALF*ABSCL
                  CL(1) = HALF*CL(1)
                  CL(2) = HALF*CL(2)
                  DX = ABS(DZNR) + ABS(DZNI)
                  IF (DX+ABSCL.NE.DX) THEN
                     DZNR = DZ0R + CL(1)
                     DZNI = DZ0I + CL(2)
                  ELSE
                     IF (FN.GE.E*XN**2) THEN
                        IER = 2
                        IF (IERS.NE.1) THEN
                           CALL X04AAF(0,NERR)
                           WRITE (REC,FMT=99997)
                           CALL X04BAF(NERR,REC(1))
                           CALL X04BAF(NERR,REC(2))
                        END IF
                        GO TO 300
                     END IF
                     CONTIN = .FALSE.
                     GO TO 260
                  END IF
               ELSE
                  ISPIR = ISPIR + 1
                  IF (SPIRAL) THEN
                     C(1) = CSPIR(1)*DZNR - CSPIR(2)*DZNI
                     C(2) = CSPIR(2)*DZNR + CSPIR(1)*DZNI
                  ELSE
                     SPIRAL = .TRUE.
                     CSPIR(1) = -ONEPQT/XN
                     CSPIR(2) = ONE
                     ABSCL = LOWERB/XN**2
                     CTEMP(1) = CDIR(1)/ABDIR
                     CTEMP(2) = CDIR(2)/ABDIR
                     C(1) = CTEMP(1)*LOWERB
                     C(2) = CTEMP(2)*LOWERB
                  END IF
                  DZNR = C(1)
                  DZNI = C(2)
               END IF
            ELSE
               STARTD = .TRUE.
               IF (RATIO.GT.GAMA .AND. (STARTD .OR. SPIRAL .OR.
     *             LOWERB.LE.GAMA*G)) THEN
                  RATIO = GAMA/RATIO
                  CDIR(1) = CDIR(1)*RATIO
                  CDIR(2) = CDIR(2)*RATIO
                  ABDIR = ABDIR*RATIO
               END IF
               G = FEJER
               CL(1) = CDIR(1)
               CL(2) = CDIR(2)
               ABSCL = ABDIR
               F0 = FN
               DZ0R = DZNR
               DZ0I = DZNI
               DZNR = DZ0R + CL(1)
               DZNI = DZ0I + CL(2)
            END IF
            OVFLOW = .FALSE.
            UNFLOW = .FALSE.
            CALL C02AFY(DZNR,DZNI,N,DU,CF,CF1,CF2,E,DEFLAT)
            FN = A02ABF(CF(1),CF(2))
            IF (OVFLOW) THEN
               SCBYEX = -EXPDEP
               GO TO 60
            END IF
            IF (FN.LE.E .OR. UNFLOW) THEN
               IF (UNFLOW) THEN
                  IER = 3
                  IF (IERS.NE.1) THEN
                     CALL X04AAF(0,NERR)
                     WRITE (REC,FMT=99998)
                     CALL X04BAF(NERR,REC(1))
                     CALL X04BAF(NERR,REC(2))
                  END IF
                  UNF = .TRUE.
                  GO TO 300
               END IF
               CONTIN = .FALSE.
               GO TO 260
            END IF
            UNF = UNFLOW .OR. UNF
            IF (FN.GE.F0 .AND. STARTD) THEN
               RATIO = BIGONE*THETA
               GO TO 260
            END IF
            OVFLOW = .FALSE.
            CALL C02AFW(CF1(1),CF1(2),CF(1),CF(2),CR(1),CR(2),OVFLOW)
            IF (OVFLOW) THEN
               UNF = .TRUE.
               CONTIN = .FALSE.
               GO TO 260
            END IF
            CF2(1) = X2N1*CF2(1)
            CF2(2) = X2N1*CF2(2)
            CTEMP(1) = X2N*CF1(1)
            CTEMP(2) = X2N*CF1(2)
            CALL C02AFX(CF2(1),CF2(2),CTEMP(1),CTEMP(2),CF(1),CF(2),C,
     *                  CF1)
            FEJER = APXABS(C(1),C(2))
            CR(1) = XN2N*CR(1)
            CR(2) = XN2N*CR(2)
            CTEMP(1) = (C(1)*CR(1)-C(2)*CR(2)) + XN1
            CTEMP(2) = C(2)*CR(1) + C(1)*CR(2)
            CALL A02ACF(C(1),C(2),CTEMP(1),CTEMP(2),CDIR(1),CDIR(2))
            ABDIR = APXABS(CDIR(1),CDIR(2))
            RATIO = ABDIR/G
            FEJER = MIN(RTN*ABDIR,FEJER)
            DX = ABS(DZNR) + ABS(DZNI)
            IF (DX+ABDIR.EQ.DX) THEN
               CONTIN = .FALSE.
               GO TO 260
            END IF
            GO TO 260
         END IF
         DO 280 I = 1, N - 1
            DU(1,I) = DEFLAT(1,I)
            DU(2,I) = DEFLAT(2,I)
  280    CONTINUE
         Z(1,N) = DZNR
         Z(2,N) = DZNI
         N = N - 1
         CAUCHY = .FALSE.
         GO TO 60
      END IF
      OVFLOW = .FALSE.
      UNFLOW = .FALSE.
      IF (N.EQ.2) THEN
         CALL C02AFX(DU(1,0),DU(2,0),DU(1,1),DU(2,1),DU(1,2),DU(2,2),
     *               CTEMP,C)
         Z(1,1) = C(1)
         Z(2,1) = C(2)
         Z(1,2) = CTEMP(1)
         Z(2,2) = CTEMP(2)
      ELSE IF (N.EQ.1) THEN
         CALL A02ACF(-DU(1,1),-DU(2,1),DU(1,0),DU(2,0),Z(1,1),Z(2,1))
      ELSE
         OVF = OVF .OR. OVFLOW
         UNF = UNF .OR. UNFLOW
         OVFLOW = SAVO
         UNFLOW = SAVU
         IF (OVF) R = FINITY*FINITY
         IF (UNF) R = TINY*TINY
      END IF
  300 RETURN
99999 FORMAT (' ** C02AFF cannot evaluate p(z) near some of its zeros ',
     *       'without overflow.',/' ** If this message occurs please c',
     *       'ontact NAG.')
99998 FORMAT (' ** C02AFF cannot evaluate p(z) near some of its zeros ',
     *       'without underflow.',/' ** If this message occurs please ',
     *       'contact NAG.')
99997 FORMAT (' ** The method has failed. This error is very unlikely ',
     *       'to occur.',/' ** Please contact NAG.')
      END
      DOUBLE PRECISION FUNCTION C02AGR(X,EXP)
      DOUBLE PRECISION                 X
      INTEGER                          EXP
      DOUBLE PRECISION                 C02AGY
      INTEGER                          C02AGX
      EXTERNAL                         C02AGY, C02AGX
      C02AGR = C02AGY(X,EXP-C02AGX(X))
      RETURN
      END
      LOGICAL FUNCTION C02AGS(X)
      DOUBLE PRECISION        ZERO
      PARAMETER               (ZERO=0.0D0)
      DOUBLE PRECISION        X
      C02AGS = (X+ZERO) .EQ. ZERO
      RETURN
      END
      INTEGER FUNCTION C02AGX(DX)
      DOUBLE PRECISION        ONE, ZERO
      PARAMETER               (ONE=1.0D0,ZERO=0.0D0)
      DOUBLE PRECISION        DX
      DOUBLE PRECISION        DEPS, DPNEWL, DPNEWU, FACT, TEMP
      INTEGER                 DBASE, MNEXP, MXEXP, NEWL, NEWU
      DOUBLE PRECISION        A, ABSX
      INTEGER                 E
      LOGICAL                 FIRST
      DOUBLE PRECISION        C02AGY, X02AJF, X02AKF, X02ALF
      INTEGER                 X02BHF, X02BKF, X02BLF
      EXTERNAL                C02AGY, X02AJF, X02AKF, X02ALF, X02BHF,
     *                        X02BKF, X02BLF
      INTRINSIC               ABS, DBLE, LOG
      COMMON                  /CC02AG/DPNEWL, DPNEWU, DEPS, TEMP, FACT,
     *                        DBASE, MNEXP, MXEXP, NEWL, NEWU
      SAVE                    /CC02AG/, FIRST
      DATA                    FIRST/.TRUE./
      IF (FIRST) THEN
         FIRST = .FALSE.
         DBASE = X02BHF()
         DEPS = X02AJF()
         MNEXP = X02BKF()
         MXEXP = X02BLF()
         DPNEWL = X02AKF()
         DPNEWU = X02ALF()
         NEWU = MXEXP - 1
         NEWL = MNEXP - 1
         TEMP = DBLE(DBASE)*(ONE-DEPS)
         FACT = DPNEWU/TEMP
      END IF
      IF (DX.NE.ZERO) THEN
         ABSX = ABS(DX)
         E = LOG(ABSX)/LOG(DBLE(DBASE))
         IF (E.GE.MXEXP) THEN
            E = MXEXP - 1
         ELSE IF (E.LT.MNEXP) THEN
            E = MNEXP
         END IF
         A = ABSX/C02AGY(ONE,E)
   20    IF (A.GE.ONE) THEN
            E = E + 1
            A = A/DBASE
            GO TO 20
         ELSE IF (A.LT.ONE/DBASE) THEN
            E = E - 1
            A = A*DBASE
            GO TO 20
         END IF
      ELSE
         E = 0
      END IF
      C02AGX = E
      RETURN
      END
      DOUBLE PRECISION FUNCTION C02AGY(DX,EXP)
      DOUBLE PRECISION                 ONE
      PARAMETER                        (ONE=1.0D0)
      DOUBLE PRECISION                 DX
      INTEGER                          EXP
      DOUBLE PRECISION                 DEPS, DPNEWL, DPNEWU, FACT, TEMP
      INTEGER                          DBASE, MNEXP, MXEXP, NEWL, NEWU
      DOUBLE PRECISION                 DPE, DSC, POWER
      INTEGER                          E
      INTRINSIC                        MOD
      COMMON                           /CC02AG/DPNEWL, DPNEWU, DEPS,
     *                                 TEMP, FACT, DBASE, MNEXP, MXEXP,
     *                                 NEWL, NEWU
      SAVE                             /CC02AG/
      E = EXP
      DSC = DX
   20 IF (E.GT.NEWU) THEN
         DSC = DSC*FACT
         E = E - NEWU
         GO TO 20
      END IF
   40 IF (E.LT.NEWL) THEN
         DSC = DSC*DPNEWL
         E = E - NEWL
         GO TO 40
      END IF
      IF (E.EQ.0) THEN
         DPE = ONE
      ELSE
         IF (E.LT.0) THEN
            E = -E
            POWER = ONE/DBASE
         ELSE
            POWER = DBASE
         END IF
         DPE = ONE
   60    IF (MOD(E,2).EQ.1) DPE = DPE*POWER
         E = E/2
         IF (E.GT.0) THEN
            POWER = POWER*POWER
            GO TO 60
         END IF
      END IF
      C02AGY = DSC*DPE
      RETURN
      END
      SUBROUTINE E01DAF(MX,MY,X,Y,F,PX,PY,LAMDA,MU,C,WRK,IFAIL)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01DAF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      INTEGER           IFAIL, MX, MY, PX, PY
      DOUBLE PRECISION  C(MX*MY), F(MX*MY), LAMDA(MX+4), MU(MY+4),
     *                  WRK((MX+6)*(MY+6)), X(MX), Y(MY)
      INTEGER           IDUM1, IDUM2, IE, IER2, IERR, IWGHT, IWRK,
     *                  IXKNOT, IXROW, IXU, IYKNOT, IYROW, IYU, JERROR,
     *                  LWRK, MDIST, NE, NREC, NXKNTS, NXU, NYKNTS, NYU
      CHARACTER*80      REC(2)
      INTEGER           P01ABF
      EXTERNAL          P01ABF
      EXTERNAL          E01DAX, E01DAY, E01DAZ
      IERR = 0
      NREC = 0
      IF (MX.LT.4) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) MX
      ELSE IF (MY.LT.4) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) MY
      ELSE
         NXKNTS = MX - 4
         NYKNTS = MY - 4
         NE = MX*MY
         NXU = 2*(MX+NXKNTS+1)
         NYU = 2*(MY+NYKNTS+1)
         IE = 1
         IXU = IE + NE
         IYU = IXU + NXU
         IXKNOT = IYU + NYU
         IYKNOT = IXKNOT + 8
         IXROW = IYKNOT + 8
         IYROW = IXROW + MX
         IWRK = IYROW + MY
         MDIST = 0
         IWGHT = 1
         WRK(IWGHT) = ONE
         CALL E01DAX(X(1),X(MX),MX,X,WRK(IWGHT),0,1,IDUM1,IDUM2,MDIST,
     *               JERROR)
         IF (JERROR.EQ.0 .AND. MDIST.EQ.MX) THEN
            CALL E01DAX(Y(1),Y(MY),MY,Y,WRK(IWGHT),0,1,IDUM1,IDUM2,
     *                  MDIST,JERROR)
            IF (JERROR.EQ.0 .AND. MDIST.EQ.MY) THEN
               CALL E01DAY(4,NXKNTS,MX,X,WRK(IWGHT),0,1,1,MX,LAMDA(5),
     *                     MX)
               CALL E01DAY(4,NYKNTS,MY,Y,WRK(IWGHT),0,1,1,MY,MU(5),MY)
               LWRK = (MX+6)*(MY+6)
               CALL E01DAZ(4,NXKNTS,X(1),X(MX),4,NYKNTS,Y(1),Y(MY),MX,X,
     *                     MY,Y,F,MY,1,MY*MX,LAMDA(5),MX,MU(5),MY,C,MY,
     *                     1,MY*MX,WRK(IXU),NXU,WRK(IYU),NYU,WRK(IE),1,
     *                     MX,NE,WRK(IXKNOT),8,WRK(IYKNOT),8,WRK(IXROW),
     *                     MX,WRK(IYROW),MY,WRK(IWRK),LWRK-IWRK+1,IER2)
               IF (IER2.NE.0) THEN
                  IERR = 3
                  NREC = 2
                  WRITE (REC,FMT=99996)
               END IF
            ELSE
               IERR = 2
               NREC = 1
               WRITE (REC,FMT=99997)
            END IF
         ELSE
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997)
         END IF
         IF (IERR.EQ.0) THEN
            LAMDA(1) = X(1)
            LAMDA(2) = X(1)
            LAMDA(3) = X(1)
            LAMDA(4) = X(1)
            LAMDA(MX+1) = X(MX)
            LAMDA(MX+2) = X(MX)
            LAMDA(MX+3) = X(MX)
            LAMDA(MX+4) = X(MX)
            MU(1) = Y(1)
            MU(2) = Y(1)
            MU(3) = Y(1)
            MU(4) = Y(1)
            MU(MY+1) = Y(MY)
            MU(MY+2) = Y(MY)
            MU(MY+3) = Y(MY)
            MU(MY+4) = Y(MY)
            PX = MX + 4
            PY = MY + 4
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
99999 FORMAT (1X,'** On entry, MX.lt.4: MX =',I16)
99998 FORMAT (1X,'** On entry, MY.lt.4: MY =',I16)
99997 FORMAT (1X,'** On entry, the X or the Y mesh points are not in s',
     *       'trictly ascending order.')
99996 FORMAT (1X,'** An intermediate set of linear equations is singul',
     *       'ar',/4X,'- the data is too ill-conditioned to compute B-',
     *       'spline coefficients.')
      END
      SUBROUTINE E01DAL(N,IBANDW,UFCTR,LUFCTR,NSETS,NTZERO,THETA,IT1,
     *                  IT2,LTHETA)
      INTEGER           IBANDW, IT1, IT2, LTHETA, LUFCTR, N, NSETS,
     *                  NTZERO
      DOUBLE PRECISION  THETA(LTHETA), UFCTR(LUFCTR)
      DOUBLE PRECISION  S
      INTEGER           I, IB, IREV, ISET, IT, ITREF, IU, IU0, IU1,
     *                  JMAX, JREV, JREVMX, JTEST, Q
      Q = N - NTZERO
      JTEST = N - IBANDW
      IU0 = -JTEST*(JTEST+1) - 2*N
      IU1 = 2*N + 1
      ITREF = 1 + (Q-1)*IT1 - IT2
      DO 80 ISET = 1, NSETS
         ITREF = ITREF + IT2
         IB = ITREF + IT1
         JMAX = Q
         I = Q + 1
         DO 60 IREV = 1, Q
            I = I - 1
            IF (I+IBANDW.LE.Q) JMAX = JMAX - 1
            IB = IB - IT1
            S = THETA(IB)
            IF (I.LE.JTEST) IU = (IBANDW-1)*(I-1) + JMAX
            IF (I.GT.JTEST) IU = (IU0+(IU1-I)*I)/2 + JMAX
            IT = ITREF - (Q-JMAX)*IT1
            IF (I.EQ.JMAX) GO TO 40
            JREVMX = JMAX - I
            DO 20 JREV = 1, JREVMX
               S = S - UFCTR(IU)*THETA(IT)
               IU = IU - 1
               IT = IT - IT1
   20       CONTINUE
   40       THETA(IT) = S/UFCTR(IU)
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
      SUBROUTINE E01DAM(NORDER,NKNOTS,XMIN,XMAX,LAMBDA,LLMBDA,JINTVL,
     *                  KNOT,LKNOT)
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           JINTVL, LKNOT, LLMBDA, NKNOTS, NORDER
      DOUBLE PRECISION  KNOT(LKNOT), LAMBDA(LLMBDA)
      INTEGER           J, K, N2
      N2 = 2*NORDER
      J = JINTVL - NORDER
      DO 20 K = 1, N2
         J = J + 1
         IF (J.LT.1) KNOT(K) = XMIN
         IF (J.GE.1 .AND. J.LE.NKNOTS) KNOT(K) = LAMBDA(J)
         IF (J.GT.NKNOTS) KNOT(K) = XMAX
   20 CONTINUE
      RETURN
      END
      SUBROUTINE E01DAN(NKNOTS,LAMBDA,LLMBDA,X,JINTVL)
      DOUBLE PRECISION  X
      INTEGER           JINTVL, LLMBDA, NKNOTS
      DOUBLE PRECISION  LAMBDA(LLMBDA)
      INTEGER           JTEMP
      EXTERNAL          E01DAP
      JTEMP = 0
      IF (NKNOTS.EQ.0) GO TO 20
      IF (X.LT.LAMBDA(1)) GO TO 20
      JTEMP = NKNOTS
      IF (X.GE.LAMBDA(NKNOTS)) GO TO 20
      JTEMP = JINTVL
      IF (JTEMP.EQ.NKNOTS) JTEMP = NKNOTS - 1
      IF (JTEMP.EQ.0) JTEMP = 1
      IF ( .NOT. (LAMBDA(JTEMP).LE.X .AND. X.LT.LAMBDA(JTEMP+1)))
     *    CALL E01DAP(NKNOTS,LAMBDA,X,JTEMP)
   20 JINTVL = JTEMP
      RETURN
      END
      SUBROUTINE E01DAP(NXDATA,XDATA,X,JINTVL)
      DOUBLE PRECISION  X
      INTEGER           JINTVL, NXDATA
      DOUBLE PRECISION  XDATA(NXDATA)
      INTEGER           I, JMAX, JMIN, JTEMP, MXDATA, POWER
      EXTERNAL          E01DAU
      JMAX = NXDATA
      JMIN = 1
      IF (X.GE.XDATA(JINTVL)) GO TO 60
      DO 20 I = 1, NXDATA
         IF (I.EQ.1) POWER = 1
         IF (I.GT.1) POWER = 2*POWER
         JTEMP = JINTVL - POWER
         IF (JTEMP.LT.1) JTEMP = 1
         IF (X.GE.XDATA(JTEMP)) GO TO 40
   20 CONTINUE
   40 JMIN = JTEMP
      JMAX = JINTVL
      IF (POWER.GT.1) JMAX = JMIN + POWER/2
      GO TO 120
   60 DO 80 I = 1, NXDATA
         IF (I.EQ.1) POWER = 1
         IF (I.GT.1) POWER = 2*POWER
         JTEMP = JINTVL + POWER
         IF (JTEMP.GT.NXDATA) JTEMP = NXDATA
         IF (X.LT.XDATA(JTEMP)) GO TO 100
   80 CONTINUE
  100 JMAX = JTEMP
      JMIN = JINTVL
      IF (POWER.GT.1) JMIN = JTEMP - POWER/2
  120 MXDATA = JMAX - JMIN + 1
      CALL E01DAU(MXDATA,XDATA(JMIN),X,JTEMP)
      JINTVL = JTEMP + JMIN - 1
      RETURN
      END
      SUBROUTINE E01DAQ(N,NSETS,NUFCTR,UFCTR,LUFCTR,THETA,IT1,IT2,
     *                  LTHETA,SFRSS,ISF1,LSFRSS,SRSS,ISR1,LSRSS)
      INTEGER           ISF1, ISR1, IT1, IT2, LSFRSS, LSRSS, LTHETA,
     *                  LUFCTR, N, NSETS, NUFCTR
      DOUBLE PRECISION  SFRSS(LSFRSS), SRSS(LSRSS), THETA(LTHETA),
     *                  UFCTR(LUFCTR)
      DOUBLE PRECISION  ONE, ZERO
      INTEGER           ISET, ITREF
      EXTERNAL          F06FBF
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
      CALL F06FBF(NUFCTR,ZERO,UFCTR,1)
      IF (NSETS.LE.0) GO TO 40
      ITREF = 1 - IT2
      DO 20 ISET = 1, NSETS
         ITREF = ITREF + IT2
         CALL F06FBF(N,ZERO,THETA(ITREF),IT1)
   20 CONTINUE
      CALL F06FBF(NSETS,ZERO,SFRSS,ISF1)
      CALL F06FBF(NSETS,ONE,SRSS,ISR1)
   40 CONTINUE
      RETURN
      END
      SUBROUTINE E01DAR(N,IBANDW,IXNZST,XROW,NSETS,YROW,IY1,LYROW,
     *                  LASTCL,UFCTR,LUFCTR,THETA,IT1,IT2,LTHETA,WRK,
     *                  LWRK)
      INTEGER           IBANDW, IT1, IT2, IXNZST, IY1, LASTCL, LTHETA,
     *                  LUFCTR, LWRK, LYROW, N, NSETS
      DOUBLE PRECISION  THETA(LTHETA), UFCTR(LUFCTR), WRK(LWRK),
     *                  XROW(N), YROW(LYROW)
      DOUBLE PRECISION  ONE, SCALE, ZERO
      INTEGER           IT, IU, IU0, IU1, J, NMBW
      EXTERNAL          DAXPY, DCOPY
      INTRINSIC         MAX, MIN
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
      NMBW = N - IBANDW
      IU0 = -NMBW**2 - 3*N + IBANDW
      IU1 = 2*N + 3
      LASTCL = MAX(LASTCL,MIN(IXNZST+IBANDW-1,N))
      CALL DCOPY(NSETS,YROW,IY1,WRK,1)
      DO 40 J = IXNZST, LASTCL
         IF (XROW(J).EQ.ZERO) GO TO 40
         IT = 1 + (J-1)*IT1
         IF (J.LE.NMBW) IU = IBANDW*(J-1) + 1
         IF (J.GT.NMBW) IU = (IU0+(IU1-J)*J)/2
         IF (UFCTR(IU).NE.ZERO) GO TO 20
         CALL DCOPY(LASTCL-J+1,XROW(J),1,UFCTR(IU),1)
         CALL DCOPY(NSETS,WRK,1,THETA(IT),IT2)
         GO TO 60
   20    SCALE = -XROW(J)/UFCTR(IU)
         CALL DAXPY(MIN(IBANDW,LASTCL-J+1),SCALE,UFCTR(IU),1,XROW(J),1)
         CALL DAXPY(NSETS,SCALE,THETA(IT),IT2,WRK,1)
   40 CONTINUE
   60 RETURN
      END
      SUBROUTINE E01DAS(N,IBANDW,UFCTR,LUFCTR,IFAIL)
      INTEGER           IBANDW, IFAIL, LUFCTR, N
      DOUBLE PRECISION  UFCTR(LUFCTR)
      DOUBLE PRECISION  ZERO
      INTEGER           IERROR, IU, J, JTEST, NP2
      INTRINSIC         MIN
      DATA              ZERO/0.0D+0/
      IERROR = 1
      NP2 = N + 2
      JTEST = NP2 - IBANDW
      IU = 1 - MIN(IBANDW,N+1)
      DO 20 J = 1, N
         IF (J.LE.JTEST) IU = IU + IBANDW
         IF (J.GT.JTEST) IU = IU + NP2 - J
         IF (UFCTR(IU).EQ.ZERO) GO TO 40
   20 CONTINUE
      IERROR = 0
   40 IFAIL = IERROR
      RETURN
      END
      SUBROUTINE E01DAT(NORDER,KNOT,LKNOT,X,BASIS)
      DOUBLE PRECISION  X
      INTEGER           LKNOT, NORDER
      DOUBLE PRECISION  BASIS(NORDER), KNOT(LKNOT)
      DOUBLE PRECISION  BNOW, BPREV, ONE, ZERO
      INTEGER           I, KL, KR
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
      IF (NORDER.EQ.1) BASIS(1) = ONE
      IF (NORDER.EQ.1) GO TO 60
   20 KL = NORDER + 1
      KR = 2*NORDER
      BNOW = ZERO
      DO 40 I = 2, NORDER
         KL = KL - 1
         KR = KR - 1
         BPREV = BNOW
         BNOW = BASIS(KL-1)/(KNOT(KR)-KNOT(KL))
         BASIS(KL) = (X-KNOT(KL))*BNOW + (KNOT(KR+1)-X)*BPREV
   40 CONTINUE
      BASIS(1) = (KNOT(NORDER+1)-X)*BNOW
   60 RETURN
      END
      SUBROUTINE E01DAU(NXDATA,XDATA,X,J)
      DOUBLE PRECISION  X
      INTEGER           J, NXDATA
      DOUBLE PRECISION  XDATA(NXDATA)
      INTEGER           JHI, JLO, MIDDLE
      LOGICAL           BELOW
      JLO = 1
      JHI = NXDATA
   20 MIDDLE = (JLO+JHI)/2
      IF (JHI-JLO.LE.1) GO TO 40
      BELOW = X .LT. XDATA(MIDDLE)
      IF (BELOW) JHI = MIDDLE
      IF ( .NOT. BELOW) JLO = MIDDLE
      GO TO 20
   40 J = JLO
      RETURN
      END
      SUBROUTINE E01DAV(NORDER,KNOT,LKNOT,X,BASIS)
      DOUBLE PRECISION  X
      INTEGER           LKNOT, NORDER
      DOUBLE PRECISION  BASIS(NORDER), KNOT(LKNOT)
      INTEGER           IKNOT, JORDER
      EXTERNAL          E01DAT
      IKNOT = NORDER + 1
      DO 20 JORDER = 1, NORDER
         IKNOT = IKNOT - 1
         CALL E01DAT(JORDER,KNOT(IKNOT),2*JORDER,X,BASIS)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE E01DAW(NORDER,XMIN,XMAX,M,X,NSETS,F,IF1,IF2,LF,NKNOTS,
     *                  LAMBDA,LLMBDA,C,IC1,IC2,LC,UFCTR,LUFCTR,KNOT,
     *                  LKNOT,XROW,LXROW,WRK,LWRK,IFAIL)
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           IC1, IC2, IF1, IF2, IFAIL, LC, LF, LKNOT,
     *                  LLMBDA, LUFCTR, LWRK, LXROW, M, NKNOTS, NORDER,
     *                  NSETS
      DOUBLE PRECISION  C(LC), F(LF), KNOT(LKNOT), LAMBDA(LLMBDA),
     *                  UFCTR(LUFCTR), WRK(LWRK), X(M), XROW(LXROW)
      INTEGER           I, IERROR, IF, JINTVL, JNTVL0, LASTCL, NUFCTR
      EXTERNAL          E01DAL, E01DAM, E01DAN, E01DAQ, E01DAR, E01DAS,
     *                  E01DAV
      NUFCTR = (NORDER*(M+NKNOTS+1))/2
      CALL E01DAQ(M,NSETS,NUFCTR,UFCTR,LUFCTR,C,IC1,IC2,LC,WRK(1),1,
     *            LWRK,WRK(1),1,LWRK)
      LASTCL = 0
      JINTVL = 0
      CALL E01DAM(NORDER,NKNOTS,XMIN,XMAX,LAMBDA,LLMBDA,JINTVL,KNOT,
     *            LKNOT)
      IF = 1 - IF1
      DO 20 I = 1, M
         IF = IF + IF1
         JNTVL0 = JINTVL
         CALL E01DAN(NKNOTS,LAMBDA,LLMBDA,X(I),JINTVL)
         IF (JINTVL.NE.JNTVL0) CALL E01DAM(NORDER,NKNOTS,XMIN,XMAX,
     *                              LAMBDA,LLMBDA,JINTVL,KNOT,LKNOT)
         CALL E01DAV(NORDER,KNOT,LKNOT,X(I),XROW(JINTVL+1))
         CALL E01DAR(M,NORDER,JINTVL+1,XROW,NSETS,F(IF),IF2,LF-IF+1,
     *               LASTCL,UFCTR,LUFCTR,C,IC1,IC2,LC,WRK,LWRK)
   20 CONTINUE
      CALL E01DAS(M,NORDER,UFCTR,LUFCTR,IERROR)
      IF (IERROR.EQ.0) CALL E01DAL(M,NORDER,UFCTR,LUFCTR,NSETS,0,C,IC1,
     *                             IC2,LC)
      IFAIL = IERROR
      RETURN
      END
      SUBROUTINE E01DAX(XMIN,XMAX,M,X,W,IW1,LW,IFNZWT,MNZWT,MDNZWT,
     *                  IFAIL)
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           IFAIL, IFNZWT, IW1, LW, M, MDNZWT, MNZWT
      DOUBLE PRECISION  W(LW), X(M)
      DOUBLE PRECISION  XPREV, ZERO
      INTEGER           IDNZWT, IERROR, INZWT, IW, R, RSTART
      DATA              ZERO/0.0D+0/
      IERROR = 0
      INZWT = 0
      IDNZWT = 0
      IW = 1 - IW1
      DO 20 R = 1, M
         IW = IW + IW1
         IF (W(IW).NE.ZERO) GO TO 40
   20 CONTINUE
      IFNZWT = M + 1
      GO TO 100
   40 RSTART = R
      IFNZWT = RSTART
      INZWT = 1
      IDNZWT = 1
      IERROR = 1
      IW = 1 + (RSTART-2)*IW1
      DO 60 R = RSTART, M
         IW = IW + IW1
         IF ((X(R).LT.XMIN .OR. X(R).GT.XMAX) .AND. W(IW).NE.ZERO)
     *       GO TO 120
   60 CONTINUE
      IERROR = 0
      IF (RSTART.EQ.M) GO TO 100
      IERROR = 0
      XPREV = X(RSTART)
      RSTART = RSTART + 1
      IW = 1 + (RSTART-2)*IW1
      DO 80 R = RSTART, M
         IW = IW + IW1
         IF (W(IW).EQ.ZERO) GO TO 80
         IF (X(R).LT.XPREV) IERROR = 2
         INZWT = INZWT + 1
         IF (X(R).GT.XPREV) IDNZWT = IDNZWT + 1
         XPREV = X(R)
   80 CONTINUE
  100 MNZWT = INZWT
      IF (IERROR.EQ.0) MDNZWT = IDNZWT
  120 IFAIL = IERROR
      RETURN
      END
      SUBROUTINE E01DAY(NORDER,NKNOTS,M,X,W,IW1,LW,IFNZWT,MDNZWT,LAMBDA,
     *                  LLMBDA)
      INTEGER           IFNZWT, IW1, LLMBDA, LW, M, MDNZWT, NKNOTS,
     *                  NORDER
      DOUBLE PRECISION  LAMBDA(LLMBDA), W(LW), X(M)
      DOUBLE PRECISION  ALPHA, C1, DIST, HALF, ONE, PPJ, RI1, RM, RQ,
     *                  SI, XI, XINEXT, ZERO
      INTEGER           I, I1, I1L, IP1, ISTRT, IW, J, K, P
      INTRINSIC         INT, MOD
      DATA              ZERO, HALF, ONE/0.0D+0, 0.5D+0, 1.0D+0/
      IF (NKNOTS.LE.0) GO TO 120
      RM = MDNZWT
      RQ = NKNOTS + NORDER
      P = NORDER/2
      IF (MOD(NORDER,2).NE.0) GO TO 20
      DIST = (RM-ONE)/(RQ-ONE)
      C1 = ONE - DIST
      GO TO 40
   20 DIST = RM/RQ
      C1 = HALF
   40 I1L = 0
      I = IFNZWT
      XINEXT = X(I)
      IW = 1 + (IFNZWT-1)*IW1
      DO 100 J = 1, NKNOTS
         PPJ = P + J
         SI = C1 + PPJ*DIST
         I1 = INT(SI)
         RI1 = I1
         ALPHA = SI - RI1
         ISTRT = I1L + 1
         I1L = I1
         DO 80 K = ISTRT, I1
            XI = XINEXT
            IP1 = I + 1
            DO 60 I = IP1, M
               XINEXT = X(I)
               IW = IW + IW1
               IF (W(IW).NE.ZERO .AND. XINEXT.NE.XI) GO TO 80
   60       CONTINUE
   80    CONTINUE
         LAMBDA(J) = (ONE-ALPHA)*XI + ALPHA*XINEXT
  100 CONTINUE
  120 CONTINUE
      RETURN
      END
      SUBROUTINE E01DAZ(NXORDR,NXKNTS,XMIN,XMAX,NYORDR,NYKNTS,YMIN,YMAX,
     *                  MX,X,MY,Y,F,IF1,IF2,LF,XLAM,LXLAM,YLAM,LYLAM,C,
     *                  IC1,IC2,LC,XUFCTR,LXU,YUFCTR,LYU,E,IE1,IE2,LE,
     *                  XKNOT,LXKNOT,YKNOT,LYKNOT,XROW,LXROW,YROW,LYROW,
     *                  WRK,LWRK,IFAIL)
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IC1, IC2, IE1, IE2, IF1, IF2, IFAIL, LC, LE, LF,
     *                  LWRK, LXKNOT, LXLAM, LXROW, LXU, LYKNOT, LYLAM,
     *                  LYROW, LYU, MX, MY, NXKNTS, NXORDR, NYKNTS,
     *                  NYORDR
      DOUBLE PRECISION  C(LC), E(LE), F(LF), WRK(LWRK), X(MX),
     *                  XKNOT(LXKNOT), XLAM(LXLAM), XROW(LXROW),
     *                  XUFCTR(LXU), Y(MY), YKNOT(LYKNOT), YLAM(LYLAM),
     *                  YROW(LYROW), YUFCTR(LYU)
      INTEGER           IERROR
      EXTERNAL          E01DAW
      CALL E01DAW(NXORDR,XMIN,XMAX,MX,X,MY,F,IF1,IF2,LF,NXKNTS,XLAM,
     *            LXLAM,E,IE1,IE2,LE,XUFCTR,LXU,XKNOT,LXKNOT,XROW,LXROW,
     *            WRK,LWRK,IERROR)
      IF (IERROR.NE.0) GO TO 20
      CALL E01DAW(NYORDR,YMIN,YMAX,MY,Y,MX,E,IE2,IE1,LE,NYKNTS,YLAM,
     *            LYLAM,C,IC2,IC1,LC,YUFCTR,LYU,YKNOT,LYKNOT,YROW,LYROW,
     *            WRK,LWRK,IERROR)
   20 IFAIL = IERROR
      RETURN
      END
      SUBROUTINE E02DFF(MX,MY,PX,PY,X,Y,LAMDA,MU,C,FF,WRK,LWRK,IWRK,
     *                  LIWRK,IFAIL)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02DFF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      INTEGER           IFAIL, LIWRK, LWRK, MX, MY, PX, PY
      DOUBLE PRECISION  C((PX-4)*(PY-4)), FF(MX*MY), LAMDA(PX), MU(PY),
     *                  WRK(LWRK), X(MX), Y(MY)
      INTEGER           IWRK(LIWRK)
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IDUMMY, IERR, JERROR, NIWRK, NREC, NWRK, NWRK1,
     *                  NWRK2, NXDIST, NXKNTS, NYDIST, NYKNTS
      CHARACTER*80      REC(2)
      INTEGER           P01ABF
      EXTERNAL          P01ABF
      EXTERNAL          E01DAX, E02DFW, E02DFZ
      INTRINSIC         MIN
      IERR = 0
      NREC = 0
      IF (PX.LT.8 .OR. PY.LT.8 .OR. MX.LT.1 .OR. MY.LT.1) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99999) PX, PY, MX, MY
      ELSE
         NXKNTS = PX - 8
         NYKNTS = PY - 8
         NWRK1 = MX*4 + PX
         NWRK2 = MY*4 + PY
         NWRK = MIN(NWRK1,NWRK2)
         IF (NWRK1.LE.NWRK2) THEN
            NIWRK = MX + PX - 4
         ELSE
            NIWRK = MY + PY - 4
         END IF
         IF (LWRK.LT.NWRK) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99998) NWRK, LWRK
         ELSE IF (LIWRK.LT.NIWRK) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997) NIWRK, LIWRK
         ELSE
            XMIN = LAMDA(4)
            XMAX = LAMDA(PX-3)
            YMIN = MU(4)
            YMAX = MU(PY-3)
            CALL E02DFW(PX,XMIN,XMAX,LAMDA,PX,JERROR)
            IF (JERROR.NE.0) THEN
               IERR = 3
               NREC = 1
               WRITE (REC,FMT=99996) 'LAMDA'
            ELSE
               CALL E02DFW(PY,YMIN,YMAX,MU,PY,JERROR)
               IF (JERROR.NE.0) THEN
                  IERR = 3
                  NREC = 1
                  WRITE (REC,FMT=99996) 'MU'
               ELSE
                  WRK(1) = ONE
                  NXDIST = 0
                  CALL E01DAX(XMIN,XMAX,MX,X,WRK(1),0,1,IDUMMY,IDUMMY,
     *                        NXDIST,JERROR)
                  IF (NXDIST.NE.MX .OR. JERROR.NE.0) THEN
                     IERR = 4
                     NREC = 2
                     WRITE (REC,FMT=99995)
                  ELSE
                     NYDIST = 0
                     CALL E01DAX(YMIN,YMAX,MY,Y,WRK(1),0,1,IDUMMY,
     *                           IDUMMY,NYDIST,JERROR)
                     IF (NYDIST.NE.MY .OR. JERROR.NE.0) THEN
                        IERR = 4
                        NREC = 2
                        WRITE (REC,FMT=99994)
                     ELSE
                        IF (NWRK1.LE.NWRK2) THEN
                           CALL E02DFZ(4,NYKNTS,YMIN,YMAX,MU,PY,4,
     *                                 NXKNTS,XMIN,XMAX,LAMDA,PX,C,1,
     *                                 PY-4,(PX-4)*(PY-4),MY,Y,MX,X,FF,
     *                                 1,MY,MY*MX,WRK,LWRK,IWRK,LIWRK)
                        ELSE
                           CALL E02DFZ(4,NXKNTS,XMIN,XMAX,LAMDA,PX,4,
     *                                 NYKNTS,YMIN,YMAX,MU,PY,C,PY-4,1,
     *                                 (PY-4)*(PX-4),MX,X,MY,Y,FF,MY,1,
     *                                 MX*MY,WRK,LWRK,IWRK,LIWRK)
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
99999 FORMAT (1X,'** On entry, either PX.lt.8, PY.lt.8, MX.lt.1 or MY.',
     *       'lt.1:',/4X,'PX =',I13,', PY =',I13,', MX =',I13,', MY =',
     *       I13,'.')
99998 FORMAT (1X,'** On entry, LWRK.lt.',I8,': LWRK =',I16)
99997 FORMAT (1X,'** On entry, LIWRK.lt.',I8,': LIWRK =',I16)
99996 FORMAT (1X,'** On entry, the knots in ',A,' are not in non-decre',
     *       'asing order.')
99995 FORMAT (1X,'** On entry, the restriction  LAMDA(4).le.X(1).lt. .',
     *       '.. .lt.X(MX).le.LAMDA(PX-3)',/4X,'is violated.')
99994 FORMAT (1X,'** On entry, the restriction  MU(4).le.Y(1).lt. ... ',
     *       '.lt.Y(MY).le.MU(PY-3)',/4X,'is violated.')
      END
      SUBROUTINE E02DFV(NKNOTS,LAMBDA,LLMBDA,XMAX,X,JINTVL)
      DOUBLE PRECISION  X, XMAX
      INTEGER           JINTVL, LLMBDA, NKNOTS
      DOUBLE PRECISION  LAMBDA(LLMBDA)
      INTEGER           JTEMP
      EXTERNAL          E01DAP
      JTEMP = 0
      IF (NKNOTS.EQ.0) GO TO 40
      IF (X.LT.LAMBDA(1)) GO TO 40
      JTEMP = NKNOTS
      IF (X.EQ.XMAX) THEN
   20    IF (LAMBDA(JTEMP).EQ.XMAX) THEN
            JTEMP = JTEMP - 1
            IF (JTEMP.GT.0) GO TO 20
         END IF
         GO TO 40
      END IF
      IF (X.GE.LAMBDA(NKNOTS)) GO TO 40
      JTEMP = JINTVL
      IF (JTEMP.EQ.NKNOTS) JTEMP = NKNOTS - 1
      IF (JTEMP.EQ.0) JTEMP = 1
      IF ( .NOT. (LAMBDA(JTEMP).LE.X .AND. X.LT.LAMBDA(JTEMP+1)))
     *    CALL E01DAP(NKNOTS,LAMBDA,X,JTEMP)
   40 JINTVL = JTEMP
      RETURN
      END
      SUBROUTINE E02DFW(NKNOTS,XMIN,XMAX,LAMBDA,LLMBDA,IFAIL)
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           IFAIL, LLMBDA, NKNOTS
      DOUBLE PRECISION  LAMBDA(LLMBDA)
      INTEGER           IERROR, J
      IERROR = 1
      IF (NKNOTS.LT.8) GO TO 60
      IF (XMAX.LE.XMIN) GO TO 60
      IF (LLMBDA.LT.NKNOTS) GO TO 60
      IF (NKNOTS.EQ.8) GO TO 40
      IF (LAMBDA(5).LE.XMIN) GO TO 60
      IF (XMAX.LE.LAMBDA(NKNOTS-4)) GO TO 60
      DO 20 J = 2, NKNOTS
         IF (LAMBDA(J).LT.LAMBDA(J-1)) GO TO 60
   20 CONTINUE
   40 IERROR = 0
   60 IFAIL = IERROR
      RETURN
      END
      SUBROUTINE E02DFX(NORDER,Q,LAMBDA,LLMBDA,XMAX,M,X,JINTVL,NINDEX,
     *                  INDEX)
      DOUBLE PRECISION  XMAX
      INTEGER           LLMBDA, M, NINDEX, NORDER, Q
      DOUBLE PRECISION  LAMBDA(LLMBDA), X(M)
      INTEGER           INDEX(Q), JINTVL(M)
      INTEGER           I, IMAX, IP, J, K, L, LAST, NKNOTS, NOW, P
      EXTERNAL          E02DFV
      INTRINSIC         MIN
      NKNOTS = Q - NORDER
      J = 0
      DO 20 I = 1, M
         CALL E02DFV(NKNOTS,LAMBDA,LLMBDA,XMAX,X(I),J)
         JINTVL(I) = J
   20 CONTINUE
      K = 0
      NOW = -NORDER
      IMAX = 0
      DO 60 I = 1, M
         LAST = NOW
         NOW = JINTVL(I)
         IF (NOW.EQ.LAST) GO TO 60
         P = MIN(NOW-LAST,NORDER)
         L = NOW + NORDER - P
         DO 40 IP = 1, P
            K = K + 1
            L = L + 1
            INDEX(K) = L
   40    CONTINUE
         IMAX = IMAX + P
   60 CONTINUE
      NINDEX = IMAX
      RETURN
      END
      SUBROUTINE E02DFY(NXORDR,NXKNTS,XMIN,XMAX,XLAM,LXLAM,NYORDR,
     *                  NYKNTS,YMIN,YMAX,YLAM,LYLAM,C,IC1,IC2,LC,MX,X,
     *                  MY,Y,S,IS1,IS2,LS,YBASIS,LYBSIS,XBASIS,YDEP,
     *                  LYDEP,IYINT,IYBSI,LIYBSI)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IC1, IC2, IS1, IS2, LC, LIYBSI, LS, LXLAM,
     *                  LYBSIS, LYDEP, LYLAM, MX, MY, NXKNTS, NXORDR,
     *                  NYKNTS, NYORDR
      DOUBLE PRECISION  C(LC), S(LS), X(MX), XBASIS(NXORDR),
     *                  XLAM(LXLAM), Y(MY), YBASIS(LYBSIS), YDEP(LYDEP),
     *                  YLAM(LYLAM)
      INTEGER           IYBSI(LIYBSI), IYINT(MY)
      DOUBLE PRECISION  T
      INTEGER           IC, IND, IS, ISREF, IX, IXBAS, IY, IYBAS, J,
     *                  JSTOP, JSTRT, JXINT, JYINT, LKNOT, NINDEX
      EXTERNAL          E01DAV, E02DFV, E02DFX
      INTRINSIC         MAX
      CALL E02DFX(NYORDR,NYKNTS+NYORDR,YLAM(5),LYLAM-NYORDR,YMAX,MY,Y,
     *            IYINT,NINDEX,IYBSI)
      LKNOT = 2*MAX(NXORDR,NYORDR)
      IYBAS = 1 - NYORDR
      DO 20 IY = 1, MY
         IYBAS = IYBAS + NYORDR
         JYINT = IYINT(IY)
         CALL E01DAV(NYORDR,YLAM(JYINT+1),LKNOT,Y(IY),YBASIS(IYBAS))
   20 CONTINUE
      ISREF = 1 - IS1 - IS2
      JXINT = 0
      DO 120 IX = 1, MX
         ISREF = ISREF + IS1
         CALL E02DFV(NXKNTS,XLAM(5),LXLAM-NXORDR,XMAX,X(IX),JXINT)
         CALL E01DAV(NXORDR,XLAM(JXINT+1),LKNOT,X(IX),XBASIS)
         DO 60 IND = 1, NINDEX
            J = IYBSI(IND)
            IC = 1 + (JXINT-1)*IC1 + (J-1)*IC2
            T = ZERO
            DO 40 IXBAS = 1, NXORDR
               IC = IC + IC1
               T = T + XBASIS(IXBAS)*C(IC)
   40       CONTINUE
            YDEP(J) = T
   60    CONTINUE
         IS = ISREF
         IYBAS = 0
         DO 100 IY = 1, MY
            IS = IS + IS2
            JSTRT = IYINT(IY) + 1
            JSTOP = IYINT(IY) + NYORDR
            T = ZERO
            DO 80 J = JSTRT, JSTOP
               IYBAS = IYBAS + 1
               T = T + YBASIS(IYBAS)*YDEP(J)
   80       CONTINUE
            S(IS) = T
  100    CONTINUE
  120 CONTINUE
      RETURN
      END
      SUBROUTINE E02DFZ(NXORDR,NXKNTS,XMIN,XMAX,XLAM,LXLAM,NYORDR,
     *                  NYKNTS,YMIN,YMAX,YLAM,LYLAM,C,IC1,IC2,LC,MX,X,
     *                  MY,Y,S,IS1,IS2,LS,WRK,LWRK,IWRK,LIWRK)
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IC1, IC2, IS1, IS2, LC, LIWRK, LS, LWRK, LXLAM,
     *                  LYLAM, MX, MY, NXKNTS, NXORDR, NYKNTS, NYORDR
      DOUBLE PRECISION  C(LC), S(LS), WRK(LWRK), X(MX), XLAM(LXLAM),
     *                  Y(MY), YLAM(LYLAM)
      INTEGER           IWRK(LIWRK)
      INTEGER           IXBAS, IYBAS, IYBSI, IYDEP, IYINT, QY
      EXTERNAL          E02DFY
      QY = NYKNTS + NYORDR
      IYBAS = 1
      IXBAS = IYBAS + MY*NYORDR
      IYDEP = IXBAS + NXORDR
      IYINT = 1
      IYBSI = IYINT + MY
      CALL E02DFY(NXORDR,NXKNTS,XMIN,XMAX,XLAM,LXLAM,NYORDR,NYKNTS,YMIN,
     *            YMAX,YLAM,LYLAM,C,IC1,IC2,LC,MX,X,MY,Y,S,IS1,IS2,LS,
     *            WRK(IYBAS),MY*NYORDR,WRK(IXBAS),WRK(IYDEP),QY,
     *            IWRK(IYINT),IWRK(IYBSI),QY)
      RETURN
      END
      SUBROUTINE F06ECF( N, ALPHA, X, INCX, Y, INCY )
      ENTRY      DAXPY ( N, ALPHA, X, INCX, Y, INCY )
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   X( * ), Y( * )
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      INTEGER            I, IX, IY
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
      SUBROUTINE F06EFF( N, X, INCX, Y, INCY )
      ENTRY      DCOPY ( N, X, INCX, Y, INCY )
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   X( * ), Y( * )
      INTEGER            I, IX, IY
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
      SUBROUTINE F06FBF( N, CONST, X, INCX )
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
      DOUBLE PRECISION   X( * )
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      INTEGER            IX
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
      CHARACTER*(*)           REC(*)
      INTEGER                 I, NERR
      CHARACTER*72            MESS
      EXTERNAL                P01ABZ, X04AAF, X04BAF
      INTRINSIC               ABS, MOD
      IF (IERROR.NE.0) THEN
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABZ
      STOP
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C
C     RETURNS  B**(1-P)  
C
      DOUBLE PRECISION X02CON
C     .. Executable Statements ..
      X02AJF = EPSILON(X02CON)
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AKF()
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
      DOUBLE PRECISION X02CON
C     .. Executable Statements ..
      X02AKF = TINY(X02CON)
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02ALF()
C
C     RETURNS  (1 - B**(-P)) * B**EMAX  (THE LARGEST POSITIVE MODEL
C     NUMBER)
C
      DOUBLE PRECISION X02CON
C     .. Executable Statements ..
      X02ALF = HUGE(X02CON)
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
C     .. Executable Statements ..
      X02AMF = 10d0**(-RANGE(1d0))
      RETURN
      END
      INTEGER FUNCTION X02BHF()
C
C     RETURNS THE MODEL PARAMETER, B.
C
C     .. Executable Statements ..
      X02BHF =     RADIX(1d0)
      RETURN
      END
      INTEGER FUNCTION X02BJF()
C
C     RETURNS THE MODEL PARAMETER, p.
C
C     .. Executable Statements ..
      X02BJF =    DIGITS(1d0)
      RETURN
      END
      INTEGER FUNCTION X02BKF()
C
C     RETURNS THE MODEL PARAMETER, EMIN.
C
C     .. Executable Statements ..
      X02BKF =  MINEXPONENT(1d0)
      RETURN
      END
      INTEGER FUNCTION X02BLF()
C
C     RETURNS THE MODEL PARAMETER, EMAX.
C
C     .. Executable Statements ..
      X02BLF =  MAXEXPONENT(1d0)
      RETURN
      END
      SUBROUTINE X04AAF(I,NERR)
      INTEGER           I, NERR
      INTEGER           NERR1
      SAVE              NERR1
      DATA              NERR1/0/
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
      INTEGER           NOUT
      CHARACTER*(*)     REC
      INTEGER           I
      INTRINSIC         LEN
      IF (NOUT.GE.0) THEN
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
99999 FORMAT (A)
      END

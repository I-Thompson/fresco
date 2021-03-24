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

*****FRXX3***************************************************************
      SUBROUTINE ERWIN(W,ECM,COEF,IEX,FORMF,NF,FORMFR,INHOMG,CH,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,NICH,N,H,M,MD,SCL,RENORM,
     $  REPEAT,WREQ,BLOCKD,AL,RE,IC,SHOW,SHSMAT,SMATEL, 
     $  CUTOFF,SKIP,FED,PTYPE,F,FIM,FIMD,LOCFIL,
     $  CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
	use io
	use factorials
	use parameters
	use drier
	use searchpar, only: final
      use fresco1,   only: rmorto,hort
      IMPLICIT REAL*8(A-H,O-Z)
C
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS
C        SOLVE   (FOR NOT BLOCKD(K))
C             (COEF(K). D2/DR2 + EN(R,K)).W(R,K)
C                  + SUM(J): COUPL(R,K,J).W(R,J) + INHOMG(R,K) = 0
C           WHERE COUPL(R,K,J) = SUM(NC): FORMF(R,JF)*CLIST(K,J,NC)
C           WHERE    JF = NFLIST(K,J,NC)
C         FOR K & J <= IEX (ALL OTHER COUPLINGS MUST BE ITERATED)
C       &   WHERE EN(R,K) IS THE DIAGONAL PART OF COUPL
C
C     ASSUMED FOR ARRAYS HERE THAT N<=MAXN AND NEQS <= MAXCH
C
C  Using 'Enhanced Numerov' of Thorlacius & Cooper (JCP 72(1987) 70)
C   with 5 terms in cosh(sqrt(12T)) expansion, but only diagonal potl.
C
      PARAMETER(MINVEC=5)
C               MINIMUM VECTOR LENGTH FOR EFFICIENT OPERATION.
      SAVE IS
      COMPLEX*16 F(NFDEC),WVD(MAXN),FIMD(MAXB,MAXCH),FIM(MAXB,MAXCH)
      COMPLEX*16 MAT(2*NEQS+1,2*NEQS+2),S,ZIV,ZPV,FI(MAXB,MAXCH),
     &           V(MAXB,MAXCH),ZI(MAXB,MAXCH),ZM(MAXB,MAXCH),CZ,
     &           SE,CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI
      COMPLEX*16 FORMF(MAXM,NF),INHOMG(MAXN,MAXCH),F8,FORMFR(NF),
     &       C,W(MAXN,NICH),CLIST(MAXCH,MAXCH,MCLIST),COUPL(MAXB,MAXB)
      COMPLEX*16 ZDL,ZHPLUS,ZHMINUS
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),FIT,CUTOFF,SHOW,SKIPL,SKIPU,FEDL,FEDU,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),HOLES,SKFED,PTYPE(8,NF)
      LOGICAL WREQ,SING,REPEAT,SHSMAT,SMATEL,BLOCKD(MAXCH),ITV,
     X           SKIP(MAXCH),FED(MAXCH),DO1,DO2,DO3,BFEED,BSKIP,LOCFIL,
     X           FCWFN,FJSWTCH
      AMD1(C) = ABS(DBLE(C)) + ABS(AIMAG(C))
      complex*16 wdiag(neqs) ! AMM

C
      NM1 = N-1
      NR = 1 + 2*NEQS
      NP = NR + 1
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
      ZI2 = (0.0D0,0.5D0)
      CI = (0.0D0,1.0d0)
      CZ = (0D0,0D0)
      ONEC = (1D0,0D0)
C     AL = 0.0
      ALUN = 0.0
	WVD(:) = CZ
      NQ1 = 1+NEQS
      IEX1= IEX + 1
      IIX = IEX
      NITS = 2+IEX
      FIT = 2
!     NLIF = NICH*2 + IEX*NICH
      NLIF = NITS*NICH
!      NVREQ = 2 * NN*NLIF  + 4*MAXB*MAXCH

! AMM !!!!!!!!!!!!!!!!!!
      norto=nint(hort/h(1))
!      write(*,*)'hort,norto,rmorto=',hort,norto,rmorto
!      write(*,'(3x,"R-turns:",100(5f8.2))')rturnv(:)


       IF(REPEAT) THEN
         NITS = 1
         FIT = 1
       ENDIF
      IF(.NOT.REPEAT .AND. IEX.GE.NEQS) FIT = 3

      DO1 = REPEAT
      DO2 = .NOT.REPEAT .AND. IEX.LT.NEQS
      DO3 = .NOT.REPEAT .AND. IEX.GT.0
      IF(.NOT.REPEAT) IS =  CUTOFF
      IF(.NOT.REPEAT) IS = MIN(1*N/2, MAX(2, IS ))
       if(SHOW.ge.3) 
     *       write(6,*) 'ERWIN: CUTOFF,IS =',CUTOFF,IS
       if(SHOW.ge.3) 
     *       write(6,*) 'ERWIN: CUTVAL =',(CUTVAL(I),I=1,NEQS)
      NN = NM1 - IS + 1
      TMAX = 20.
      TMIN = -125.
        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS
      NVREQ =  NN*NLIF
!      write(48,*) 'NICH,IEX,NITS,NLIF,NN =',NICH,IEX,NITS,NLIF,NN
!#      write(48,*)  'NN,NLIF,NVREQ,NFDEC',NN,NLIF,NVREQ,NFDEC
!#      call flush(48)
      CALL CHECK(NVREQ,NFDEC,17)
      if(SHOW>1) then
	do K=1,min(3,NEQS)
	write(KO,*) 'Ch ',K,' ECM/COEF=',real(ECM(K)/COEF(K))
	enddo
      endif

C
      IF (FJSWTCH)  THEN
c                              skip nmrov , match CRCWFN to zero
c
761    DO 765 I=1,NR
       DO 765 J=1,NP
765       MAT(I,J) = 0
c
         MAT(1,1) = 1
         DO 780 IT=1,NEQS
               MAT(IT+NQ1,IT+1) = ONEC
               MAT(IT+1,IT+1)   = ONEC*CRCRAT(it)
            MAT(1 ,IT+1) = 1
            MAT(1 ,IT+NQ1) = 0
780       CONTINUE
C
      ELSE
C
      IF(SHOW.GE.2) WRITE(KO, 5) IEX,NEQS,NF,M,REPEAT,NICH,NQ1,NR,
     X                          NVREQ
      IF(SHOW.GE.2) WRITE(KO, 4) (H(K),K=1,NEQS)
4     FORMAT(' ERWIN step sizes are',10F8.5)
5     FORMAT(' ERWIN given',4I4,L4,' & requires',3I6,' i.e.',I8,'/',I8)
      CALL CHECK(NITS,MAXB,15)
C
C   If SKIP(K), don't integrate channel K, as same as last iteration!
C   If not FED(K) too, set channel to zero to start with.
C   NB. Some channels are non-zero if FED(K) & SKIP(K).
C
C    KFIRST - KLAST (incl) = range of channels integrated
C    FEDl   - FEDL  (incl) = range of channels non-zero
C     (may be larger than KFIRST->KLAST, if some SKIPed)
C
       BFEED = .FALSE.
       BSKIP = .TRUE.
       DO 5001 K=1,IEX
        BSKIP = BSKIP .AND. SKIP(K)
5001    IF(FED(K)) BFEED = .TRUE.
C     in a closely-coupled set,
C     all FED if any one is,
C and only allow any SKIP if all are.
C
       DO 5002 K=1,IEX
       IF(BFEED) FED(K) = .TRUE.
5002   IF(.NOT.BSKIP) SKIP(K) = .FALSE.
      FEDU = 0
      FEDL = NEQS+1
      SKIPU = 0
      SKIPL = NEQS+1
          DO 7 K=1,NEQS
            SMAT(K) = CZ
             IF(.NOT.REPEAT) THEN
               SKIP(K) = .FALSE.
               FED(K)  = .FALSE.
             ELSE
              FED(K) = FED(K) .AND..NOT. BLOCKD(K)
              SKIP(K) = SKIP(K) .OR. .NOT.FED(K)
               IF(FED(K)) FEDL = MIN(FEDL,K)
               IF(FED(K)) FEDU = MAX(FEDU,K)
               IF(.NOT.SKIP(K)) SKIPL = MIN(SKIPL,K)
               IF(.NOT.SKIP(K)) SKIPU = MAX(SKIPU,K)
             ENDIF
          H2C(K) = H(K)**2 / COEF(K)
C7        LL1(K) = -L(K)*(L(K)+1)/H2C(K)
C7        LL1(K) = -LL1(K)/H2C(K)
7        CONTINUE
         HOLES = 0
         SKFED = 0
      IF(REPEAT) THEN
          KFIRST = MAX(FEDL,SKIPL)
          KLAST  = MIN(FEDU,SKIPU)
            DO 701 K=KFIRST+1,KLAST-1
             IF(SKIP(K).AND.FED(K))  SKFED = SKFED + 1
701          IF(SKIP(K))  HOLES = HOLES + 1
            DO 702 K=FEDL,FEDU
702          IF(SKIP(K).AND.FED(K))  SKFED = SKFED + 1
       ELSE
           FEDL   = 1
           FEDU   = NEQS
           KFIRST = 1
           KLAST  = NEQS
         ENDIF
        ITV = IIX-KFIRST+1 .GE. MINVEC
C
         IF(SHOW.GE.2.OR.KFIRST.GT.1.AND.KFIRST.LE.IEX) THEN
          WRITE(KO,*) 'FEDL,FEDU,SKIPL,SKIPU,KFIRST,KLAST,HOLES,SKFED ='
     X            ,FEDL,FEDU,SKIPL,SKIPU,KFIRST,KLAST,HOLES,SKFED
            WRITE(KO,*) 'FED(*) =',(FED(K),K=1,NEQS)
            WRITE(KO,*) 'SKIP(*) =',(SKIP(K),K=1,NEQS)
            ENDIF
         DO 705 I=1,N
705           MAGN(I) = 1.0
C
      DO 11 IT=1,NITS
           KMAX = NEQS
           IF(IT.GE.3) KMAX = IEX
          KIMAX = MIN(KMAX,NICH)
           L1 = (IT-1) * NICH
!           IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
          DO 10 K=1,NEQS
            IF(WREQ.AND.K.LE.KIMAX.AND..NOT.(SKIP(K).AND.FED(K))) THEN
C                            THIS VECTORISES OK.
            II = L1+K - IS*NLIF
             DO 8 I=IS,NM1
8           F(II+I*NLIF) = CZ
            ENDIF
         V(IT,K)   = 0
         FI(IT,K) = 0
         ZM(IT,K)  = 0
         ZI(IT,K)  = 0
          IF(REPEAT) GO TO 10
            FIM(IT,K) = 0
            FIMD(IT,K)= 0
10       CONTINUE
11     CONTINUE

C
        DO 141 K=1,FEDU
141      IF(SHOW.GE.6) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
142      FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------
          IF(REPEAT .AND. KLAST.LT.KFIRST) GO TO 61
	  if(LOCFIL) then
	  rewind 19
	  do I=1,IS-1
	   read(19) 
	  enddo
	  endif
         DO 60 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	  if(LOCFIL) read(19) FORMFR
          DO 148 K=KFIRST,KLAST
           IF(SKIP(K)) GO TO 148
            C = CZ
	   if(LOCFIL) then
CDIR$                NOVECTOR
            DO NC=1,NCLIST(K,K)
		JF = NFLIST(K,K,NC)
            IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMFR(JF)
	    enddo
	   else
CDIR$                NOVECTOR
            DO NC=1,NCLIST(K,K)
		JF = NFLIST(K,K,NC)
            IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMF(I,JF)
	    enddo
	   endif
CDIR$ VECTOR
C          SMAT(K) = (LL1(K)*RI2 - ECM(K) + C) * H2C(K)
           SMAT(K) = -LL1(K)*RI2 + (-ECM(K) + C) * H2C(K)
            F8 = SMAT(K)
            IF(ABS(AIMAG(F8)).LT.1E-20*ABS(DBLE(F8)))
     X         SMAT(K) = DBLE(SMAT(K))
148        CONTINUE
C
       if(DO2) then
         DO 12 K=IEX1,NEQS
           IF(SKIP(K)) GO TO 12
           IF(I.ne.CUTVAL(K)) GO TO 12
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
               ZI(2,K) = SCL * EXP(T)
       	IF(SHOW.GE.3) write(KO,1302) I,2,K,ZI(2,K)
12       CONTINUE
       endif
       if(DO3) then
          DO 13 IT=3,NITS
            K = IT-2
              IF(SKIP(K)) GO TO 13
              IF(I.ne.CUTVAL(K)) GO TO 13
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT,K) = SCL * EXP(T)
       	IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT,K)
13         CONTINUE
       endif
1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)
         DO 16 IT=FIT,NITS
            IF(SKFED.GT.0) THEN
C       Approx. ZI(IT,J), if J is FED and SKIPed,
C         by F(L1+J,I) already stored, for use in DO 24 loop.
              KMAX = FEDU
              IF(IT.GE.3) KMAX = MIN(IEX,FEDU)
              KIMAX = MIN(KMAX,NICH)
               L1 = (IT-1) * NICH
!                IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
                DO 157 K=FEDL,KIMAX
 157             IF(SKIP(K).AND.FED(K)) ZI(IT,K) = F(L1+K+(I-IS)*NLIF)
             ENDIF
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
        IF(KMAX-KFIRST.LT.MINVEC-1) THEN
c AMM.................................................................
       WDIAG(:)=0.
       DO K=1,NEQS
       wdiag(k)=(ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))
       ENDDO
c.......................................................................

CDIR$ NOVECTOR
         DO 158 K=KFIRST,KMAX
158       IF(.NOT.SKIP(K)) FI(IT,K) = ZI(IT,K) *
C    X      (ONEC - SMAT(K) * R12 + ENA2*SMAT(K)**2 + ENA3*SMAT(K)**3)
     X      (ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))
CDIR$ VECTOR
        ELSE
         DO 159 K=KFIRST,KMAX
159       IF(.NOT.SKIP(K)) FI(IT,K) = ZI(IT,K) *
     X      (ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))
        ENDIF
16      CONTINUE
         IF(SHOW.GE.7) THEN
          DO 161 K=KFIRST,KLAST
161          WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    SMAT(K)/H2C(K)+ECM(K),FI(FIT,K)
165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,2E10.1)
          ENDIF
C
         IF(DRY) GO TO 40
         IF(KFIRST.GT.IEX) GO TO 28
         DO 25 K=KFIRST,IEX
            DO 17 J=FEDL,IIX
17          COUPL(K,J) = 0.0
         IF(SKIP(K)) GO TO 25
            H2 = H2C(K) * R12
           IF(I.LT.IC) GO TO 21
           DO 19 J=FEDL,IIX
            C = CZ
	   if(LOCFIL) then
CDIR$                NOVECTOR
            DO 18 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
18          IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMFR(JF)
	   else
CDIR$                NOVECTOR
            DO 181 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
181         IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMF(I,JF)
	   endif
CDIR$ VECTOR
19         COUPL(K,J) = C
C (THE DIAGONAL PART COUPL(K,K) IS NOT USED, AS CCU(K,K)=0 ALWAYS
C
21       DO 24 J=FEDL,IIX
            IF(K==J.or.NCLIST(K,J)==0) GO TO 24
            C = COUPL(K,J) * H2
c    if(i.gt.200.and.i.lt.220) write(6,*) i,k,j,c
          IF(DO2) FI(2,K) = FI(2,K) - C * ZI(2,J)
          IF(DO1) FI(1,K) = FI(1,K) - C * ZI(1,J)
            IF(DO3) THEN
            IF(.NOT.ITV) THEN
CDIR$                      NOVECTOR
             DO 22 IT=3,NITS
 22          FI(IT,K) = FI(IT,K) - C * ZI(IT,J)
CDIR$ VECTOR
            ELSE
             DO 23 IT=3,NITS
 23          FI(IT,K) = FI(IT,K) - C * ZI(IT,J)
            ENDIF
            ENDIF
24       CONTINUE
25    CONTINUE
C
28     IF(.NOT.REPEAT) THEN
         DO 29 IT=FIT,NITS
         DO 29 K=KFIRST,KLAST
29       V(IT,K) = 0.0
       ELSE
C           REPEAT: IT = 1 ONLY
         IF(KLAST-KFIRST.LT.MINVEC-1) THEN
CDIR$             NOVECTOR
         DO 314 K=KFIRST,KLAST
          IF(.NOT.SKIP(K)) V(1,K) = INHOMG(I,K) * H2C(K)
314      FI(1,K) = FI(1,K) - V(1,K) * R12
CDIR$ VECTOR
          ELSE
         DO 315 K=KFIRST,KLAST
C                                 THIS STILL VECTORISES!
          IF(.NOT.SKIP(K)) V(1,K) = INHOMG(I,K) * H2C(K)
315      FI(1,K) = FI(1,K) - V(1,K) * R12
         ENDIF
        ENDIF
         IF(KFIRST.GT.IEX) GO TO 40
         DO 35 K=KFIRST,IEX
         IF(SKIP(K)) GO TO 35
         DO 34 J=FEDL,IIX
           IF(K==J.or.NCLIST(K,J)==0) GO TO 34
            C = COUPL(K,J) * H2C(K)
          IF(DO2) V(2,K) = V(2,K) + C * FI(2,J)
          IF(DO1) V(1,K) = V(1,K) + C * FI(1,J)
          IF(DO3) THEN
            IF(.NOT.ITV) THEN
CDIR$                      NOVECTOR
            DO 33 IT=3,NITS
33          V(IT,K) = V(IT,K) + C * FI(IT,J)
CDIR$ VECTOR
            ELSE
            DO 331 IT=3,NITS
331          V(IT,K) = V(IT,K) + C * FI(IT,J)
            ENDIF
           ENDIF
34       CONTINUE
35       CONTINUE
40       DO 46 IT=FIT,NITS
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
         K = KMAX-KFIRST+1
         DO 43 K=KFIRST,KMAX
            IF(SKIP(K)) GO TO 43
            ZIV = ZI(IT,K)
            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
            ZM(IT,K) = ZIV
            ZI(IT,K) = ZPV
43            CONTINUE
46    CONTINUE

c AMM ........................................................
      r=h(1)*(i-1.)
      if ((r.lt.rmorto).and.(i.gt.is).and.(hort.gt.0)
     & .and.(mod(i-is,norto).eq.0).and.(i.lt.md)) then 
        write(*,*) '    -> re-orthogonalizing at i,r=',i,r
      CALL LQFACT(ZI,ZM,FI,WDIAG,COUPL,R12,NEQS,MAXB,MAXCH)
      endif
c ............................................................

      KMAX = KLAST
      DO 468 IT=FIT,NITS
       IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
       IF(I.EQ.M) THEN
        DO 463 K=KFIRST,KMAX
463     IF(.NOT.SKIP(K)) FIM(IT,K) = FI(IT,K)
       ELSE IF(I.EQ.MD) THEN
        DO 465 K=KFIRST,KMAX
465     IF(.NOT.SKIP(K)) FIMD(IT,K) = FI(IT,K)
       ENDIF
      IF(.NOT.WREQ) GO TO 468
       KIMAX = MIN(KMAX,NICH)
        L1 = (IT-1) * NICH
!        IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
        DO 466 K=KFIRST,KIMAX
466     IF(.NOT.SKIP(K)) F(L1+K + (I-IS)*NLIF) = FI(IT,K)
468   CONTINUE
            T = 0.0
         IF(REPEAT)  GO TO 60
         DO 47 IT=FIT,NITS
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
            DO 47 K=KFIRST,KMAX
            F8 = FI(IT,K)
            T = MAX(T, AMD1(F8))
   47       IF(T .GT. BIG ) GO TO 48
            MAGN(I) = T
         IF(I.NE.NM1) GO TO 60
 48      IF(T.EQ.0.0) GO TO 60
            MAGN(I) = T
         T = RENORM/T
         IF(SHOW.GE.1) WRITE(KO,49) T,I,RENORM
49       FORMAT(' Renormalising by',E12.4,' at',I4,' to max. WF.',E12.4)
         if(t.eq.0d0) then
          write(KO,*) 'RENORMALISING TO 0! at I=',I,' as lWF=',MAGN(I)
          write(KO,*) ' FI:',FI
          T = SMALL
          endif
         FMIN = small/T
         DO 55 IT=FIT,NITS
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
          L1 = (IT-1) * NICH
!           IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
         DO 51 K=KFIRST,KMAX
            IF(BLOCKD(K)) GO TO 51
            ZI(IT,K) = ZI(IT,K) * T
            ZM(IT,K) = ZM(IT,K) * T
            FIM(IT,K) = FIM(IT,K) * T
            FIMD(IT,K) = FIMD(IT,K) * T
           IF(.NOT.WREQ.OR.K.GT.NICH) GO TO 51
            DO 50 J=IS,I
             JJ = L1+K+(J-IS)*NLIF
             F8 = F(JJ)
50          IF(AMD1(F8).GT.FMIN) F(JJ) = F(JJ) * T
51       CONTINUE
55       CONTINUE
C
60       CONTINUE
C-----------------------------------------------------------------------
c
61    DO 65 I=1,NR
      DO 65 J=1,NP
65     MAT(I,J) = 0
         DO 80 IT=1,NQ1
         IF(IT.EQ.1 .AND..NOT.REPEAT) GO TO 75
         J = IT-1
         IF(J.GT.IEX) I = 2
         IF(J.LE.IEX) I = 2 + J
         IF(IT.EQ.1) I = 1
           KMAX = NEQS
           IF(I.GE.3) KMAX = IEX
         DO 70 K=1,KMAX
            IF(BLOCKD(K)) GO TO 70
            IF(IT.EQ.1.AND..NOT.FED(K)) GO TO 70
         IF(I.EQ.2 .AND. K.NE.IT-1) GO TO 70
            MAT(K+NQ1,IT) = FIMD(I,K)
            MAT(K+1,IT)   = FIM(I,K)
70          CONTINUE
75       MAT(1 ,IT) = 1
         IF(IT .GT. 1 .AND. REPEAT) MAT(1,IT) = 0
C           IF(IT.GT.K) GO TO 80
         MAT(1 ,IT+NQ1) = 0
80       CONTINUE
C
c
      ENDIF
c
C                    MATCH TO REQUIRED ASYMPTOTIC CONDITIONS
C
      IF (FCWFN) THEN
C                           match Numerov to CRCWFN
CDIR$   NOVECTOR
       DO 82 K=1,NEQS
       DO 82 K2=1,NEQS
           IF(SHOW.ge.4.and.ABS(CGMAT(K2,K,2)).gt.0.0) 
     *       WRITE(6,716) K,K2,CGMAT(K2,K,2),CFMAT(K2,K,2),
     *           CGMAT(K2,K,1),CFMAT(K2,K,1)
716            FORMAT(1X,2I3,': CRCWFN= ',4(D15.7,3X))
c
c james
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,2)+CI*CFMAT(K2,K,2)
         MAT(K2+1,K+NQ1)=   ZI2*ZDL*ZHPLUS
         ZHPLUS = CGMAT(K2,K,1)+CI*CFMAT(K2,K,1)
         MAT(K2+NQ1,K+NQ1)= ZI2*ZDL*ZHPLUS
c
c        MAT(K2+1,K+NQ1)=   ZI2*(CGMAT(K2,K,2)+CI*CFMAT(K2,K,2))
c        MAT(K2+NQ1,K+NQ1)= ZI2*(CGMAT(K2,K,1)+CI*CFMAT(K2,K,1))
82       CONTINUE
CDIR$ VECTOR
         MAT(1,NP)=1
         DO 85 K=1,NEQS
c james
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,2)-CI*CFMAT(K,EL,2)
         MAT(K+1,NP)=   ZI2*ZDL*ZHMINUS
         ZHMINUS = CGMAT(K,EL,1)-CI*CFMAT(K,EL,1)
         MAT(K+NQ1,NP)= ZI2*ZDL*ZHMINUS
c
c          MAT(K+1,NP)=   ZI2*(CGMAT(K,EL,2)-CI*CFMAT(K,EL,2))
c          MAT(K+NQ1,NP)= ZI2*(CGMAT(K,EL,1)-CI*CFMAT(K,EL,1))
85       CONTINUE
      ELSE
c                         match Numerov to Uncoupled Coulomb wfns
       DO 90  K=1,NEQS
         MAT(K+1,K+NQ1) = CH(1,K)
         MAT(K+NQ1 ,K+NQ1) = CH(2,K)
C                                      C(K) == MAT(K,NP) IS THE RHS.
         MAT(K+1,NP) = 0
         MAT(K+NQ1 ,NP) = 0
90       CONTINUE
         MAT(1 ,NP) = 1
c                                         note minus as ch = i/2*H+
         MAT(EL+1,NP) = - CONJG(CH(1,EL))
         MAT(EL+NQ1 ,NP) = - CONJG(CH(2,EL))
      ENDIF
         DO 100 I=1,NR
100        IF (SHOW.GE.5) WRITE(KO,105) I,(MAT(I,J),J=1,NP)
105      FORMAT(' MAT(',I2,',*) =', /(6(1P,E10.2,E9.1,1X)))
C
      CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
         IF(SING) GO TO 600
C
      T = EXP(DBLE(S) * 0.5/NEQS)
      IF(.NOT.REPEAT) RENORM = RENORM / T
       IF(RENORM.GT.BIG) RENORM = BIG
       IF(RENORM.LT.EPS) RENORM = EPS
      IF(SHOW.GE.3) WRITE(KO,130) (MAT(I,NP),I=1,NR)
130   FORMAT(/' The combining coefficients are',1P,/(2(E20.8,E15.8)))
C
c     ENDIF
c                  end of FJSWITCH loop matching CRCWFN to zero
      DO 150 K=1,NEQS
         DO 138 I=1,N
138         WVD(I) = CZ
      SMAT(K) = MAT(K+NQ1,NP)
C
      IF(ABS(SMAT(K)).GT.1E-10) THEN
        T = 0.0
       DO 145 IT=1,NQ1
         J = IT-1
         IF(J.GT.IEX) I = 2
         IF(J.LE.IEX) I = 2 + J
         IF(IT.EQ.1) I = 1
           KMAX = NEQS
           IF(I.GE.3) KMAX = IEX
         KIMAX = MIN(KMAX,NICH)
          L1 = (I-1) * NICH
!           IF(I.GE.3) L1 = (I-3)*NICH + 2*NICH
         IF(I.EQ.2 .AND. K.NE.IT-1 .OR. MAT(IT,NP).EQ.CZ
     &                 .OR. K.GT.KMAX) GO TO 145
           IF(WREQ.AND.K.LE.KIMAX) THEN
                II = L1+K - IS*NLIF
         DO 140 J=IS,NM1
140         WVD(J) = WVD(J) + MAT(IT,NP)*F(II+J*NLIF)
          ELSE
            WVD(MD) = WVD(MD) + MAT(IT,NP)*FIMD(I,K)
          ENDIF
            T = MAX(T,ABS(WVD(MD)))
145      CONTINUE
             H2C(K) = T / (ABS(WVD(MD))+1D-20)
	   if(ECM(K)>0.0) then
             AL = MAX(AL, H2C(K))
             ALUN = MAX(ALUN, H2C(K)*ABS(SMAT(K)))
           endif
!	     write(106,*) ALUN,K,WVD(MD),H2C(k),smat(K),AL
      ENDIF
         IF(.NOT.WREQ.OR.K.GT.NICH)  GO TO 150
         DO 149 I=1,NM1
149         W(I,K) = WVD(I)
150      CONTINUE
      SING = SHOW.GE.1.AND..NOT.REPEAT
      IF(SING) WRITE(KO,1157) S,AL
1157  FORMAT(' ERWIN: LOG10(DETERMINANT) =',2F11.3,',  ACCURACY LOSS ='
     X  , 1P,E10.2)
C                       NOW CHECK THE LOGARITHMIC DERIVS. MATCH REQD.
      DO 200 K=1,NEQS
        IF(SHOW.LT.3.OR..NOT.WREQ.OR.K.GT.NICH) GO TO 185
         IF(ABS(W(M,K)) .EQ. 0) W(M,K) = (1.0E-20,0.0)
        SE = W(MD,K) / W(M,K)
        INIEX = 0
        IF(K.EQ.EL) INIEX = 1
        S  =            ( -INIEX*CONJG(CH(2,K)) - SMAT(K)*CH(2,K))
     &                 /( -INIEX*CONJG(CH(1,K)) - SMAT(K)*CH(1,K)+1D-20)
      C  = S - SE
      WRITE(KO,180) SE,C
180   FORMAT(/'  ',2F10.6,' matched within ',2E12.4)
c185   IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).AND.ABS(SMAT(K)).GT.1.D-6)
185   IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).and.final  )
     & WRITE(KO,190) K,SMAT(K),L(K),JVAL(K),CORESP(K),L(EL),
     &           LOG10(MAX(H2C(K),1D0))
190   FORMAT( '   S-matrix',I5,' = '
     & , 2F10.5,' for L=',I5,', J=',F7.1,' channel on core I =',F6.1,
     & ' from L=',I5,',  Acc. loss =',F5.1,' D.')
      IF(K.EQ.EL.AND.SMATEL.AND.ABS(SMAT(K)).GT.1.D-6.and.final) THEN
         SE = LOG(SMAT(EL))
         WRITE(KO,195) K,SE*(0.,-.5)*180./3.14159,L(K),JVAL(K)
195      FORMAT( ' Elastic phase shift ',I3,' = '
     &     , 2F8.3,' deg. for the L =',I5,', J =',F7.1,' channel.')
         WRITE(45,196) ECM(EL),SE*(0.,-.5)*180./3.14159,L(K),JVAL(K)
196        format(f10.3,2f9.3,' for LJin =',i6,f6.1)
	 written(45) = .true.
         ENDIF
200   CONTINUE
      if(.not.FJSWTCH) then
C     IF(AL.GT.1D10.OR.AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
      IF(AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
     X          3.*ALUN*RE * 100.
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWIN =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(ALUN*RE.GT..03) then
  		WRITE(KO,*) AL,ALUN,RE
		write(KO,*) 'Increase lower cutoff: eg CUTL. STOP'
		stop
	  endif
	endif
      RETURN
600   IF(REPEAT) CALL ABEND(32)
C600   IF(REPEAT) STOP
CDIR$   NOVECTOR
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
CDIR$ VECTOR
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
C     STOP
      CALL ABEND(4)
      END
      SUBROUTINE KERNEL(KCOEF,NLL,ALOSS,ALOSSM,RER,
     &                  C1,C2,IC1,IC2,REPEAT, LVAL,ICTO,ICFROM,REV,
     &                  PART,DNL,A,B,NK,FPT,FI,CHNO,CUTOFF,CP,KIND,NIB,
     &                  NUMLT,LTMIN,P,Q,LTRANS,QNF,IC7,HF,HT,FFR)
	use parameters
	use io
	use factorials
	use drier
	use kcom
	use kcom2
	use trace
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KCOEF(MAXQRN,NUMLT),ALOSS(MAXQRN),DNL(NLO)
      INTEGER LVAL(MAXCH),D,C,PART(MAXCH),FPT(2,NK),FI,
     &        CP,QNF(19,MSP),CUTOFF,CHNO(MFNL,6),C1,C2
      LOGICAL REV,PRES,ODD,THERE,FFR,REPEAT,C1FR,LTRANS(MAXQRN)
      DATA ONE,TWO /1D0, 2D0/
      ODD(I) = I.NE.(I/2)*2
C
      IF(REPEAT) GO TO 7
      if(FFR) then
	if(allocated(QRLN)) deallocate(QRLN,FNL)
	allocate (QRLN(NLO,NLL,NIB),FNL(NLL,NLO))
	else
	if(allocated(QERN)) deallocate(QERN,FNC)
      	allocate (QERN(NLO,NLL,NIB),FNC(NLL,NLO))
	endif
	if(allocated(WHERE)) deallocate(WHERE,WHOL,WHOI)
	allocate (WHERE(MAXMUL,MAXQRN),WHOL(NIB),WHOI(NIB))
      Z = 0.0
      EPS = 1E-14
      DO 1    J=1,NLO
1     DNL(J)  = (J - NLC - ONE) * MLT * HF
      HFHT = HF/HT
      RINTO = HT * MR
      MAXL1 = MAXL + 1
      MINL1 = MINL + 1
      MINLRQ = MAXL1
      MAXLRQ = -1
      DO 5 IB=1,NIB
      WHOL(IB) = 0
5     WHOI(IB) = 0
      DO 6 IN=1,NK
      DO 6 L1=MINL1,MAXL1
6     WHERE(L1,IN) = 0
      REPEAT = .TRUE.
      IB = 0
7     C1FR = ICFROM.EQ.IC1
      IF(C1FR) THEN
         C = C1
         D = C2
      ELSE
         C = C2
         D = C1
      ENDIF
C     IF(ICTO.NE.PART(D).OR.ICFROM.NE.PART(C)) STOP 7
      IF(ICTO.NE.PART(D).OR.ICFROM.NE.PART(C)) CALL ABEND(32)
C       SO TO 'D' FROM 'C' CHANNEL
C  AND WITH NAME(1,IC1) LIKE D & NAME(1,IC2) LIKE P IN (D,P) REACTION
      LD = LVAL(D)
      LC = LVAL(C)
      R7 = SQRT((TWO*LD+ONE)*(TWO*LC+ONE))
      THERE = .FALSE.
      DO 8 IN=1,NK
      DO 8 ILT=1,NUMLT
 8     IF(ABS(KCOEF(IN,ILT)).GT.EPS.AND.LTRANS(IN)) THERE = .TRUE.
      IF(.NOT.THERE) GO TO 100
C
      IF(NL.EQ.0) REWIND 12
      NL = NL + 1
      DRY = DRY .OR. NL.GT.MFNL
9     FORMAT(//' ****** NOT ENOUGH ROOM FOR',I4,' NON-LOCAL FORM',
     &' FACTORS  IN MFNL ARRAY OF',I4,' *****')
	if(FFR.and..not.allocated(FNL)) write(6,*) 'FNL fail!'
	if(.not.FFR.and..not.allocated(FNC)) write(6,*) 'FNC fail!'
      DO 1111 I=1,NLL
         IF(FFR) THEN
            DO 10 J=1,NLO
10            FNL(I,J) = 0
      ELSE
             DO 11 J=1,NLO
11                 FNC(I,J) = 0
      ENDIF
1111   CONTINUE
      CHNO(NL,1) = D
      IF(.NOT.REV) CHNO(NL,1) = -D
      CHNO(NL,2) = C
      CHNO(NL,3) = 2
      IF(KIND.EQ.8) CHNO(NL,3) = IC7
      IF(.NOT.FFR) CHNO(NL,3) = CHNO(NL,3) + 10
      CHNO(NL,4) = NLL
      CHNO(NL,6) = 0
      PRES = .FALSE.
C
      DO 80 IN=1,NK
       IF(.NOT.LTRANS(IN)) GO TO 80
C
C   HERE, TO START :    FPT(1  IS NO. OF B.S. IN CHANNEL 'C1' ("DEUT")
C                   AND FPT(2  IS NO. OF B.S. IN CHANNEL 'C2' ("PROT")
         IF(C1FR) THEN
           KNP = FPT(1,IN)
           KN = FPT(2,IN)
         ELSE
           KNP = FPT(2,IN)
           KN  = FPT(1,IN)
         ENDIF
C
C   NOW, MORE PROPERLY, KNP IS NO. OF B.S. IN CHANNEL 'C' (FROM)
C                   AND KN  IS NO. OF B.S. IN CHANNEL 'D' (TO)
         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
         VINT = 0.0
         THERE = .FALSE.
         DO 16 ILT=1,NUMLT
 16        IF(ABS(KCOEF(IN,ILT)).GT.EPS) THERE = .TRUE.
         IF(.NOT.THERE) GO TO 80
      IF(LISTCC.GT.1) WRITE(KO,15) NL,D,C,IN,KNP,KN,NUMLT,LTMIN
15    FORMAT('0NL. interaction #',I3,' to',I3,' from',I3,' :',5I3)
      LN1 = LN + 1
      LNP1 = LNP + 1
      DO 60 NU1=1,LN1
         NU = NU1 - 1
         LNMNU = LN-NU
      DO 60 NUP1=1,LNP1
         NUP = NUP1 - 1
         LNPMNU = LNP-NUP
         LAMIN1 = ABS(LNPMNU-NU)+ 1
         LAMAX1 =      LNPMNU+NU + 1
         LAMIN2 = ABS(LNMNU-NUP)+ 1
         LAMAX2 =      LNMNU+NUP + 1
         LMAX = MIN(LC+LAMAX1-1, LD+LAMAX2-1)
         LMAXP = LMAX + 1
         LMINP= MAX(LC-(LAMAX1-1), LD-(LAMAX2-1)) + 1
      R2 =       FACT(2*LN+1 +1)-FACT(2*NU +1)-FACT(2*LNMNU+1 +1)
     &         + FACT(2*LNP+1+1)-FACT(2*NUP+1)-FACT(2*LNPMNU+1+1)
      R2 = EXP( R2 * 0.5D0 )
      R3 = SQRT(ONE*(2*LN+1)*(2*LNMNU+1)*(2*LNP+1)*(2*LNPMNU+1))
      DO 55 L1=LMINP,LMAXP
         L = L1-1
C
      RL = 0.0
      DO 25 LAM11=LAMIN1,LAMAX1
         LAM1 = LAM11-1
      DO 25 LAM21=LAMIN2,LAMAX2
         LAM2 = LAM21-1
      IF(ODD(LAM1+LC+L).OR.ODD(LAM2+LD+L).OR.ODD(LNPMNU+NU+LAM1)
     &   .OR. ODD(LNMNU+NUP+LAM2)) GO TO 25
      RL1 = (TWO*LAM1+ONE) * (TWO*LAM2+ONE)
      RL2 = WIG3J(LAM1+Z,LC+Z,L+Z,Z,Z,Z)
      RL3 = WIG3J(LAM2+Z,LD+Z,L+Z,Z,Z,Z)
      RL4 = WIG3J(LNPMNU+Z,NU+Z,LAM1+Z,Z,Z,Z)
      RL5 = WIG3J(LNMNU+Z,NUP+Z,LAM2+Z,Z,Z,Z)
      RLP =       RL1*RL2*RL3*RL4*RL5
      IF(ABS (RLP).LT.EPS) GO TO 25
      IQMIN = MAX(ABS(LN-LNP),ABS(LAM1-LAM2))  + 1
      IQMAX = MIN(    (LN+LNP),    (LAM1+LAM2))  + 1
      RQ = 0.0
      DO 22 IQ1=IQMIN,IQMAX
         IQ = IQ1 - 1
      RQ1 = (2*IQ+1)*(2*L+1d0)
      RQ2 = RACAH(LAM1+Z,LC+Z,LAM2+Z,LD+Z,L+Z,IQ+Z)
      RQ3 = WIG9J(LN+Z    ,IQ+Z,    LNP+Z,
     &            NU+Z,    LAM1+Z,  LNPMNU+Z,
     &            LNMNU+Z, LAM2+Z,  NUP+Z)
      RQP =       RQ1 * RQ2 * RQ3
      IF(ABS (RQP).LT.EPS) GO TO 22
      RT = 0.0
         DO 20 ILT=1,NUMLT
            LTOTAL =LTMIN + ILT - 1
         RT1 = KCOEF(IN,ILT)
         IF(ABS(RT1).LT.EPS) GO TO 20
         RT2 = RACAH(LNP+Z,LC+Z,LN+Z,LD+Z,LTOTAL+Z,IQ+Z)
         RT3 = (-1) ** (LTOTAL + L + LC - LD)
      IF(LISTCC.GT.2) THEN
      REST = RL1*RL4*RL5 * (TWO*IQ+ONE)*RQ3 * RT2 * SQRT(TWO*LC+ONE)
                   WRITE(KO,18) D,C,KN,NU,NUP,L,LD,LC,LN,LNP,LAM1,LAM2,
     &   IQ,LTOTAL,RT1,R2,RL3,RL2,R3         ,RQ2              ,REST
C    #   IQ,LTOTAL,RT1,R2,RL3,RL2,R3*(2*L+1.),RQ2*SQRT(2*LD+1.),REST
18    FORMAT('0INNER :',14I3,7F9.4/)
      IF(LISTCC.GT.3) WRITE(KO,*)  RL1,RL4,RL5,RQ3 , RT2
         ENDIF
         RT = RT +        RT1 * RT2 * RT3
20       CONTINUE
      RQ = RQ + RQP * RT
22    CONTINUE
      RL = RL + RLP * RQ
25    CONTINUE
      T = RL * R2 * R3 * R7
      IF(ABS(T).LT.EPS) GO TO 55
      PRES = .TRUE.
      IF(C1FR) THEN
      TC = T * A**NU * B**LNMNU * Q**NUP * P**LNPMNU
C    TO GET  (A*RF)**NU      * (B*RT)**(LN-NU)     IF C1FR
C          * (P*RF)**(LNP-NUP)*(Q*RT)**NUP
         NRF = LNPMNU + NU
         NRT = LNMNU + NUP
      ELSE
      TC = T * P**NU * Q**LNMNU * B**NUP * A**LNPMNU
C    TO GET  (B*RT)**NUP     * (A*RF)**(LNP-NUP)     IF NOT C1FR
C          * (Q*RT)**(LN-NU)  *(P*RF)**NU
         NRF = LNPMNU + NU
         NRT = LNMNU + NUP
      ENDIF
C
      IF(LISTCC.GT.1)WRITE(KO,30) D,C,KN,NU,NUP,
     & L,LD,LC,LN,LNP,          RL,R2,R3,R7,    T,TC
30    FORMAT('0' ,10I3,E13.6,3F10.5,   1P,2E20.9)
         IF(MAXLRQ.LT.L) MAXLRQ = L
         IF(MINLRQ.GT.L) MINLRQ = L
      IF(L .LE. MAXL.AND.L .GE. MINL) GO TO 40
      NL = NL - 1
      GO TO 100
35    FORMAT(' OVERLAP MULTIPOLE ORDER',I4,' IS REQUIRED, AND ',A3,'IMUM
     & PRE-CALCULATED IS',I4,' SO SOME CHANNELS ARE UNCOUPLED.')
40    IF(WHERE(L1,IN).NE.0) GO TO 47
         IB = IB + 1
         IF(IB.GT.NIB) IB = 1
         WHERE(L1,IN) = IB
         IF(WHOL(IB).NE.0) then
 		I = WHOL(IB); J=WHOI(IB)
 		WHERE(I,J) = 0
 		endif
         WHOL(IB) = L1
         WHOI(IB) = IN
      IF(LISTCC.GT.3) WRITE(KO,46) L1,IN,IB,NIB
46    FORMAT(' Read multipole #',I3,' of pair',I2,' to',I3,' in block of
     &',I4)
         IF(DRY) GO TO 55
C
      NREC =  (IN-1)*(MAXL1-MINL1+1) + L1 - MINL1    + 1    +  FI

!!!	write(6,*) 'NREC,IB,NLO,NLL = ',NREC,IB,NLO,NLL
      IF(     FFR) READ(11,REC=NREC) ((QRLN(J,I,IB),J=1,NLO),I=2,NLL)
!     IF(     FFR) CALL FIOR(11,QRLN(1,1,IB),NLO*(NLL  ),NREC,ISTAT)
      IF(.NOT.FFR) READ(9 ,REC=NREC) ((QERN(J,I,IB),J=1,NLO),I=2,NLL)
!     IF(.NOT.FFR) CALL FIORC( 9,QERN(1,1,IB),NLO*(NLL  )*2,NREC,ISTAT)
      IF(LISTCC.GT.10.and..not.FFR) then
!!!	write(6,*) 'QERN from NREC = ',NREC,' @ ',IN,MINL1,MAXL1,L1,FI
         	CALL DISPLR(QERN(1,1,IB),NLO,NLL,NLO,SCALR)
		call flush(6)
           	SCALI = SCALR * 1e-12
         	CALL DISPLI(QERN(1,1,IB),NLO,NLL,NLO,SCALI)
	endif
C
 47      RQ  = ONE/(MLT * HF)
         CR  = (CUTOFF-1) * HF
      DO 54 I=2,NLL
      RT = (I-1)*RINTO
         RF = RT * HFHT
      TT = TC * RT**NRT
         JMIN = MAX( NLC+2 + INT((CR-RF)*RQ),1)
      IF(FFR) THEN
         DO 50 J=JMIN,NLO
   50    FNL(I,J) = FNL(I,J) + TT * QRLN(J,I,WHERE(L1,IN))
     &                     * (RF + DNL(J))**NRF
      ELSE
         DO 51 J=JMIN,NLO
   51    FNC(I,J) = FNC(I,J) +      QERN(J,I,WHERE(L1,IN))
     &                     *(TT* (RF + DNL(J))**NRF )
      ENDIF
      IF(NLPL.LE.0.AND.JTEST.GT.4) GO TO 54
        IF(FFR) THEN
         DO 52 J=JMIN,NLO
   52    VINT = MAX(VINT,ABS(FNL(I,J)))
        ELSE
         DO 53 J=JMIN,NLO
   53    VINT = MAX(VINT,ABS(FNC(I,J)))
        ENDIF
54    CONTINUE
55    CONTINUE
60    CONTINUE
      VINT = VINT * 3.0
      ALOSS(IN) = MAX(ALOSS(IN),VINT)
      ALOSSM = MAX(ALOSSM,VINT)
80    CONTINUE
!      WRITE(48,90) D,C,CP
!      WRITE(48,99) ALOSSM,ALOSSM*RER*100.
!      written(48) = .true.
      IF(DRY) GO TO 100
      IF(.NOT.PRES) NL = NL - 1
      IF(.NOT.PRES) GO TO 100
      IF(FFR) THEN
             WRITE(12) FNL
      ELSE
              WRITE(12) FNC
      ENDIF
C
C     IF(.NOT.(NLPL.GT.0.AND.NL.LE.1)) GO TO 100
      IF(NLPL.LE.0) GO TO 100
      WRITE(KO,90) D,C,CP
  90  FORMAT(' NL interaction V-c(',I2,') c(',I2,') of Cplg',I3,' is'/)
      IF(FFR) CALL DISPLY(FNL,NLL,NLO,NLL,SCALE)
      IF(.NOT.FFR) THEN
         CALL DISPLR(FNC,NLL,NLO,NLL,SCALR)
           SCALI = SCALR * RER * 1E3
         CALL DISPLI(FNC,NLL,NLO,NLL,SCALI)
         SCALE = MAX(SCALR,SCALI)
         ENDIF
      IF(SCALE.ne.0.) then
      NLPL = NLPL - 1
       IF(FFR) THEN
         DO 95 I=2,NLL
         FNL(I,1) = 0
         DO 95 J=1,NLO
95       FNL(I,1) = FNL(I,1) + FNL(I,J) * HF *MLT
         WRITE(KO,98) (FNL(I,1),I=2,NLL)
98       FORMAT(1X,18F7.2)
       ENDIF
       endif
      WRITE(KO,99) ALOSSM,ALOSSM*RER*100.
99    FORMAT(1X,'Largest intermediate sum =',1P,E12.4,
     X ', so expect errors of',F10.4,' %')
100   CONTINUE
      IF(NL.GT.MFNL)  WRITE(KO,9) NL,MFNL
      IF(MAXLRQ.GT.MAXL) WRITE(KO,35) MAXLRQ,'MAX',MAXL
      IF(MINLRQ.LT.MINL) WRITE(KO,35) MINLRQ,'MIN',MINL
!     if(FFR) then
!	deallocate (QRLN,FNL)
!	else
!     	deallocate (QERN,FNC)
!	endif
      RETURN
      END
      SUBROUTINE QERNEL(IFL,NLN,NLM,NLO,MAXL1,NK,A,B,P,Q,FPT,NG,
     &             IC7,ICV,IP1,IREM,NNT,MINL1,RINTO,EPC,HNL,NLT,
     &             VFOLD,VCORE,RINS,RINC,NLL,NIBL,cxwf,FORML,VSP,
     &             FORMC,CSP,MINT,RIN,NLC,WID,CENTR,CENTRE,NONO,HFHT,
     &             FFREAL,JACOBN,LTRANS)
	use parameters
	use io
	use drier
	use trace
	use fresco1, only:rnl
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 JACOBN,COFFIN(4),XG(6),WG(6),C2,CI,WD,TH
      REAL*8,allocatable:: QRLN(:,:,:)
      COMPLEX*16,allocatable:: QERN(:,:,:)
      INTEGER FPT(6,NK),KLIM(MAXNLN)
      COMPLEX*16 VCOR(MAXNNU),VSUBI,DV,FFC,WC,SUM(MAXL1),WRT,WRP
      LOGICAL FFREAL,LTRANS(NG),THREV,VCC,cxwf
      REAL*8 RP4,NLOC(NLM),VR,FFR,FFR4,RCOR2,VR1,VR2,
     X       RUM(MAXL1),VSP(MAXNLR,MSP),FORML(MAXNLR,MSP)
      COMPLEX*16 CSP(MAXNLC,MSP),FORMC(MAXNLC,MSP),FFC4
      COMPLEX*16 VFOLD(MAXNLN),VCORE(MAXNLN),DV0,DV1,VSUBF
      REAL*8 PLEG(MAXNNU,MAXMUL),UK(MAXNNU),GK(MAXNNU),THM(MAXNLN),
     &       WT(MAXNNU),RN(MAXNNU),RDINT(MAXNNU),RCOR(MAXNNU)
      DATA NWW,XG(4),XG(5),XG(6),WG(4),WG(5),WG(6)/3,
     1   .2386191861D0,.6612093865D0,.9324695142D0,
     2   .4679139346D0,.3607615730D0,.1713244924D0/, PI /3.14159D0/
      CALL CCTIME(IPTIME)
      if(FFREAL) then
	allocate (QRLN(NLO,NLL,NIBL))
	else
      	allocate (QERN(NLO,NLL,NIBL))
	endif
!	write(6,*) ' IC7,ICV,FFREAL =',IC7,ICV,FFREAL
      THMN = 0.1 *PI/180.  ! at least 0.1 degres, and moreover non-zero!
      EP = EPC * 0.01
      IF(EP.EQ.0.) EP = 1.D-5
      HIRS = 1.0/RINS
      HIRC = 1.0/RINC
      RDRPM = 1.0/((NLL-1)*RINTO)**2
      NW=2*NWW
      DO 1 N=1,NWW
      NN=NW-N+1
      XG(N)=-XG(NN)
  1   WG(N)=WG(NN)
C
           MAXL = MAXL1 - 1
           IC3 = NNT/NW
           NNU = IC3 * NW
            CI = 1D0/IC3
            C1 = 0.5D0 * CI
      DO 10 I=1,NLM
10    NLOC(I) = 0.0
      THREV = IP1.LT.-1
      VCC   = IP1.LE.-3
      PA = P - A
      QB = Q - B
      PQ = P*HFHT + Q
      AB = A*HFHT + B
      PAQB2 = PA*QB*2.
      PAQB  = PA*HFHT+QB
      AB2 = A*B*2
      PQ2 = P*Q*2
      C = 0.5 * JACOBN
      NREC = IFL
      NBLOCK = (MAXL1-MINL1+1 + NIBL-1) / NIBL
      DO 150 IGG=1,NG
      IF(.NOT.LTRANS(IGG)) GO TO 150
        IMINL1 = MINL1
         THMAX = PI
         XTI = 0.0
      DO 120 IBL=1,NBLOCK
        IMAXL1 = MIN(IMINL1 + NIBL - 1, MAXL1)
        IMAXL  = IMAXL1 - 1
        IMINL  = IMINL1 - 1
      DO 100 I=2,NLL
         RP = (I-1) * RINTO
       IF(IBL.EQ.1) THEN
         IF(I.LE.3) GO TO 11
         IF(KLIM(I-1).EQ.NNU) THMAX = MIN(THMAX*1.333,PI)
         IF(KLIM(I-1).LT.NNU.AND..NOT.THREV)
     X                             THMAX = ACOS(UK(KLIM(I-1))+1.)
         IF(KLIM(I-1).LT.NNU.AND.THREV)
     X                        THMAX = PI - ACOS(UK(KLIM(I-1))+1.)
	 X = THMAX
	 THMAX = max(THMAX,THMN)  ! at least non-zero!
 11      KLIM(I) = NNU/2
         THM(I) = THMAX
         XTL = XTI
         XTI = 0.0
       ELSE
        THMAX = THM(I)
       ENDIF
            C2 = C1
               K = 0
            DO 14  J=1,IC3
            DO 13 NN=1,NW
               K = K + 1
               X = C1 * XG(NN) + C2
               TH = (3D0*X*X + 1D0)*X*THMAX*0.25D0
               IF(THREV) TH = PI - TH
               UK(K) =  COS(TH) - 1D0
               GK(K) = C1 * WG(NN)
               WT(K) =  SIN(TH) * (9D0 *X*X+1D0 )*THMAX*0.25D0
               PLEG(K,1) = 1D0
               if(MAXMUL>1) PLEG(K,2) = UK(K) + 1D0
13             CONTINUE
14           C2 = C2 + CI
               DO 15 L=2,IMAXL
               DO 15 K=1,NNU
               WD = UK(K)  + 1D0
15         PLEG(K,L+1) = ((2*L-1)*WD*PLEG(K,L-1+1) -(L-1)*PLEG(K,L-2+1))
     &                       / DBLE(L)
      VSUBF = 0.0
      DO 25 K=1,NNU
25    VCOR(K) = 0.0
         IF(IREM.NE.0.AND.IC7.ne.1) then
	   if(.NOT.FFREAL) VSUBF = FFC(RP*HIRS,VFOLD,NLN)   ! Vopt subtracted in post
	   if(     FFREAL) VSUBF = FFR(RP*HIRS,VFOLD,NLN)
	 endif
         DO 29 JJ=1,NLO
      IF(FFREAL) THEN
         DO 27 L1=IMINL1,IMAXL1
27       QRLN(JJ,I,L1-IMINL) = 0.0
      ELSE
         DO 28 L1=IMINL1,IMAXL1
28       QERN(JJ,I,L1-IMINL) = 0.0
      ENDIF
29    CONTINUE
C
         G = 4./3.
         G = 1.
         DO 85 J=1,NLM
            DNL =     (J - NLM/2 - 1) * HNL  +  CENTRE
            RT = RP*HFHT
            RDCM = RT    + DNL
            PQR = (P * DNL + PQ * RP)**2
            ABR = (A * DNL + AB * RP)**2
            PAQBR = (PA * DNL + (PAQB)*RP)**2
            RDRP = RDCM*RP
            DO 30 L1=IMINL1,IMAXL1
               IF(FFREAL) RUM(L1)      = 0.
30             IF(.NOT.FFREAL) SUM(L1) = 0
            IF(RDCM.LE.0.01) GO TO 60
            DO 31 K=1,NNU
               RDINT2 =     PQR + PQ2*RDRP*UK(K)
               RN2  =       ABR + AB2*RDRP*UK(K)
               RDINT(K) = SQRT(ABS(RDINT2)) * RIN
31             RN(K) =    SQRT(ABS(RN2))    * RIN
               IF(IREM.EQ.0) GO TO 32
            IF(IC7.ne.0.AND..NOT.FFREAL) VSUBI= FFC(RDCM*HIRS,VFOLD,NLN)   ! prior
            IF(IC7.ne.0.AND.     FFREAL) VSUBI= FFR(RDCM*HIRS,VFOLD,NLN)
            IF(IC7.eq.0)                 VSUBI= VSUBF			   ! post
	    if(IC7.eq.2) VSUBI = VSUBI - VSUBF			   	   ! prior - post
               IF(VCC) VSUBI = 0.0
                  DO 315 K=1,NNU
                  RCOR2 = PAQBR + PAQB2*RDRP * UK(K)
 315              RCOR(K) = SQRT(ABS(RCOR2)) * HIRC
32             DO 50 K=1,NNU
		   VCOR(K) = 0.0
                 IF(IREM.ne.0) then
 		   if(IC7.ne.2) VCOR(K)=FFC(RCOR(K),VCORE,NLN)
                   VCOR(K)=VCOR(K) - VSUBI
	         endif
               DO 38 IN=1,NK
                   IG  = FPT(3,IN)
                   IF(IG.NE.IGG) GO TO 38
                   IFT = FPT(2,IN)
                   IFP = FPT(1,IN)
		   if(.not.cxwf) then
                   WRT= FFR4(RN(K),FORML(1,IFT),MINT)
                   WRP= FFR4(RDINT(K),FORML(1,IFP),MINT)
		   else
                   WRT= FFC4(RN(K),FORMC(1,IFT),MINT)
                   WRP= FFC4(RDINT(K),FORMC(1,IFP),MINT)
		   endif
         IF(.NOT.FFREAL) THEN
C                            Complex form factors
                   IF(NONO.GT.0) THEN
                      DV = WRT * WRP
                      ELSE
		  DV0=0.; DV1=0.;
	 if(.not.cxwf) then
            IF(ICV.ne.1) DV0=FFR4(RDINT(K),VSP(1,IFP),MINT) * WRT  ! with post vertex function
            IF(ICV.ne.0) DV1=FFR4(RN(K)   ,VSP(1,IFT),MINT) * WRP  ! with prior
	 else
            IF(ICV.ne.1) DV0=FFC4(RDINT(K),CSP(1,IFP),MINT) * WRT  ! with post vertex function
            IF(ICV.ne.0) DV1=FFC4(RN(K)   ,CSP(1,IFT),MINT) * WRP  ! with prior
	 endif
		if(IC7<2)  DV = DV1 + DV0    ! only 1 in non-zero
		if(IC7==2) DV = DV1 - DV0    ! prior-post
                IF(VCC) DV = 0.0
                IF(IREM.ne.0) DV = DV + VCOR(K) *WRT*WRP
                      ENDIF
               WC= DV * WT(K)
                 IF(IBL.EQ.1) THEN
                   W = WC
                   IF(ABS(W).GT.EP*XTL) KLIM(I) = MAX(KLIM(I),K)
                   XTI = MAX(XTI,ABS(W))
                 ENDIF
                WC = WC * GK(K)
C
               DO 35 L1=IMINL1,IMAXL1
               SUM(L1) = SUM(L1) + WC * PLEG(K,L1)
35             CONTINUE
            ELSE
C                        Real form factors
                   IF(NONO.GT.0) THEN
                      VR = WRT * WRP
                      ELSE
         IF(ICV.ne.1) VR0=FFR4(RDINT(K),VSP(1,IFP),MINT) * WRT
         IF(ICV.ne.0) VR1=FFR4(RN(K)   ,VSP(1,IFT),MINT) * WRP
                if(IC7<2)  VR = VR1 + VR0    ! only 1 in non-zero
                if(IC7==2) VR = VR1 - VR0    ! prior-post
                IF(VCC) VR = 0.0
                IF(IREM.NE.0) VR = VR + DBLE(VCOR(K)) *WRT*WRP
                      ENDIF
               WD= VR * WT(K)
                IF(IBL.EQ.1) THEN
                   W = WD
                   IF(ABS(W).GT.EP*XTL) KLIM(I) = MAX(KLIM(I),K)
                   XTI = MAX(XTI,ABS(W))
                 ENDIF
                WD = WD * GK(K)
               DO 36 L1=IMINL1,IMAXL1
36             RUM(L1) = RUM(L1) + WD * PLEG(K,L1)
            ENDIF
38      CONTINUE
C
50          CONTINUE
              C2 = RDRP * RDRPM
               IF(FFREAL) THEN
                 DO 57 L1=IMINL1,IMAXL1
57                NLOC(J)=NLOC(J) +    (RUM(L1-IMINL))**2 * C2
               ELSE
                 DO 58 L1=IMINL1,IMAXL1
58                NLOC(J)=NLOC(J) + ABS(SUM(L1-IMINL))**2 * C2
               ENDIF
60          CALL SPLINT( DNL/(HNL*NLT) + NLC    ,NLO,IJ,NP,COFFIN)
            G = 2. - G
            DO 82 M=1,NP
            WD  = (COFFIN(M) * G / NLT) * C * RDRP
            JJ = IJ + M - 1
         IF(FFREAL) THEN
            DO 79 L1=IMINL1,IMAXL1
79          QRLN(JJ,I,L1-IMINL) = QRLN(JJ,I,L1-IMINL) + RUM(L1) * WD
         ELSE
            DO 80 L1=IMINL1,IMAXL1
80          QERN(JJ,I,L1-IMINL) = QERN(JJ,I,L1-IMINL) + SUM(L1) * WD
         ENDIF
82       CONTINUE
85            CONTINUE
100      CONTINUE
              DO 90 L1=IMINL1,IMAXL1
              NREC = (IGG-1)*(MAXL1-MINL1+1) + L1 - MINL1  + 1 + IFL
                 IF(DRY) GO TO 90
            IF(FFREAL) THEN
             WRITE(11,REC=NREC) ((QRLN(JJ,I,L1-IMINL),JJ=1,NLO),I=2,NLL)
!           CALL FIOW(11,QRLN(1,1,L1-IMINL),NLO*(NLL  ),NREC,ISTAT)
            ELSE
             WRITE(9 ,REC=NREC) ((QERN(JJ,I,L1-IMINL),JJ=1,NLO),I=2,NLL)
!           CALL FIOWC( 9,QERN(1,1,L1-IMINL),NLO*(NLL  )*2,NREC,ISTAT)
!!!	   write(6,*) 'WRITE 9, NREC,L1,NLO,NLL = ',NREC,L1,NLO,NLL
!!!         	CALL DISPLR(QERN(1,1,L1-IMINL),NLO,NLL,NLO,SCALR)
!!!           	SCALI = SCALR * 1e-12
!!!         	CALL DISPLI(QERN(1,1,L1-IMINL),NLO,NLL,NLO,SCALI)
            ENDIF
90         CONTINUE
120    IMINL1 = IMAXL1 + 1
150    CONTINUE
      IFL = IFL + NG * (MAXL1-MINL1+1)
      CALL CCTIME(IPTIME)
       PT = IPTIME * 0.01
      IF(LISTCC.GE.3) WRITE(KO,992) PT,(THM(I)*180/PI,I=2,NLL)
992   FORMAT(/' Theta - maxima (after',F6.2,' secs) are'/(1X,20F6.1))
      RP4 = 0.0
      DO 210 J=1,NLM
210   RP4 = MAX(NLOC(J),RP4)
      VR = 100.*EP * 0.1
      IF(ABS(RP4).LT.1.E-30) WRITE(KO,992) PT,(THM(I),I=2,NLL)
      IF(ABS(RP4).LT.1.E-30) GO TO 219
      IFT = 10000
      IFP = -10000
      DO 215 J=1,NLM
            DNL = (J - NLC - 1) * HNL
      NLOC(J) = SQRT( NLOC(J) / RP4) * 100.
      IF(NLOC(J).LE.VR) GO TO 215
         IFP = MAX(IFP,J)
         IFT = MIN(IFT,J)
215   CONTINUE
      IF(NLOC(1).GT.VR)
     X IFT = 1 - MAX(INT(LOG(NLOC(1)/VR)/LOG(NLOC(2)/NLOC(1))),0)
      IF(NLOC(NLM).GT.VR)
     X IFP =NLM+MAX(INT(LOG(NLOC(NLM)/VR)/LOG(NLOC(NLM-1)/NLOC(NLM))),0)
      WID  = (IFP - IFT) * HNL
      CENTR = ((IFP + IFT) * 0.5 - NLM/2-1) * HNL + CENTRE
219   WRITE(KO,220) WID,CENTR,PT,(NLOC(J),J=1,NLM)
220   FORMAT('0RECOMMENDED NON-LOCAL WIDTH IS GREATER THAN',F8.2,' FM.',
     &     ',  RECOMMENDED CENTRATION ',F6.2,',',
     &' after',F8.2,' secs.'/'0Relative non-local usages (rms) are'
     &                       /(1X,15F8.3))
      if(WID>rnl*1.5) then
	write(KO,221)  rnl,WID
  221 	format(' Non local width RNL',f8.3,' is too small: INCREASE TO',
     x              ' >',f8.3,', so STOP!')
	  stop
	endif
      if(FFREAL) then
	deallocate (QRLN)
	else
      	deallocate (QERN)
	endif
      RETURN
      END
      SUBROUTINE SOURCE(INHOMG,PSI,N,H,NCH,IEX,FORMF,NF,FORMFR,CUTOFF,
     &  ICUTC,SIMPLE,ITNL,NLN,NLO,MR,NL,EMPTY,SAME,SH,WAVES,LVAL,
     &  MLT,SKIP,FED,NLC, CHNO,MLM,NICH,NAXICH,PTYPE,LOCFIL,
     &  CLIST,NCLIST,NFLIST,GAM,M1GAM,M2GAM,EXCIT,PART) ! Modified by MGR & AMoro
!     &  CLIST,NCLIST,NFLIST)
	use parameters
	use io
       use fresco1, only: rela ! AMoro
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 PSI(MAXN,NICH),INHOMG(MAXN,MAXCH),FORMF(MAXM,MLOC),
     &      FORMFR(NF),
     &      CLIST(MAXCH,MAXCH,MCLIST),S,ECF(MAXNLN),FNC(NLN,NLO),CI,PH
      INTEGER C,D,WAVES,CUTOFF,CHNO(MFNL,6),LVAL(MAXCH),PTYPE(8,NF),
     & 	    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH)
      REAL*8 EPS,EXF(NLN),FNL(NLN,NLO)
      REAL*8 H(NCH)
      LOGICAL ITNL,EMPTY(MAXCH),SAME(MAXCH),SKIP(MAXCH),NREV,NFOR,SH,
     &        SIMPLE(MAXCH),FFR,FED(MAXCH),CP,LOCFIL
!MGR & AMoro ------------------------------------------------------------
      INTEGER M1GAM, M2GAM
      REAL*8 GAM(M1GAM,M2GAM)
      INTEGER PART(NCH,3),EXCIT(NCH,3)
!MGR & AMoro ------------------------------------------------------------     
      DATA EPS / 1E-12 /, CI / (0.0,1.0) /
C
      DO 3 C=1,NCH
        SKIP(C) = SIMPLE(C)
        FED(C) = .FALSE.
        DO 3 D=1,NCH
        DO 3 NC=1,NCLIST(C,D)
	 IF = NFLIST(C,D,NC)
         IF((C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX)
     X      .AND.mod(PTYPE(6,IF),2).EQ.0) GO TO 3
            CP = ABS(CLIST(C,D,NC)).GT.EPS
        SKIP(C) = SKIP(C) .AND. (SAME(D).OR..NOT.CP)
        FED(C) = FED(C) .OR. CP.AND..NOT.EMPTY(D)
3     CONTINUE
      ICH = 0
      DO 4 INL=1,NL
            IF(MOD(CHNO(INL,3),10).LE.1) GO TO 4
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         SKIP(D) = SKIP(D) .AND. SAME(C)
         SKIP(C) = SKIP(C) .AND.(SAME(D) .OR. CHNO(INL,1).LT.0)
         FED(D) = FED(D) .OR. .NOT.EMPTY(C)
         FED(C) = FED(C) .OR. .NOT.EMPTY(D)  .AND.CHNO(INL,1).GT.0
        IF(.NOT.EMPTY(C)) ICH = MAX(C,ICH)
        IF(.NOT.EMPTY(D).AND.CHNO(INL,1).GT.0) ICH = MAX(D,ICH)
4     CONTINUE
      NAXICH = MAX(ICH,NAXICH)
!      CALL CHECK(ICH,MAXICH,-25)
	if(ICH>MAXICH) then
	write(KO,7) ICH,MAXICH+1,ICH
7	format(//'  ***** Channels up to',i4,' have unexpected flux,'/
     X		 '        probably from long-range Coulomb couplings.'/
     X		 '        FLUX IN CHANNELS ',i4,' TO',i4,' IS IGNORED'/)
	SH=.true.
	endif
      IF(SH)WRITE(KO,8) (C,EMPTY(C),SAME(C),SIMPLE(C),
     X                 SKIP(C),FED(C),C=1,NCH)
8     FORMAT(' #,EMPTY,SAME,SIMPLE,SKIP,FED =',12(1X,I2,5L1,';') )
C
      DO 10 C=1,NCH
      IF(SKIP(C)) GO TO 10
      DO 9 I=1,N
9     INHOMG(I,C) = 0.0
10    CONTINUE
C
      if(.not.LOCFIL) then
      DO 15 D=1,NCH
      DO 15 C=1,NCH
      DO 15 NC=1,NCLIST(C,D)
	 JF = NFLIST(C,D,NC)
C        IF(C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX) GO TO 15
         IF((C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX)
     X      .AND.mod(PTYPE(6,JF),2).EQ.0) GO TO 15
         IF(D.gt.NICH) go to 15
         S = CLIST(C,D,NC)
         IF(ABS(S).LT.EPS.OR.SKIP(C).OR.EMPTY(D)) GO TO 15
C     IF(SH) WRITE(KO,11) C,D,JF,S,IEX
      IF(ABS(WAVES).GE.2) WRITE(41,11) C,D,JF,S,IEX
11    FORMAT('#ZR to',I3,' fr',I3,' by',I3,' of',2E14.3,'(IEX=',I3,')')
         DO 12 I=max(1,ICUTC),N
      INHOMG(I,C) = INHOMG(I,C) + S * FORMF(I,JF) * PSI(I,D)
      IF(ABS(WAVES).GE.2.AND.ABS(WAVES).le.3)
     X   WRITE(41,13) (I-1)*H(D) , INHOMG(I,C),PSI(I,D),FORMF(I,JF)
12    CONTINUE
13     FORMAT(F8.4,1P,6E12.3)
      IF(ABS(WAVES).GE.2.AND.ABS(WAVES).le.3) write(41,*) '&'
15    CONTINUE
	else
      rewind 19
      IS=max(1,ICUTC)
	 do I=1,IS-1
	 read(19)
	 enddo
      DO 26 I=IS,N
	 read(19) FORMFR
      DO 25 D=1,NCH
      DO 25 C=1,NCH
      DO 25 NC=1,NCLIST(C,D)
	 JF = NFLIST(C,D,NC)
         IF((C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX)
     X      .AND.mod(PTYPE(6,JF),2).EQ.0) GO TO 25
         IF(D.gt.NICH) go to 25
         S = CLIST(C,D,NC)
         IF(ABS(S).LT.EPS.OR.SKIP(C).OR.EMPTY(D)) GO TO 25
         IF(ABS(WAVES).GE.2.and.I==IS) WRITE(41,11) C,D,JF,S,IEX
      INHOMG(I,C) = INHOMG(I,C) + S * FORMFR(JF) * PSI(I,D)
25    CONTINUE
26    CONTINUE
	endif
C
      IF(.NOT.ITNL.OR.NL.EQ.0) GO TO 51
      HI = 1.0 / DBLE(MR)
       REWIND 12
      DO 50 INL=1,NL
         FFR = CHNO(INL,3).LT.10
       NLL = CHNO(INL,4)
          IF(FFR) READ(12) ((FNL(I,J),I=1,NLL),J=1,NLO)
          IF(.NOT.FFR) READ(12) ((FNC(I,J),I=1,NLL),J=1,NLO)
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         PH = CI ** (LVAL(C) - LVAL(D))
            IF(MOD(CHNO(INL,3),10).LE.1) GO TO 50
         NFOR = SKIP(D) .OR. EMPTY(C)
         NREV = SKIP(C) .OR. EMPTY(D) .OR. CHNO(INL,1).LT.0 .OR. C.EQ.D
         IF(SH) WRITE(KO,28) INL,D,C,NICH,.NOT.NFOR,.NOT.NREV,FFR,NLL
28       FORMAT(' NL coupling #',I3,' to',I3,' from Ch.',I3,'<=',i3,
     &         ', Forw =',L2,', Reverse=',L2,', Real =',L2,' NLL=',i4)
         IF(NFOR.AND.NREV) GO TO 50
            SQH = SQRT(H(C)/H(D))
C        IF(SH.AND.FFR) CALL DISPLY(FNL,NLL,NLO,NLN,SCALE)
         IF(NFOR .OR. C.GT.NICH) GO TO 32
         DO 31 J=1,MLM
            JJ = (J - NLC*MLT  - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF,ICUTC)
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1)/MR+2,2)
            KMAX = MIN(IMAX/MR,NLL-2)
            DO 31 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(C) * SQH
              Y = Q * P2 * 0.5  * H(C) * SQH
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 29 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
!            INHOMG(I,D) = INHOMG(I,D) + (V*PH) * PSI(I+JJ,C)
!MGR & AMoro -------------------------------------------------------------
!             write(*,*) 'Inserting gamma line 834'
             if ((rela .eq. 'na').or.(rela .eq. '3d') .or.(rela.eq.'c')
     & ) then
!              write(*,*) 'rela ',rela
             INHOMG(I,D) = INHOMG(I,D) + (V*PH) * PSI(I+JJ,C)
     &        *GAM(PART(D,1),EXCIT(D,1)) !MGR include gammas
             else
!             write(*,*) 'rela ',rela
            INHOMG(I,D) = INHOMG(I,D) + (V*PH) * PSI(I+JJ,C)
             endif
!-----------------------------------------------------------------------    
29           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 30 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
!            INHOMG(I,D) = INHOMG(I,D) + (S*PH) * PSI(I+JJ,C)
!MGR & AMoro --------------------------------------------------------------
             

             if ((rela .eq. 'na').or.(rela .eq. '3d') .or.(rela.eq.'og')
     &       ) then

             INHOMG(I,D) = INHOMG(I,D) + (S*PH) * PSI(I+JJ,C)
     &        *GAM(PART(D,1),EXCIT(D,1)) !MGR include gammas
             else
             INHOMG(I,D) = INHOMG(I,D) + (S*PH) * PSI(I+JJ,C)
             endif
!-----------------------------------------------------------------------             
30           I = I + MR
         ENDIF
31      CONTINUE
32       IF(NREV .OR. D.GT.NICH) GO TO 50
         DO 42 J=1,MLM
            JJ =-(J - NLC*MLT - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF,ICUTC)
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1+JJ)/MR+2,2)
            KMAX = MIN((IMAX+JJ)/MR,NLL-2)
            DO 42 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(D) / SQH
              Y = Q * P2 * 0.5  * H(D) / SQH
          IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 38 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
!            INHOMG(I,C) = INHOMG(I,C) + (V*CONJG(PH)) * PSI(I+JJ,D)
!MGR & AMoro -------------------------------------------------------------
             if ((rela .eq. 'na').or.(rela .eq. '3d') .or.(rela.eq.'c')
     &        ) then
             INHOMG(I,C) = INHOMG(I,C) + (V*CONJG(PH)) * PSI(I+JJ,D)
     &        *GAM(PART(C,1),EXCIT(C,1)) !MGR include gammas
             else
             INHOMG(I,C) = INHOMG(I,C) + (V*CONJG(PH)) * PSI(I+JJ,D)
             endif
!-----------------------------------------------------------------------             
38           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 40 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
!            INHOMG(I,C) = INHOMG(I,C) + S * CONJG(PH) * PSI(I+JJ,D)
!MGR & AMoro ---------------------------------------------------------------
             if ((rela .eq. 'na').or.(rela .eq. '3d') .or.(rela.eq.'c')
     &       ) then
!             write(*,*) 'Inserting gamma line 886'             
             INHOMG(I,C) = INHOMG(I,C) + S * CONJG(PH) * PSI(I+JJ,D)
     &       *GAM(PART(C,1),EXCIT(C,1)) !MGR include gammas
             else
             INHOMG(I,C) = INHOMG(I,C) + S * CONJG(PH) * PSI(I+JJ,D)
             endif
!-----------------------------------------------------------------------             
40           I = I + MR
         ENDIF
42    CONTINUE
50    CONTINUE
51    DO 80 C=1,NCH
80    IF(ABS(WAVES).GE.2.AND.ABS(WAVES).le.3.AND.
     X    .NOT.SKIP(C).AND.FED(C)) WRITE(KO,90) C, (INHOMG(I,C),I=1,N-1)
90    FORMAT(//' For channel no.',I3,', The Inhomogeneous driving',
     & ' Terms are'// 1P,(5(E12.3,' +i*',E9.2,',') ) )
      RETURN
      END
      SUBROUTINE RENO(W,PHI,N,H,NCH,CUTOFF,NLN,NLO,MR,MLM,IC7, NF,
     &      LB,NL,EMPTY,MLT,NLC,SH,CLIST,NCLIST,NFLIST,CHNO,NPOST,
     &   NPRIOR,INHOMG,COEF,SIMPLE,FORMF,FORMFR,LOCFIL,ECMR,LVAL)
	use parameters
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 W(MAXN,MAXCH),PHI(N,2),INHOMG(N),S,ENS,FNC(NLN,NLO),
     &          ECF(NLN),FORMF(MAXM,NF),CI,PH,CLIST(MAXCH,MAXCH,MCLIST),
     & 		FORMFR(NF)
      REAL*8 EPS,FNL(NLN,NLO),EXF(NLN)
      REAL*8 COEF(MAXCH),H(MAXCH),ECMR(MAXCH)
      LOGICAL EMPTY(MAXCH),NFOR,NREV,SIMPLE(MAXCH),NPRIOR,FFR,NPOST,
     & 		LOCFIL
      INTEGER CHNO(MFNL,6),CUTOFF,D,C,SH,LVAL(MAXCH),
     & 	    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH)
C

      EPS = 1.E-20
      CI = (0.0,1.0)
      HI = 1.0 / DBLE(MR)
      LAST1 = 0
      LAST2 = 0
      DO 10 C=1,NCH
      SIMPLE(C) = .TRUE.
10    EMPTY(C) = .FALSE.
      IF(NL.EQ.0) GO TO 51
      REWIND 12
      DO 50 INL=1,NL
         FFR = CHNO(INL,3).LT.10
       NLL = CHNO(INL,4)
          IF(FFR) READ(12) ((FNL(I,J),I=1,NLL),J=1,NLO)
          IF(.NOT.FFR) READ(12) ((FNC(I,J),I=1,NLL),J=1,NLO)
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         PH = CI ** (LVAL(C) - LVAL(D))
            IF(MOD(CHNO(INL,3),10).EQ.2) GO TO 50
         NFOR = EMPTY(C) .OR. CHNO(INL,3).EQ.1-IC7
         NREV = EMPTY(D) .OR. CHNO(INL,3).EQ.IC7 .OR. CHNO(INL,1) .LT. 0
         IF(SH.GE.7) WRITE(KO,22) INL,CHNO(INL,1),C,.NOT.NFOR,
     X                                          .NOT.NREV,CHNO(INL,3)
22       FORMAT(' NONO NL coupling #',I3,' to',I3,' from Ch.',I3,', Forw
     & =',      L2,', Reverse=',L2,' of CHNO =',I3)
         IF(NFOR.AND.NREV) GO TO 50
            SQH = SQRT(H(C)/H(D))
         IF(SH.GE.8.AND.FFR) CALL DISPLY(FNL(1,1),NLL,NLO,NLN,SCALE)
         IF(NFOR) GO TO 32
            SIMPLE(D) = .FALSE.
                  IF(LAST1.NE.LB+C) READ(8,REC=LB+C) (PHI(I,1),I=1,N)
                  LAST1 = LB+C
             JJ = N/10
           DO 23 I=1,N,JJ
23         IF(ABS(DBLE(PHI(I,1)))+ABS(AIMAG(PHI(I,1))).GT.EPS) GO TO 25
           EMPTY(C) = .TRUE.
            GO TO 32
25       DO 31 J=1,MLM
            JJ = (J - NLC*MLT  - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  CUTOFF
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1)/MR+2,2)
            KMAX = MIN(IMAX/MR,NLL-2)
            DO 31 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(C) * SQH
              Y = Q * P2 * 0.5  * H(C) * SQH
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 29 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            W(I,D) = W(I,D) - V * PH * PHI(I+JJ,1)
29           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 30 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            W(I,D) = W(I,D) - S * PH * PHI(I+JJ,1)
30           I = I + MR
         ENDIF
31    CONTINUE
32       IF(NREV) GO TO 50
            SIMPLE(C) = .FALSE.
            IF(LAST2.NE.LB+D) READ(8,REC=LB+D) (PHI(I,2),I=1,N)
            LAST2 = LB+D
             JJ = N/10
           DO 33 I=1,N,JJ
33         IF(ABS(DBLE(PHI(I,2)))+ABS(AIMAG(PHI(I,2))).GT.EPS) GO TO 35
           EMPTY(D) = .TRUE.
            GO TO 50
35       DO 42 J=1,MLM
            JJ =-(J - NLC*MLT - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  CUTOFF
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1+JJ)/MR+2,2)
            KMAX = MIN((IMAX+JJ)/MR,NLL-2)
            DO 42 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(D) / SQH
              Y = Q * P2 * 0.5  * H(D) / SQH
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 39 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            W(I,C) = W(I,C) - V * CONJG(PH) * PHI(I+JJ,2)
39           I = I + MR
         ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 40 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            W(I,C) = W(I,C) - CONJG(S * PH)* PHI(I+JJ,2)
40           I = I + MR
         ENDIF
42    CONTINUE
50    CONTINUE
51    DO 53 C=1,NCH
53    IF(IC7.EQ.1.AND.NPOST.AND..NOT.SIMPLE(C))
     #    WRITE(8,REC=NCH+C) (W(I,C),I=1,N)
      IF(NPRIOR.AND.IC7.EQ.0) then
      IF(SH.GE.6) WRITE(KO,55) IC7,(SIMPLE(C),C=1,NCH)
55    FORMAT('0RENO',I3,' : sources still simple are ',   40L2)
C56    FORMAT('0RENO',I3,' : empty channels were      ',   40L2)
      DO 70 C=1,NCH
         IF(SIMPLE(C)) GO TO 70
C		Mistake in original theory:
C	   	Should be previous PHI, not current W
         READ(8,REC=C) (PHI(I,1),I=1,N)
      INHOMG(:) = 0.0
      IMIN = 1+ CUTOFF + 3 + MR
      IMAX = N - 5 - 3        - MR
      T = COEF(C) /(H(C)**2  * 12.0)
       if(LOCFIL) then
	 rewind 19
	 do i=1,IMIN-1
	 read(19) 
	 enddo
	 endif
      DO 60 I=IMIN,IMAX
         ENS = - LVAL(C)*(LVAL(C)+1)*COEF(C)/(H(C)*(I-1))**2
     X         - ECMR(C)
	if(LOCFIL) then
	 read(19) FORMFR
         DO 58 NC=1,NCLIST(C,C)
	   JF = NFLIST(C,C,NC)
58       ENS = ENS + CLIST(C,C,NC) * FORMFR(JF)
	else
         DO 59 NC=1,NCLIST(C,C)
	   JF = NFLIST(C,C,NC)
59       ENS = ENS + CLIST(C,C,NC) * FORMF(I,JF)
	endif
60    INHOMG(I) = - (ENS * PHI(I,1) +
     & T * (-PHI(I+2,1)+16.*PHI(I+1,1)-30.*PHI(I,1)+
     & 	                16.*PHI(I-1,1)-PHI(I-2,1)))
      WRITE(8,REC=NCH+C) (INHOMG(I),I=1,N)
70    CONTINUE
	endif
      RETURN
      END
      SUBROUTINE CHECK(I,LIM,KIND)
	use io
	use drier
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*19  LI(30)
      CHARACTER*6 WHO(30)
      DATA LI /
     &'EXCITATION PAIRS   ','MASS PARTITIONS    ','BOUND STATES       ',
     &'COUPLING TYPES     ','CHANNELS           ','NON-LOCAL FORMS.   ',
     &'C.M. RADIAL POINTS ','CM INTERPOLATN PTS.','NL INTERPOLATN PTS ',
     &'COUPLED CHANNELS   ','L-VALS & MULTIPOLES','ANGULAR QUADRTR PTS',
     &'KERNELS/COUPLING   ','MATRIX ROWS:ERWIN  ','BLOCKS OF EQUATIONS',
     &'COUPLED PAIRS OF WF','                   ','P(L,M) VALUES.     ',
     &'                   ','PLACES IN FAM()    ','TOTAL NO. STATES   ',
     &'PROJECTILE M VALUES','IN YSIG ARRAY,     ','LOCAL FORMS        ',
     &'INTERMEDIATE CHANS.','MRXY in EXTERN2    ','MM in EXTERN2      ',
     &'arrays in EXTERN1  ','                   ','COUPLINGS/PW-PAIR  '/
      DATA WHO /'MXX   ','MXP   ','MSP   ','MAXCPL','MAXCH ','MFNL  ',
     &          'MAXN  ','MAXNLN','MAXNLO','      ','LMAX1 ','MAXNNU',
     &          'MAXQRN','MAXNR ','MAXB  ','MPAIR ','NFDEC ','PL-DEC',
     &          '      ','MAXF  ','MXPEX ','MSPIN ','MXYSIG','MLOC  ',
     &          'MAXICH','MRXY  ','MM    ','there ','      ','MCLIST'/
      IF(I.LE.LIM) RETURN
      WRITE(KO,833) LIM,LI(ABS(KIND)),I,WHO(ABS(KIND))
833   FORMAT(' ****** THERE IS ONLY ROOM FOR',I8,1X,A19,' BUT',I8,
     & ' ARE REQUIRED,  SO INCREASE PARAMETER ',A6,' !!')
      DRY = .TRUE.
	write(KO,20)
20	format(/'  *** INTERNAL ERROR IN ALGORITHM TO CALCULATE',
     X ' PARAMETERS!! ***',/,'  *** Please report this (along with',
     X ' a copy of the input file)',/,
     X '      to Ian@kernz.org ***')
      IF(KIND.LT.0) RETURN
C     STOP
      CALL ABEND(8)
      END
      SUBROUTINE CCTIME(I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /TIMER/ START
      T = SECOND()*100
C
C
      I = T - START
      START = T
      RETURN
      END
      SUBROUTINE DISPLR(A,M,N,MA,DMAX)
      use io
      REAL*8 EMAX
      COMPLEX*16 A(MA,N)
      REAL*8 DMAX
      CHARACTER CHARS(11),SIGNS(3),LINE(132),PLUS,CAPI
      LOGICAL P
      DATA  CHARS     / ' ','1',' ','3',' ','5',' ','7',' ','9','A' /
      DATA  SIGNS    / '.','0',' ' /
      DATA PLUS,CAPI / '+','I' /
      MP= MIN(M,130)
      EMAX = 0
      DO 10 I=1,N
      DO 10 J=1,MP
10    EMAX = MAX(EMAX, ABS(DBLE(A(J,I))) )
      DMAX = EMAX
      IF(EMAX.EQ.0) GOTO 100
      SCALE = 10/EMAX
C20   FORMAT(// 20('*') //)
      IF(ABS(LOG10(EMAX)).LE.4) THEN
      WRITE(KO,12) EMAX
       ELSE
      WRITE(KO,13) EMAX
       ENDIF
12    FORMAT(' Real part: full scale =',F12.5/ )
13    FORMAT(' Real part: full scale =',1P,E12.4/ )
      DO 50 I=1,N
      DO 30 J=1,MP
         ANO= SCALE * ABS(DBLE(A(J,I)))
         NO = ANO
         LINE(J) = CHARS(NO+1)
         IF(NO.LT.1 .AND. ANO .GE. 0.10) LINE(J) = PLUS
30    CONTINUE
      WRITE(KO,35) (LINE(K),K=1,MP),CAPI
35    FORMAT(' I', 131A1)
      P = .FALSE.
      DO 40 J=1,MP
         IS = SIGN(1D0,DBLE(A(J,I)))
         P = P .OR. IS.LT.0
         P = P .OR. IS.EQ.0
40       LINE(J) = SIGNS(IS+2)
      IF(P) WRITE(KO,45)   (LINE(J),J=1,MP)
45    FORMAT('+ ',130A1)
50    CONTINUE
C     WRITE(KO,20)
      RETURN
100   WRITE(KO,105)
105   FORMAT(/' Real part everywhere zero???'/)
      RETURN
      END
      SUBROUTINE DISPLI(A,M,N,MA,DMAX)
      use io
      REAL*8 EMAX
      COMPLEX*16 A(MA,N)
      REAL*8 DMAX
      CHARACTER CHARS(11),SIGNS(3),LINE(132),PLUS,CAPI
      LOGICAL P
      DATA  CHARS     / ' ','1',' ','3',' ','5',' ','7',' ','9','A' /
      DATA  SIGNS    / '.','0',' ' /
      DATA PLUS,CAPI / '+','I' /
      MP= MIN(M,130)
      EMAX = 0
      DO 10 I=1,N
      DO 10 J=1,MP
10    EMAX = MAX(EMAX, ABS(AIMAG(A(J,I))) )
C      DMAX (ON INPUT) IS A SMALL FACTOR OF SCALE OF REAL PART
      IF(EMAX.LE.DMAX) GOTO 100
      DMAX = EMAX
      SCALE = 10/EMAX
C20   FORMAT(// 20('*') //)
      IF(ABS(LOG10(EMAX)).LE.4) THEN
      WRITE(KO,12) EMAX
       ELSE
      WRITE(KO,13) EMAX
       ENDIF
12    FORMAT(' Imaginary part: full scale =',F12.5/ )
13    FORMAT(' Imaginary part: full scale =',1P,E12.4/ )
      DO 50 I=1,N
      DO 30 J=1,MP
         ANO= SCALE * ABS(AIMAG(A(J,I)))
         NO = ANO
         LINE(J) = CHARS(NO+1)
         IF(NO.LT.1 .AND. ANO .GE. 0.10) LINE(J) = PLUS
30    CONTINUE
      WRITE(KO,35) (LINE(K),K=1,MP),CAPI
35    FORMAT(' I', 131A1)
      P = .FALSE.
      DO 40 J=1,MP
         IS = SIGN(1D0,DBLE(AIMAG(A(J,I))))
         P = P .OR. IS.LT.0
         P = P .OR. IS.EQ.0
40       LINE(J) = SIGNS(IS+2)
      IF(P) WRITE(KO,45)   (LINE(J),J=1,MP)
45    FORMAT('+ ',130A1)
50    CONTINUE
C     WRITE(KO,20)
      RETURN
100   WRITE(KO,13) EMAX
      RETURN
      END
*****FRXX3B***************************************************************
      SUBROUTINE ERWINCC(ECM,COEF,FORMF,NF,FORMFR,CH,REPEAT,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,N,H,M,MD,SCL,RENORM,
     $  AL,RE,IC,SHOW,SHSMAT,SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     $  CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
	use io
	use factorials
	use parameters
	use drier
	use searchpar, only: final
      use fresco1,   only: rmorto,hort
      IMPLICIT REAL*8(A-H,O-Z)
C
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS BY EXACT CC
C        SOLVE  
C             (COEF(K). D2/DR2 + EN(R,K)).W(R,K)
C                  + SUM(J): COUPL(R,K,J).W(R,J) + INHOMG(R,K) = 0
C           WHERE COUPL(R,K,J) = SUM(NC): FORMF(R,JF)*CLIST(K,J,NC)
C           WHERE    JF = NFLIST(K,J,NC)
C       &   WHERE EN(R,K) IS THE DIAGONAL PART OF COUPL
C
C     ASSUMED FOR ARRAYS HERE THAT N<=MAXN AND NEQS <= MAXCH
C
C  Using 'Enhanced Numerov' of Thorlacius & Cooper (JCP 72(1987) 70)
C   with 5 terms in cosh(sqrt(12T)) expansion, but only diagonal potl.
C
      COMPLEX*16 FIMD(MAXB,MAXCH),FIM(MAXB,MAXCH)
      COMPLEX*16 MAT(2*NEQS,2*NEQS+1),S,ZIV,ZPV,FI(MAXB,MAXCH),
     &           V(MAXB,MAXCH),ZI(MAXB,MAXCH),ZM(MAXB,MAXCH),CZ,SE,
     &           CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI,ZDL,ZHPLUS,ZHMINUS
      COMPLEX*16 FORMF(MAXM,NF),FORMFR(NF),F8,
     &       C,CLIST(MAXCH,MAXCH,MCLIST),COUPL(MAXB,MAXB)
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),CUTOFF,SHOW,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),PTYPE(8,NF)
      LOGICAL SING,SHSMAT,SMATEL,FCWFN,FJSWTCH,REPEAT,LOCFIL
      AMD1(C) = ABS(DBLE(C)) + ABS(AIMAG(C))
      complex*16 wdiag(neqs)  ! AMM

C
      NM1 = N-1
      NR = 2*NEQS
      NP = NR + 1
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
      ZI2 = (0.0D0,0.5D0)
      CI = (0.0D0,1.0d0)
      CZ = (0D0,0D0)
      ONEC = (1D0,0D0)
      RAD = 45d0/atan(1d0)
C     AL = 0.0
      ALUN = 0.0
      IS =  CUTOFF
      IS = MIN(1*N/2, MAX(2, IS ))
       if(SHOW.ge.3) write(6,*) 'ERWINCC: CUTOFF,IS =',CUTOFF,IS
       if(SHOW.ge.3) write(6,*) 'ERWINCC: CUTVAL =',(CUTVAL(I),I=1,NEQS)
      NN = NM1 - IS + 1
      TMAX = 20.
      TMIN = -125.
        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS
C
      IF (FJSWTCH)  THEN
c                              skip nmrov , match CRCWFN to zero
c
761    DO 765 I=1,NR
       DO 765 J=1,NP
765       MAT(I,J) = 0
c
         DO 780 IT=1,NEQS
               MAT(IT+NEQS,IT) = ONEC
               MAT(IT,IT)   = ONEC*CRCRAT(it)
780       CONTINUE
C
      ELSE
C
      IF(SHOW.GE.2) WRITE(KO, 5) NF,M,NEQS,NR,REPEAT
      IF(SHOW.GE.2) WRITE(KO, 4) (H(K),K=1,NEQS)
4     FORMAT(' ERWINCC step sizes are',10F8.5)
5     FORMAT(' ERWINCC given',4I4,L4,', requires',3I6,';',I8,'/',I8,L3)
      CALL CHECK(NEQS,MAXB,15)
	if(REPEAT) go to 61
C
C
          DO 7 K=1,NEQS
            SMAT(K) = CZ
          H2C(K) = H(K)**2 / COEF(K)
7        CONTINUE
C
         DO 705 I=1,N
705           MAGN(I) = 1.0
C
      DO 11 IT=1,NEQS
          DO 10 K=1,NEQS
         V(IT,K)   = 0
         FI(IT,K) = 0
         ZM(IT,K)  = 0
         ZI(IT,K)  = 0
            FIM(IT,K) = 0
            FIMD(IT,K)= 0
10       CONTINUE
11     CONTINUE

C
        DO 141 K=1,NEQS
141      IF(SHOW.GE.6) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
142      FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------
	  if(LOCFIL) then
	  rewind 19
	  do I=1,IS-1
	   read(19) 
	  enddo
	  endif
!        write(*,*)'IS=',IS
!        write(*,*)'CUTVAL=',CUTVAL(1:NEQS)
         DO 60 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	  if(LOCFIL) read(19) FORMFR
          DO 144 K=1,NEQS
             C = CZ
	     if(LOCFIL) then
             DO NC=1,NCLIST(K,K)
	 	  JF = NFLIST(K,K,NC)
             IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMFR(JF)
	     ENDDO
	     else
             DO NC=1,NCLIST(K,K)
	 	  JF = NFLIST(K,K,NC)
             IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMF(I,JF)
	     ENDDO
	     endif
           SMAT(K) = -LL1(K)*RI2 + (-ECM(K) + C) * H2C(K)
144        CONTINUE
C
          DO 13 IT=1,NEQS
            K = IT
              IF(I.ne.CUTVAL(K)) GO TO 13
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT,K) = SCL * EXP(T)
       	IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT,K)
13         CONTINUE
1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)



c AMM....................................................................
       WDIAG(:)=0.
       DO K=1,NEQS
       wdiag(k)=(ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))
       ENDDO
c ......................................................................

         DO 158 IT=1,NEQS
         DO 158 K=1,NEQS
158       FI(IT,K) = ZI(IT,K) *
     X      (ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))

         IF(SHOW.GE.7) THEN
          DO 161 K=1,NEQS
161          WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    SMAT(K)/H2C(K)+ECM(K),FI(1,K)
165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,20E10.1)
          ENDIF
C
         DO 20 K=1,NEQS
            DO 17 J=1,NEQS
17          COUPL(K,J) = 0.0
           IF(I.LT.IC) GO TO 20
           DO 19 J=1,NEQS
            C = CZ
	    if(K==J) go to 19
	   if(LOCFIL) then
            DO 18 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
18          IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMFR(JF)
	   else
            DO 181 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
181          IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMF(I,JF)
	   endif
19         COUPL(K,J) = C * H2C(K)
20	 CONTINUE
C THE DIAGONAL PART COUPL(K,K) IS ZERO
C
!	 CALL POTWF(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	 CALL POTWF2(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)   ! works better even on SUNs!

         DO 46 IT=1,NEQS
         DO 43 K=1,NEQS
            ZIV = ZI(IT,K)
            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
            ZM(IT,K) = ZIV
            ZI(IT,K) = ZPV
43            CONTINUE
46    CONTINUE

c AMoro ........................................................
      norto=nint(hort/h(1))  ! AMM
      r=h(1)*(i-1.)
      if ((r.lt.rmorto).and.(i.gt.is).and.(hort.gt.0)
     & .and.(mod(i-is,norto).eq.0).and.(i.lt.md)) then 
        write(*,'(5x,"-> re-orthogonalizing at i,r=",i4,1f8.1)')i,r
      CALL LQFACT(ZI,ZM,FI,WDIAG,COUPL,R12,NEQS,MAXB,MAXCH)
      endif
c ...............................................................

      DO 468 IT=1,NEQS
       IF(I.EQ.M) THEN
        DO 463 K=1,NEQS
463     FIM(IT,K) = FI(IT,K)
       ELSE IF(I.EQ.MD) THEN
        DO 465 K=1,NEQS
465     FIMD(IT,K) = FI(IT,K)
       ENDIF
468   CONTINUE
            T = 0.0
         DO 47 IT=1,NEQS
            DO 47 K=1,NEQS
            F8 = FI(IT,K)
            T = MAX(T, AMD1(F8))
   47       IF(T .GT. BIG ) GO TO 48
            MAGN(I) = T
         IF(I.NE.NM1) GO TO 60
 48      IF(T.EQ.0.0) GO TO 60
            MAGN(I) = T
         T = RENORM/T
         IF(SHOW.GE.1) WRITE(KO,49) T,I,RENORM
49       FORMAT(' Renormalising by',E12.4,' at',I4,' to max. WF.',E12.4)
         if(t.eq.0d0) then
          write(KO,*) 'RENORMALISING TO 0! at I=',I,' as lWF=',MAGN(I)
          write(KO,*) ' FI:',FI
          T = SMALL
          endif
         FMIN = small/T
         DO 55 IT=1,NEQS
         DO 51 K=1,NEQS
            ZI(IT,K) = ZI(IT,K) * T
            ZM(IT,K) = ZM(IT,K) * T
            FIM(IT,K) = FIM(IT,K) * T
            FIMD(IT,K) = FIMD(IT,K) * T
51       CONTINUE
55       CONTINUE
C
60       CONTINUE
C-----------------------------------------------------------------------
c
61    DO 65 I=1,NR
      DO 65 J=1,NP
65     MAT(I,J) = 0d0
         DO 70 IT=1,NEQS
         DO 70 K=1,NEQS
            MAT(K+NEQS,IT) = FIMD(IT,K)
            MAT(K,IT)   = FIM(IT,K)
70          CONTINUE
C
	ENDIF
C                    MATCH TO REQUIRED ASYMPTOTIC CONDITIONS
C
      IF (FCWFN) THEN
C                           match Numerov to CRCWFN
       DO 82 K=1,NEQS
       DO 82 K2=1,NEQS
           IF(SHOW.ge.4.and.ABS(CGMAT(K2,K,2)).gt.0.0) 
     *       WRITE(6,716) K,K2,CGMAT(K2,K,2),CFMAT(K2,K,2),
     *           CGMAT(K2,K,1),CFMAT(K2,K,1)
716            FORMAT(1X,2I3,': CRCWFN= ',4(D15.7,3X))
c
c james
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,2)+CI*CFMAT(K2,K,2)
         MAT(K2,K+NEQS)=   ZI2*ZDL*ZHPLUS
         ZHPLUS = CGMAT(K2,K,1)+CI*CFMAT(K2,K,1)
         MAT(K2+NEQS,K+NEQS)= ZI2*ZDL*ZHPLUS
c
c        MAT(K2,K+NEQS)=   ZI2*(CGMAT(K2,K,2)+CI*CFMAT(K2,K,2))
c        MAT(K2+NEQS,K+NEQS)= ZI2*(CGMAT(K2,K,1)+CI*CFMAT(K2,K,1))
82       CONTINUE
         DO 85 K=1,NEQS
c james
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,2)-CI*CFMAT(K,EL,2)
         MAT(K,NP)=   ZI2*ZDL*ZHMINUS
         ZHMINUS = CGMAT(K,EL,1)-CI*CFMAT(K,EL,1)
         MAT(K+NEQS,NP)= ZI2*ZDL*ZHMINUS
c
c          MAT(K,NP)=   ZI2*(CGMAT(K,EL,2)-CI*CFMAT(K,EL,2))
c          MAT(K+NEQS,NP)= ZI2*(CGMAT(K,EL,1)-CI*CFMAT(K,EL,1))
85       CONTINUE
      ELSE
c                         match Numerov to Uncoupled Coulomb wfns
       DO 90  K=1,NEQS
         MAT(K,K+NEQS) = CH(1,K)
         MAT(K+NEQS ,K+NEQS) = CH(2,K)
C                                      C(K) == MAT(K,NP) IS THE RHS.
         MAT(K,NP) = 0d0
         MAT(K+NEQS ,NP) = 0d0
90       CONTINUE
c                                         note minus as ch = i/2*H+
         MAT(EL,NP) = - CONJG(CH(1,EL))
         MAT(EL+NEQS ,NP) = - CONJG(CH(2,EL))
      ENDIF
         DO 100 I=1,NR
100        IF (SHOW.GE.5) WRITE(KO,105) I,(MAT(I,J),J=1,NP)
105      FORMAT(' MAT(',I2,',*) =', /(6(1P,E10.2,E9.1,1X)))
C
      CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
         IF(SING) GO TO 600
C
      T = EXP(DBLE(S) * 0.5/NEQS)
      IF(.NOT.REPEAT) RENORM = RENORM / T
       IF(RENORM.GT.BIG) RENORM = BIG
       IF(RENORM.LT.EPS) RENORM = EPS
      IF(SHOW.GE.3) WRITE(KO,130) (MAT(I,NP),I=1,NEQS)
130   FORMAT(/' The combining coefficients are',1P,/(2(E20.8,E15.8)))
C
c     ENDIF
c                  end of FJSWITCH loop matching CRCWFN to zero
      DO 150 K=1,NEQS
      SMAT(K) = MAT(K+NEQS,NP)
C
      IF(ABS(SMAT(K)).GT.1E-10) THEN
        T = 0.0
	WVD = 0.0
       DO 149 IT=1,NEQS
         IF(MAT(IT,NP).EQ.CZ) GO TO 149
            WVD = WVD + MAT(IT,NP)*FIMD(IT,K)
            T = MAX(T,ABS(WVD))
149      CONTINUE
             H2C(K) = T / (ABS(WVD)+1D-20)
           if(ECM(K)>0.0) then
             AL = MAX(AL, H2C(K))
             ALUN = MAX(ALUN, H2C(K)*ABS(SMAT(K)))
	   endif
!	     write(106,*) ALUN,K,WVD,H2C(k),smat(K),AL
      ENDIF
150      CONTINUE
      IF(SHOW.ge.1) WRITE(KO,1157) S,AL
1157  FORMAT(' ERWINCC: LOG10(DETERMINANT) =',2F11.3,', ACCURACY LOSS ='
     X  , 1P,E10.2)
C                       NOW CHECK THE LOGARITHMIC DERIVS. MATCH REQD.
      DO 200 K=1,NEQS
      IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).and.final)
     & WRITE(KO,190) K,SMAT(K),L(K),JVAL(K),CORESP(K),L(EL),
     &           LOG10(MAX(H2C(K),1D0))
190   FORMAT( '   S-matrix',I5,' = '
     & , 2F10.5,' for L=',I5,', J=',F7.1,' channel on core I =',F6.1,
     & ' from L=',I5,',  Acc. loss =',F5.1,' D.')
      IF(K.EQ.EL.AND.SMATEL.AND.ABS(SMAT(K)).GT.1.D-6.and.final) THEN
         SE = LOG(SMAT(EL))
         WRITE(KO,195) K,SE*(0.,-.5)*RAD,L(K),JVAL(K)
195      FORMAT( ' Elastic phase shift ',I3,' = '
     &     , 2F8.3,' deg. for the L =',I5,', J =',F7.1,' channel.')
         WRITE(45,196) ECM(EL),SE*(0.,-.5)*RAD,L(K),JVAL(K)
196        format(f10.3,2f9.3,' for LJin =',i6,f6.1)
	 written(45) = .true.
         ENDIF
200   CONTINUE
      if(.not.FJSWTCH) then
C     IF(AL.GT.1D10.OR.AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
      IF(AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
     X          3.*ALUN*RE * 100.
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWINCC =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(ALUN*RE.GT..03) then
                WRITE(KO,*) AL,ALUN,RE
                write(KO,*) 'Increase lower cutoff: eg CUTL. STOP'
                stop
          endif

	endif
      RETURN
600   continue
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      CALL ABEND(4)
      END
	SUBROUTINE POTWF(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	IMPLICIT NONE
	REAL*8 R12
	INTEGER NEQS,MAXB,K,J,IT,MAXCH
	COMPLEX*16 COUPL(MAXB,MAXB),ZI(MAXB,MAXCH),FI(MAXB,MAXCH),
     X		   V(MAXB,MAXCH),SE,ZIV
C					SCALAR MATRIX MULTIPLIES
C					NO SKIPS IF ZERO COUPLING

         DO 24 K=1,NEQS
             DO 22 IT=1,NEQS
	     SE = 0.0
	     DO J=1,NEQS
             SE = SE + COUPL(K,J) * ZI(IT,J)
	     ENDDO
 22          FI(IT,K) = FI(IT,K) - SE * R12
 24       CONTINUE
         DO 34 K=1,NEQS
            DO 33 IT=1,NEQS
	    ZIV = 0.
            DO J=1,NEQS
		ZIV = ZIV + COUPL(K,J) * FI(IT,J)
	    ENDDO
33          V(IT,K) = ZIV
34       CONTINUE
	return
	end
	SUBROUTINE POTWF2(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	IMPLICIT NONE
	REAL*8 R12
	INTEGER NEQS,MAXB,K,J,IT,MAXCH
	COMPLEX*16 COUPL(MAXB,MAXB),ZI(MAXB,MAXCH),FI(MAXB,MAXCH),
     X		   V(MAXB,MAXCH),C,ZERO
	PARAMETER(ZERO = (0d0,0d0))
C					VECTORISED MATRIX MULTIPLIES
C					ALLOWS SKIPS IF ZERO COUPLING

         DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
             C = COUPL(K,J) * R12
	      if(C/=ZERO) FI(1:NEQS,K) = FI(1:NEQS,K) - C * ZI(1:NEQS,J)
 24       CONTINUE
         DO 34 K=1,NEQS
	    V(:,K)  = ZERO
            DO 34 J=1,NEQS
	      C = COUPL(K,J)
	      if(C/=ZERO) V(1:NEQS,K) = V(1:NEQS,K) + C * FI(1:NEQS,J)
34       CONTINUE
	return
	end
****DISPX**************************************************************
      SUBROUTINE DISPX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGJ,PART,
     X           EXCIT,JEX,ITC,K,RMASS,PEL,EXL,LVAL,FUSL,OUTJ,
     X           JVAL,JPROJ,JTARG,JTOTAL,PARITY,
     X           SMATS,CHANS,JEND,SCALE,ITCM,IF1,IF2)
	use io
	use parameters
	use searchpar, only: final
	use trace, only: smatl
	use fresco1, only: rela,rener ! AMoro
      IMPLICIT REAL*8(A-H,O-Z)
C
      INTEGER PEL,EXL,PARITY,SMATS,CHANS,C,ITC(MXP,MXX),EL,
     X        PART(NCH),EXCIT(NCH),LVAL(NCH)
      COMPLEX*16 SMAT(NCH)
      REAL*8 XS(MAXCH),SIGJ(0:MXPEX),K(MXP,MXX),SCSIG(0:MXPEX)
      REAL*8 JCOEF,FUSL(1+NFUS)
      REAL*8 RMASS(MXP),JEX(6,MXP,MXX),JTOTAL
      LOGICAL JEND
      REAL*8 JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH)
      CHARACTER*1 BLANK,SLASH,HATCH,PSIGN(3),SCALE(0:MXPEX),
     x		CHF(1:1+NFUS)
C
      DATA  SLASH,HATCH /'/','#'/
      DATA PSIGN / '-','?','+' /, BLANK / ' ' /, Z / 0D0 /
C
      CHF(1) = 'f'
      do i=1,NFUS
      CHF(1+i) = char(ichar('c')+i-1)
      enddo
      FUSL(1) = 0.0
      SIGJ(0) = 0.0
      IFT = 1
      DO 690 C=1,NCH
         IC = PART(C)
         IA = EXCIT(C)
         IT = ITC(IC,IA)
        IF(DBLE(K(IC,IA)**2) .LT. 0.0) GOTO 690
      AMDSQS = ABS(SMAT(C))**2
      IF(C.EQ.EL) AMDSQS = 1 - AMDSQS
!@@
         IF(RMASS(IC).lt.1e-5) then
             S =          RMASS(PEL)*amu/ (HBC*K(PEL,EXL))
!@@@         S = K(IC,IA)*RMASS(PEL)*amu/ (HBC*K(PEL,EXL))
         ELSE IF(RMASS(PEL).lt.1e-5) then
             S = HBC*K(IC,IA)/(RMASS(IC)*amu)
         ELSE
             S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
         ENDIF
! AMoro & MGR
         if ((rela.eq.'3d').or.(rela.eq.'na').or.(rela.eq.'c').or.
     &    (rela.eq.'30').or.(rela.eq.'n0')) then
!         write(0,*) 'Rescaling frxx3.f'
         S=S*(RENER(PEL,EXL)/RENER(IC,IA))
     &      *(RMASS(IC)/RMASS(PEL))
!         write(*,*)'frxx3: xs scaled by',
!     & (RENER(PEL,EXL)/RENER(IC,IA))
!     &      *(RMASS(IC)/RMASS(PEL))
         endif

!	 S is v_out/v_in but for photons this is just
!	      c/v_in for emission so since
!	      v_in/c = h*K_in/(RMASS_in*c) 
!	             = (hc)*K_in/(RMASS_in*c^2)
!		     = (hc)*K_in/(RMASS_in*amu)
!	 S = RMASS(PEL)*amu/ (HBC*K(PEL,EXL))	or as before
!	 
!        S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
!@@
      XS(C) = JCOEF* AMDSQS * S
      IF(C.EQ.EL) THEN
         SIGJ(0) = SIGJ(0) + XS(C)
         FUSL(1) = FUSL(1) + XS(C)
        ELSE
         SIGJ(IT) = SIGJ(IT) + XS(C)
         FUSL(1) = FUSL(1) - XS(C)
        ENDIF
      IFT = IF2
  690 CONTINUE
C
	if(SMATL>=2) 
     X WRITE(38,1445) JTOTAL,PSIGN(PARITY+2),LVAL(EL),JIN
     X       , (SIGJ(IT),IT=0,ITCM),(FUSL(I),I=IF1,IFT)
 1445 FORMAT(F7.1,A1,I6,I2,10G12.4,:,/(16X,10G12.4))
	OUTJ = SIGJ(0)
	if(SMATL>=2) then
       	written(38) = .true.
       	call flush(38)
	endif
	SCSIG(:) = SIGJ(:)
      IF(final.and.(CHANS+SMATS.GE.1.OR.JEND)) THEN
      DO 710 IT=0,ITCM
      SCALE(IT) = BLANK
      IF(ABS(SCSIG(IT)).GT.0.010) GO TO 710
         SCSIG(IT) = SCSIG(IT) * 1000.0
         SCALE(IT)= SLASH
      IF(ABS(SCSIG(IT)).GT.0.010) GO TO 710
         SCSIG(IT) = SCSIG(IT) * 1000.0
         SCALE(IT)= HATCH
  710 CONTINUE
       WRITE(KO,1450) JTOTAL,PSIGN(PARITY+2),MOD(LVAL(EL),100),JIN
     X    , (SCSIG(IT),SCALE(IT),IT=0,ITCM),(FUSL(I),CHF(I),I=IF1,IFT)
       ENDIF
 1450 FORMAT(' Reaction Xsec',F8.1,A1,'/',I2,' @',I2,' =',F8.3,A1,',',
     X   ' Out:', 9(F8.3,A1),:,/(11X,'Xsec',23X,10(F8.3,A1)))
      if(SMATL>1) write(KO,*)
      CALL FLUSH(KO)
C
      RETURN
      END
      function lnbl(s)
      character*70 s
      l=len(s)
      do 1 i=l,1,-1
      if(s(i:i).ne.' ') go to 5
1     continue
      i=0
      i=1
5     lnbl=i
      end
      SUBROUTINE FILEMV(KIN,KOUT)
      CHARACTER*133 ALINE
	write(51,*) ' Append all file ',KIN,' to ',KOUT
      WRITE(KIN,713)
  713    FORMAT('EOF')
      REWIND KIN
  777 READ(KIN,'(a)',END=778) ALINE
      IF(ALINE(1:3).eq.'EOF') GO TO 778
      do 1 l=133,1,-1
      if(ALINE(l:l).ne.' ') go to 5
1     continue
      l=0
5     WRITE(KOUT,'(a)') ALINE(1:l)
      GO TO 777
  778 CALL FLUSH(KOUT)
      REWIND KIN
      end
	subroutine rewop(ifile)
	integer ifile
	logical op
	character*8 name
	inquire(ifile,opened=op)
	if(.not.op) then
	   if(ifile<10) then
		write(name,'(''fort.'',i1)') ifile
	   else if(ifile<100) then
		write(name,'(''fort.'',i2)') ifile
	   else if(ifile<1000) then
		write(name,'(''fort.'',i3)') ifile
	   endif
	   open(ifile,form='formatted',file=name)
	endif
	rewind ifile
	return
	end

	subroutine openif(ifile)
	integer ifile
	logical op
	character*8 name
	inquire(ifile,opened=op)
	if(.not.op) then
	   if(ifile<10) then
		write(name,'(''fort.'',i1)') ifile
	   else if(ifile<100) then
		write(name,'(''fort.'',i2)') ifile
	   else if(ifile<1000) then
		write(name,'(''fort.'',i3)') ifile
	   endif
	   open(ifile,form='formatted',file=name)
	endif
!!!!!!	rewind ifile
	return
	end

c *** -----------------------------------------------------------
c *** Re-orthogonalize solutions by LQ factorization (AMoro)
c *** -----------------------------------------------------------
      subroutine LQFACT(ZI,ZM,FI,WDIAG,COUPL,R12,NEQS,MAXB,MAXCH)
      implicit none
      integer   :: neqs,maxb,maxch,it,itp,k
      complex*16:: zi(maxb,maxch),zm(maxb,maxch),fi(maxb,maxch)
      complex*16:: coupl(maxb,maxb),c
      complex*16:: wdiag(neqs),w(maxb,maxch)
      real*8    :: r12
      complex*16,parameter:: zero=(0d0,0d0)
c ... for ztrsm
      complex*16 alpha
      integer ldb
      character diag,side,transa,uplo
c ... for ZGELQF
      complex*16:: LT(maxb,maxch),work(2*maxb),tau(min(maxb,maxch))
      integer   :: nrhs, info,lwork
      character*1  trans
      integer   :: ipiv(neqs)
!      EXTERNAL  ZGETRF, ZGETRS

c ... LQ factorization: ZI=LT*Q 
c     (Q=orthogonal matrix; LT= triangular LOWER matrix)
      LT=ZI
!      LWORK=-1
!      call ZGELQF(neqs,neqs, LT, MAXB, TAU, WORK, LWORK, INFO )
!      LWORK=WORK(1)
      LWORK=2*maxb
      call ZGELQF(neqs,neqs, LT, MAXB, TAU, WORK, LWORK, INFO )
      if (info.ne.0)  then
         write(*,*)'zgelqf failed with exit code=',info
      endif
 
c ... Calculate new values of ZI, ZM by solving:
c       Z(old)= LT*Z(new)  
c     Solve   A*X  = B
c     and make: Z(new)= X 
      
c ... X*op( A ) = alpha*B
      side  ='L' ! A acts on the left
      uplo  ='L' ! R is lower triangle matrix
      transa='N' ! no transpose
      diag  ='N' !??????????
      alpha =1d0

      if (1>2) then 

      write(*,*)'Re[zi(old)] '
      do k=1,min(5,neqs)
        write(*,'(5x,i3,50g14.5)') k,
     &  (real(zi(it,k)), it=1,min(5,neqs))
      enddo
      write(*,*)' '


      write(*,*)'|zi x zip| (before QR)'
        do it=1,min(neqs,5)
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(zi(it,:),zi(itp,:)))
     &      / sqrt(abs(dot_product(zi(it,:),zi(it,:))))
     &      / sqrt(abs(dot_product(zi(itp,:),zi(itp,:)))),
     &   itp=1,min(neqs,5)) 
        enddo
      endif
      
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,NEQS,NEQS,ALPHA,LT,MAXB,ZI,MAXB)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,NEQS,NEQS,ALPHA,LT,MAXB,ZM,MAXB)

      if (1>2) then !
      write(*,*)'|zi x zip| (after QR)'
        do it=1,min(neqs,5)
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(zi(it,:),zi(itp,:)))
     &      / sqrt(abs(dot_product(zi(it,:),zi(it,:))))
     &      / sqrt(abs(dot_product(zi(itp,:),zi(itp,:)))), 
     &    itp=1,min(neqs,5)) 
        enddo
      endif


c Recalculate FI(:,:) 
      FI=zero
      do 24 k=1,neqs
       do 24 it=1,neqs
       if (it.eq.k) fi(k,k)=fi(k,k)+wdiag(k)
       c=coupl(k,it)*r12
       if (c/=zero) fi(1:neqs,k)=fi(1:neqs,k)-c*zi(1:neqs,it)
24    continue
      end subroutine        



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


******FRXX4*************************************************************
*****CPAIR**************************************************************
      SUBROUTINE CPAIR(CP,IC1,IC2,KIND,QQ,IP3,REV,LOCAL,ICOM,ICOR,
     &           LOCF,NG,GPT,XA,XB,XP,XQ,IREM,FILE,CHNO,VREAL,KLT,
     &           NCH,LVAL,JVAL,JPROJ,JTARG,JTOTAL,PART,EXCIT,JEX,COPY,
     &           BAND,QNF,AFRAC,CUTOFF,ALOSSM,RERR,FPT,NLL,NK,BE,
     &           CLIST,NFLIST,NCLIST,DNL,EPS,NFORML1,NFORML2,
     &           HP,ITC,NEX2,ICFROM,ICTO,MATRIX,MEK,CPSO,SPINTR,K,
     & 		 POTCAP,
     X           LTRANS,NLOC,LCALL,MCALL,FORML,VFOLD,VCORE,NIB,MASS,
     X		 FORMC)
	use parameters
        use factorials
	use trace
	use kcom
	use drier
	use io
	use fresco1, only: cxwf,ccbins,sumccbins,sumkql
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER CP,QQ,Q,ICOM(MAXCPL,2),ICOR(MAXCPL,2),GPT(2,MAXQRN),NG(2),
     &        FILE,CHNO(MFNL,6),LVAL(NCH),PART(NCH),EXCIT(MAXCH,3),
     &        COPY(2,MXP,MXX),QNF(19,MSP),CUTOFF,C1,C2,C2P,PARITY,
     &        FPT(6,MAXQRN),ITC(MXP,MXX),MATRIX(6,MPAIR),C2LAST,
     &        NCLIST(MAXCH,MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),NC,QC,LA
      INTEGER PARITYP,BAND(2,MXP,MXX)
      INTEGER KMMI(3),KMMF(3),NFS(MLOC),POTCAP(MLOC,2)
      LOGICAL REV,LOCAL,VREAL,REPEAT,LCL,MCL,R1DONE,REO,FAIL3,FRAC,
     X        C1FR,LTRANS(MAXQRN),LCALL,MCALL,CPSO,CPLD(NCH,NCH,0:4),
     X        LCLA,MCLA
      REAL*8 JVAL(NCH),JPROJ(NCH),JTARG(NCH),JTOTAL,JEX(6,MXP,MXX),JJI,
     &       KG,JN,JNP,DNL(MAXNLO),ICORE,ICOREP,KCORE,KCOREP,
     &       JA,HP(MXP),JAMIN,JAMAX,JI,JJN,JCOM,JCOMP
      REAL*8 AFRAC(MXPEX,MXPEX,2,MSP),MEK(MPAIR),SPINTR(2,MPAIR),
     X       NLOC(NLM),FORML(NFORML1,NFORML2,2),VSOCC(NCH,NCH,2),
     X	     COUPI(NCH,NCH),COUPF(NCH,NCH),JNPP,BE(MSP,4),MASS(4,MXP+1)
      REAL*8 K(MXP,MXX)
      COMPLEX*16 CLIST(MAXCH,MAXCH,MCLIST),C6,CH,
     X           FORMC(MAXNLC,MSP,2),
     X          VFOLD(NLN),VCORE(NLN),PH12!,FORMF(MAXM,MLOC)
      REAL*8,allocatable:: MCG(:,:,:,:),KCOEF(:,:),ALOSS(:),THM(:)
      COMPLEX*16,allocatable:: FNC(:,:)
      INTEGER,allocatable:: NM(:,:)
      CHARACTER PSIGN(3)
      DATA PSIGN / '-','?','+' /
      DATA ONE,TWO /1D0, 2D0/
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-NINT(X)).GT.1D-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
      SL2(JN,LN,SN) = JN*(JN+1.) - LN*(LN+1.) - SN*(SN+1.)
      Z = 0.0
      Q = ABS(QQ)
      ICH = 0
      allocate(KCOEF(MAXQRN,KLT),ALOSS(MAXQRN),THM(MAXNLN))
      allocate(FNC(NLN,NLO))
      GO TO (90,100,30,30,50,50,50,50,90,100),KIND
C FOR KIND =  1   2  3  4  5  6  7  8  9, 10
C
C    SINGLE-PARTICLE INELASTIC EXCITATIONS (KIND=3 PROJ;  KIND=4 TARG)
C    .................................................................
C
30    IF(ICFROM.EQ.ICTO) GO TO 200
      IN = KIND-2
      IC1 = ICOM(CP,IN)
      IC2 = ICOR(CP,IN)
      IG1=0

      DO 45 IG=1,NK

         KN  = FPT(2,IG)
         KNP = FPT(1,IG)
         KM  = IG + LOCF
         IK = FPT(3,IG)
         KSO= FPT(4,IG)
         QC = FPT(5,IG)
         LA = FPT(6,IG)
         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
         SN = 0.5*QNF(10,KN)
         JN = 0.5*QNF(11,KN)
         JNP= 0.5*QNF(11,KNP)
         IACORE = QNF(3,KN)
         IACOREP= QNF(3,KNP)
         ICORE  = JEX(IN,IC2,IACORE)
         ICOREP = JEX(IN,IC2,IACOREP)
         KCORE  = JEX(IN+2,IC2,IACORE)
         KCOREP = JEX(IN+2,IC2,IACOREP)

         IF(LISTCC.GE.2)then
          if(ccbins)then
           WRITE(KO,309) IG,KNP,KN,KM,IK,KSO,QC,LA
            IF(LISTCC.GE.3) WRITE(KO,306) 'l',LN,JN,ICORE
            IF(LISTCC.GE.3) WRITE(KO,306) 'r',LNP,JNP,ICOREP
          else
           WRITE(KO,304) IG,KNP,KN,KM,IK,KSO
            IF(LISTCC.GE.3) WRITE(KO,303) 'l',LN,JN
            IF(LISTCC.GE.3) WRITE(KO,303) 'r',LNP,JNP
          endif
         ENDIF
304       format('  SP pair#',i8,':',2i4,'@',i8,'; K=',i2,', Vso#',i2)
309       format('  SP pair#',i8,':',2i4,'@',i8,'; K=',i2,', Vso#',i2,
     &           ', QC=',i2,', lambda=',i2)
303       format('     ',a1,'hs state lj =',i3,f5.1)
306       format('     ',a1,'hs state ljI =',i3,2f5.1)       
!       ENDIF
       
      COUPF(:,:) = 0.
      
!***********************************************************************
      IF(KSO==0) then
       COUPI(:,:) = 0.
      R13 = WIG3J(IK+Z,LN+Z,LNP+Z,Z,Z,Z)
      DO 1025 C1=1,NCH
	C2LAST=NCH
	if(REV.and..not.CPSO.and..not.cxwf) C2LAST=C1
      DO 1025 C2=1,C2LAST
        IF(IC1.NE.PART(C1).OR.IC1.NE.PART(C2)) GO TO 1025
         IF(.NOT.REV   .AND.EXCIT(C1,IN+1).LT.EXCIT(C2,IN+1)) GO TO 1025
         IF(MOD(IP3,10)==2.AND.EXCIT(C1,IN+1)/=EXCIT(C2,IN+1))goto 1025
!
!  Short cut if only couplings to or from the ground state, & NOT gs diag:
         IF(MOD(IP3,10).eq.3 .and.
     X     .not.((1.eq.EXCIT(C2,IN+1).or.1.eq.EXCIT(C1,IN+1)).and.
     X     .not.(1.eq.EXCIT(C2,IN+1).and.1.eq.EXCIT(C1,IN+1)))) GOTO1025
!
!  Short cut if only couplings to or from the ground state or diag:
         IF(MOD(IP3,10).eq.4 .and.
     X     .not.(1.eq.EXCIT(C2,IN+1).or.1.eq.EXCIT(C1,IN+1)
     X           .or.EXCIT(C2,IN+1).eq.EXCIT(C1,IN+1)) ) GO TO 1025
!
!  Short cut if only couplings to or from any bound state or diag:
         IF(MOD(IP3,10).eq.5 .and.
     X     .not.(BE(KNP,1)>0.0.or.BE(KN,1)>0.0
     X           .or.EXCIT(C2,IN+1).eq.EXCIT(C1,IN+1)) ) GO TO 1025
!
!  Short cut if to remove couplings between bins with different J values in projectile
!	 IF(LISTCC>0.and.MOD(IP3,10).eq.6 .and. KIND==3) then
! 307	 format(' KNP,KN, BEP,BE, JP1,JP2, test:',2i4,2f8.4,2f5.1,L4)
!  	 endif
         IF(MOD(IP3,10).eq.6 .and. KIND==3 .and.
     X          (BE(KNP,1)<0.0.and.BE(KN,1)<0.0
     X           .and.ABS(JPROJ(C1)-JPROJ(C2)).GT.1E-5) ) GO TO 1025
!
!  Short cut if to remove couplings between bins with different J values in target
         IF(MOD(IP3,10).eq.6 .and. KIND==4 .and.
     X          (BE(KNP,1)<0.0.and.BE(KN,1)<0.0
     X           .and.ABS(JTARG(C1)-JTARG(C2)).GT.1E-5) ) GO TO 1025

         IF(KIND.EQ.3.AND.ABS(JTARG(C1)-JTARG(C2)).GT.1E-5) GO TO 1025
         IF(KIND.EQ.3.AND.ABS(JVAL(C1)-JVAL(C2)).GT.1E-5) GO TO 1025
         IF(KIND.EQ.4.AND.ABS(JPROJ(C1)-JPROJ(C2)).GT.1E-5) GO TO 1025
	 INO = 3-IN ! other nucleus
         IF(EXCIT(C1,INO+1).ne.EXCIT(C2,INO+1)) GO TO 1025 ! should be same state, not just same spin
C      COUPLING FROM FPT(2 = KNP = C2
C               TO   FPT(1 = KN  = C1
           IB1 = EXCIT(C1,1)
           IB2 = EXCIT(C2,1)
         IF(COPY(IN,IC1,IB1 ).NE.0) IB1  = COPY(IN,IC1,IB1 )
         IF(COPY(IN,IC1,IB2 ).NE.0) IB2  = COPY(IN,IC1,IB2 )
         REO = MOD(IP3,10).eq.1.AND.EXCIT(C1,IN+1).EQ.EXCIT(C2,IN+1)
	 L12 = LVAL(C2)+LVAL(C1)

	 if(mod(IK+L12,2).eq.1) go to 1025
         IF(REO.and.IK.gt.0) GO TO 1025
          IF(REO) GO TO 1025  	! Uncomment this line, to allow removal
!				of reorientation for monopoles.
	 if(QNF(10,KN).ne.QNF(10,KNP)) go to 1025

	 R1DONE=.false.
      S = 0.0
       DO 1023 IACORE=1,NEX2
       R6 = AFRAC(ITC(IC1,IB1),ITC(IC2,IACORE),IN,KN)
         IF(ABS(R6) .LT. EPS) GO TO 1023
       R6P= AFRAC(ITC(IC1,IB2),ITC(IC2,IACORE),IN,KNP)
        IF(LISTCC>3) then
!	  write(KO,*) ' IACORE, KN,KNP,R6,R6P =',IACORE, KN,KNP,R6,R6P
          write(KO,310) 'TO',ITC(IC1,IB1),ITC(IC2,IACORE),IN,KN,R6
          write(KO,310) 'FROM',ITC(IC1,IB2),ITC(IC2,IACORE),IN,KNP,R6P
310       format('  ',a4,' AFRAC(',4i3,') =',f10.5)
	endif
         IF(ABS(R6P).LT.EPS) GO TO 1023
          IF(COPY(IN,IC2,IACORE).NE.0) GO TO 1023
	if(.not.R1DONE) then
         ICORE = JEX(IN,IC2,IACORE)
         R14 = WIG3J(IK+Z,LVAL(C1)+Z,LVAL(C2)+Z,Z,Z,Z)
     &          / SQRT(2 * IK + 1.0)
          IF(ABS(R14) .LT. EPS) GO TO 1025
         R12=SQRT((2*LN+1.)*(2*LVAL(C1)+1.)*(2*LNP+1.)*(2*LVAL(C2)+1.))
          R1 = R12 * R13 * R14
           LTMIN = MAX(ABS(LVAL(C1)-LN), ABS(LVAL(C2)-LNP))
           LTMAX = MIN(    LVAL(C1)+LN ,     LVAL(C2)+LNP )
	 R1DONE=.true.
	 endif
C            IF(QNF(12,KNP).GT.1) STOP 'SP/DP'
             IF(QNF(12,KNP).GT.1) CALL ABEND(32)
         DO 1022 LTOTAL=LTMIN,LTMAX
         R21 = (-1)**(LTOTAL+IK) * (2*IK+1)
         R22 = RACAH(LN+Z,LNP+Z,LVAL(C1)+Z,LVAL(C2)+Z,IK+Z,LTOTAL+Z)
         R2 = R6 * R6P * R21 * R22
      IF(KIND.EQ.3) THEN
C          KIND = 3 ; COUPLING BETWEEN JPROJ(C1) & JPROJ(C2)
          FMIN = MAX(ABS(SN-ICORE),ABS(LN-JPROJ(C1)),ABS(LNP-JPROJ(C2)),
     &               ABS(LTOTAL-JVAL(C1)),ABS(LTOTAL-JVAL(C2)))
          FMAX = MIN(    SN+ICORE ,    LN+JPROJ(C1) ,    LNP+JPROJ(C2) ,
     &                   LTOTAL+JVAL(C1) ,    LTOTAL+JVAL(C2) )
!         DO 1019 F=FMIN,FMAX
        NF = NINT(FMAX-FMIN)
        DO 1019 IIF=0,NF
        F=FMIN+IIF
           R31 = SQRT((2*JN+1)*(2*JNP+1)) * (2*F+1)
     &           * RACAH(LN+Z,SN,JPROJ(C1),ICORE,JN,F)
     &           * RACAH(LNP+Z,SN,JPROJ(C2),ICORE,JNP,F)
           R32 = SQRT((2*JPROJ(C1)+1)*(2*JPROJ(C2)+1)) * (2*LTOTAL+1)
     &           * RACAH(LVAL(C1)+Z,LN +Z,JVAL(C1),F,LTOTAL+Z,JPROJ(C1))
     &           * RACAH(LVAL(C2)+Z,LNP+Z,JVAL(C2),F,LTOTAL+Z,JPROJ(C2))
           R3 = R31 * R32
           T = R1 * R2 * R3
           S = S + T
      IF(LISTCC.GE.2) WRITE(KO,2185) C1,C2,KNP,KN,KM,IG,
     X                     IK,KSO,LN,LVAL(C1),LVAL
     &(C2), ICORE,F,LTOTAL,R6,R6P,R1,R2,R3,T,S
2185   FORMAT(/ 1X,4I3,2i5,5i3,2F5.1,I4,2F8.4,4F9.5,F10.5)
1019         CONTINUE
      ELSE
C          KIND = 4 ; COUPLING BETWEEN JTARG(C1) & JTARG(C2)
       JAMIN = MAX(ABS(JVAL(C1)-JN),ABS(JTOTAL-ICORE),ABS(JVAL(C2)-JNP))
       JAMAX = MIN(    JVAL(C1)+JN ,    JTOTAL+ICORE ,    JVAL(C2)+JNP )
!       DO 1020 JA=JAMIN,JAMAX
       NJA = NINT(JAMAX-JAMIN)
       DO 1020 IJA=0,NJA
       JA=JAMIN+IJA
       R31 = SQRT((2*JTARG(C1)+1)*(2*JTARG(C2)+1)) * (2*JA+1)
     &      * RACAH(JVAL(C1),JN ,JTOTAL,ICORE,JA,JTARG(C1))
     &      * RACAH(JVAL(C2),JNP,JTOTAL,ICORE,JA,JTARG(C2))
          SAMIN = MAX(ABS(JPROJ(C1)-SN),ABS(LTOTAL-JA))
          SAMAX = MIN(    JPROJ(C1)+SN ,    LTOTAL+JA )
!          DO 1020 SA=SAMIN,SAMAX
          NSA=NINT(SAMAX-SAMIN)
          DO 1020 ISA=0,NSA
          SA=SAMIN+ISA
          R32 = SQRT((2*JVAL(C1)+1)*(2*JN+1) * (2*LTOTAL+1)*(2*SA+1))
     &          * WIG9J(LVAL(C1)+Z, LN+Z, LTOTAL+Z,
     &                  JPROJ(C1),  SN,   SA,
     &                  JVAL(C1),   JN,   JA)
          R33 = SQRT((2*JVAL(C2)+1)*(2*JNP+1) * (2*LTOTAL+1)*(2*SA+1))
     &          * WIG9J(LVAL(C2)+Z, LNP+Z,LTOTAL+Z,
     &                  JPROJ(C2),  SN,   SA,
     &                  JVAL(C2),   JNP,  JA)
           R3 = R31 * R32 * R33
           T = R1 * R2 * R3
           S = S + T
      IF(LISTCC.GE.2) WRITE(KO,195) C1,C2,KNP,KN,KM,IG,
     X IK,LN,LVAL(C1),LVAL(C2), 
     X ICORE,JA,SA,LTOTAL,R6,R6P,R1,R2,R3,T,S
195   FORMAT(/ 1X,10I3,3F5.1,I4,2F8.4,4F9.5,F10.5)
1020       CONTINUE
      ENDIF
1022  CONTINUE
1023  CONTINUE
                                       COUPI(C1,C2) = S
      if(REV.and..not.cxwf.and.C1/=C2) COUPI(C2,C1) = S

      IF(LISTCC.GE.1.AND.ABS(S).GT.EPS)
     & WRITE(KO,24) C1,C2,KNP,KN,KM,IG,IK,KSO,LN,LVAL(C1),LVAL(C2),R1,S
24    FORMAT(/ 1X,4I4,2I5,3I2,2I5,F10.5,F12.5)
1025  CONTINUE
	COUPF(:,:) = COUPI(:,:)


!  QC>0 CODE DELETED IN THIS VERSION
!***********************************************************************
!***********************************************************************
      ELSE if(KSO==1) then   	! - Vso(opt)
	DO 26 C1=1,NCH
	DO 26 C2=1,NCH
26	COUPF(C1,C2) = COUPI(C1,C2) * SL2(JVAL(C2),LVAL(C2),JPROJ(C2))

!***********************************************************************
      ELSE if(KSO==2) then   	! + Vso(frag) lv*s gives L.s and l.s

	 VSOCC(:,:,:) = 0.0
        DO 29 KNPP=1,MSP
         LNPP= QNF(9,KNPP)
         JNPP= 0.5*QNF(11,KNPP)
	 if(LNP/=LNPP) go to 29  
! This next line is a crude orthogonality test
! Assumes all jJ bins for same l,I have the same mid-point energies.
! and assumes no jJ-dependent potentials in the bins.
         if(abs(BE(KNP,1)-BE(KNPP,1))>0.01) go to 29
      	   IF(LISTCC.GE.2) WRITE(KO,*) 'SO state pair:',KNP,KNPP
      	   IF(LISTCC.GE.3) WRITE(KO,*) 'lhs state lj =',LNP,real(JNP)
      	   IF(LISTCC.GE.3) WRITE(KO,*) 'rhs state lj =',LNPP,real(JNPP)

	DO 27 C2P=1,NCH
         IB2P = EXCIT(C2P,1)
         IF(COPY(IN,IC1,IB2P).NE.0) IB2P = COPY(IN,IC1,IB2P)
	DO 27 C2=1,NCH
         IB2 = EXCIT(C2,1)
         IF(COPY(IN,IC1,IB2 ).NE.0) IB2  = COPY(IN,IC1,IB2 )
         DO 27 IACORE=1,NEX2
          IF(IC1.NE.PART(C2)) GO TO 27
          IF(ABS(JVAL(C2P)-JVAL(C2)).GT.1E-5) GO TO 27
          R6P= AFRAC(ITC(IC1,IB2P),ITC(IC2,IACORE),IN,KNP)
             IF(ABS(R6P).LT.EPS) GO TO 27
          R6PP = AFRAC(ITC(IC1,IB2 ),ITC(IC2,IACORE),IN,KNPP)
             IF(ABS(R6PP) .LT. EPS) GO TO 27
             IF(COPY(IN,IC2,IACORE).NE.0) GO TO 27
         ICORE = JEX(IN,IC2,IACORE)

	 if(C2==C2P.and.KNP==KNPP) then
 	   VSOCC(C2P,C2,2) = SL2(JNP,LNP,SN) 		! l.s
    	   if(LISTCC.ge.2) write(KO,261) IG,KNP,KNPP,IK,C2P,C2,
     X			VSOCC(C2P,C2,2)
261     format(1x,i5,2i3,i2,2i4,' l.s VSOCC =',f12.6)
	   endif
								! now do L.s
	 LVA = LVAL(C2)
	 if(LVAL(C2P).ne.LVA) go to 27
!	if(LISTCC>4)  write(KO,*) 
!     x	  ' C2P,C2,ICORE,LNP,JVAL(C2) = ',C2P,C2,ICORE,LNP,JVAL(C2)
	
	JJN2min = nint(2*abs(LVA-SN))
	JJN2max = nint(2*(LVA+SN))
	if(LISTCC>4) then
	  write(KO,*) '      JJN from ',real(JJN2min)/2,real(JJN2max)/2
	  endif
	DO JJN2 = JJN2min,JJN2max,2
	 JJN = JJN2*0.5
	 JI2min = max(nint(2*abs(LNP-ICORE)),nint(2*abs(JVAL(C2)-JJN)),
     x	                                     nint(2*abs(JVAL(C2P)-JJN)))
	 JI2max = min(nint(2*   (LNP+ICORE)),nint(2*   (JVAL(C2)+JJN)),
     x	                                     nint(2*   (JVAL(C2P)+JJN)))
	if(LISTCC>4) then
	  write(KO,*) ' JJN = ',real(JJN)
	  write(KO,*) '      JI from ',real(JI2min)/2,real(JI2max)/2
	  endif
      	 VSO = SL2(JJN,LVA,SN) 
	DO JI2 = JI2min,JI2max,2
	 JI = JI2*0.5
	if(LISTCC>4) write(KO,*) ' JI = ',real(JI)
	 R1 = sqrt((JJN2+1.)*(2*JPROJ(C2)+1.))* 
     X		RACAH(JI,SN,JVAL(C2),LVA+Z,JPROJ(C2),JJN)
	 R2 = sqrt((JI2 +1.)*(2*JNPP+1.))* 
     X		RACAH(SN,LNPP+Z,JPROJ(C2),ICORE,JNPP,JI)
	 R1P= sqrt((JJN2+1.)*(2*JPROJ(C2P)+1.))* 
     X		RACAH(JI,SN,JVAL(C2P),LVA+Z,JPROJ(C2P),JJN)
	 R2P= sqrt((JI2 +1.)*(2*JNP+1.))* 
     X		RACAH(SN,LNP+Z,JPROJ(C2P),ICORE,JNP,JI)
	 T = R6P*R6PP * R1*R1P*R2*R2P * VSO  
	 VSOCC(C2P,C2,1) = VSOCC(C2P,C2,1) + T
	if(LISTCC.ge.3) write(KO,268) IG,KNP,KNPP,IK,C2P,C2,JJN,JI,
     X			R6P,R6PP,R1,R1P,R2,R2P,VSO,T,VSOCC(C2P,C2,1)
268      format(1x,i5,2i3,i2,2i4,f7.1,f5.1,2f6.3,6f8.4,
     X			' L.s VSOCC=',f9.5)
	 enddo
	 enddo
    	if(LISTCC.ge.2) write(KO,269) IG,KNP,KNPP,IK,C2P,C2,
     X			VSOCC(C2P,C2,1)
269     format(1x,i5,2i3,i2,2i4,' L.s VSOCC =',f12.6)
27	continue
29	continue
	DO 31 C1=1,NCH
	DO 31 C2=1,NCH
	DO 31 C2P=1,NCH
31	COUPF(C1,C2)=COUPF(C1,C2) 
     x                  + COUPI(C1,C2P) * VSOCC(C2P,C2,1)*(-XQ)  !L.s
     x                  + COUPI(C1,C2P) * VSOCC(C2P,C2,2)*  XB   !l.s

!***********************************************************************
      ELSE if(KSO==3) then   	! + Vso(core)

	 VSOCC(:,:,:) = 0.0
        DO 38 KNPP=1,MSP
         LNPP= QNF(9,KNPP)
         JNPP= 0.5*QNF(11,KNPP)
	 if(LNP/=LNPP) go to 38  
! This next line is a crude orthogonality test
! Assumes all jJ bins for same l,I have the same mid-point energies.
! and assumes no jJ-dependent potentials in the bins.
         if(abs(BE(KNP,1)-BE(KNPP,1))>0.01) go to 38
           IF(LISTCC.GE.2) WRITE(KO,*) 'SO state pair:',KNP,KNPP
           IF(LISTCC.GE.3) WRITE(KO,*) 'lhs state lj =',LNP,real(JNP)
           IF(LISTCC.GE.3) WRITE(KO,*) 'rhs state lj =',LNPP,real(JNPP)

	DO 37 C2P=1,NCH
         IB2P = EXCIT(C2P,1)
         IF(COPY(IN,IC1,IB2P).NE.0) IB2P = COPY(IN,IC1,IB2P)
	DO 37 C2=1,NCH
         IB2 = EXCIT(C2,1)
         IF(COPY(IN,IC1,IB2 ).NE.0) IB2  = COPY(IN,IC1,IB2 )
         DO 37 IACORE=1,NEX2
          IF(IC1.NE.PART(C2)) GO TO 37
          IF(ABS(JVAL(C2P)-JVAL(C2)).GT.1E-5) GO TO 37
          R6P  = AFRAC(ITC(IC1,IB2P),ITC(IC2,IACORE),IN,KNP)
             IF(ABS(R6P).LT.EPS) GO TO 37
          R6PP = AFRAC(ITC(IC1,IB2 ),ITC(IC2,IACORE),IN,KNPP)
             IF(ABS(R6PP) .LT. EPS) GO TO 37
             IF(COPY(IN,IC2,IACORE).NE.0) GO TO 37
         ICORE = JEX(IN,IC2,IACORE)
	 LVA = LVAL(C2)
	 if(LVAL(C2P).ne.LVA) go to 37
	
								! now l.I
	JI2min = nint(2*abs(LNP-ICORE))
	JI2max = nint(2*(LNP+ICORE))
	if(LISTCC>4) then
	  write(KO,*) '      JI from ',real(JI2min)/2,real(JI2max)/2
	  endif
          IF(ABS(JPROJ(C2P)-JPROJ(C2))<1E-5) then
	DO JI2 = JI2min,JI2max,2
	 JI = JI2*0.5
	if(LISTCC>4) write(KO,*) ' JI = ',real(JI)
      	 VSO = SL2(JI,LNP,ICORE) 
	 R1PP = sqrt((JI2+1.)*(2*JNPP+1.))* 
     X		RACAH(SN,LNPP+Z,JPROJ(C2),ICORE,JNPP,JI)
	 R1P= sqrt((JI2+1.)*(2*JNP+1.))* 
     X		RACAH(SN,LNP+Z,JPROJ(C2P),ICORE,JNP,JI)
	 T = R6P*R6PP * R1P*R1PP * VSO  
	 VSOCC(C2P,C2,2) = VSOCC(C2P,C2,2) + T
	if(LISTCC.ge.3) write(KO,338) IG,KN,KNP,IK,C2P,C2,JI,
     X			R6P,R6PP,R1P,R1PP,VSO,T,VSOCC(C2P,C2,2)
338      format(1x,i5,2i3,i2,2i4,f7.1,2f6.3,4f8.4,' l.I VSOCC=',f9.5)
	 enddo
	 endif
    	if(LISTCC.ge.2) write(KO,339) IG,KNP,KNPP,IK,C2P,C2,
     X			VSOCC(C2P,C2,2)
339     format(1x,i5,2i3,i2,2i4,' l.I VSOCC =',f12.6)
34	continue


								! now L.I
	JJI2min = nint(2*abs(LVA-ICORE))
	JJI2max = nint(2*(LVA+ICORE))
	if(LISTCC>4) then
	  write(KO,*) '      JJI from ',real(JJI2min)/2,real(JJI2max)/2
	  endif
	DO JJI2 = JJI2min,JJI2max,2
	 JJI = JJI2*0.5
	if(LISTCC>4) write(KO,*) ' JJI = ',real(JJI)
      	 VSO = SL2(JJI,LVA,ICORE) 
	 R1PP = sqrt((JJI2+1.)*(2*JPROJ(C2)+1.))* 
     X		RACAH(LVA+Z,ICORE,JVAL(C2),JNPP,JJI,JPROJ(C2))
	 R1P= sqrt((JJI2+1.)*(2*JPROJ(C2P)+1.))* 
     X		RACAH(LVA+Z,ICORE,JVAL(C2P),JNP,JJI,JPROJ(C2P))
	 T = R6P*R6PP * R1P*R1PP * VSO  
	 VSOCC(C2P,C2,1) = VSOCC(C2P,C2,1) + T
	if(LISTCC.ge.3) write(KO,368) IG,KN,KNP,IK,KSO,C2P,C2,JJI,
     X			R6P,R6PP,R1P,R1PP,VSO,T,VSOCC(C2P,C2,1)
368      format(1x,i5,2i3,2i2,2i4,f7.1,2f6.3,4f8.4,' L.I VSOCC=',f9.5)
	 enddo
    	if(LISTCC.ge.2) write(KO,369) IG,KNP,KNPP,IK,KSO,C2P,C2,
     X			VSOCC(C2P,C2,1)
369     format(1x,i5,2i3,2i2,2i4,' L.I VSOCC =',f12.6)
37	continue
38	continue


	DO 39 C1=1,NCH
	DO 39 C2=1,NCH
	DO 39 C2P=1,NCH
39	COUPF(C1,C2)=COUPF(C1,C2) 
     x                  + COUPI(C1,C2P) * VSOCC(C2P,C2,1)*  XB   !L.I
     x                  + COUPI(C1,C2P) * VSOCC(C2P,C2,2)*(-XQ)  !l.I

      ENDIF
!***********************************************************************
!		Now put COUPF into the main coupling arrays:
      NC=0
      DO 44 C1=1,NCH
	C2LAST=NCH
	if(REV.and..not.CPSO.and..not.cxwf) C2LAST=C1
      DO 44 C2=1,C2LAST
        IF(IC1.NE.PART(C1).OR.IC1.NE.PART(C2)) GO TO 44
      S = COUPF(C1,C2)
      if(S.ne.(0d0,0d0)) then
     	 PH12 = (0d0,1d0)**(LVAL(C2)-LVAL(C1))
           NC = NCLIST(C1,C2)+1
	   if(NC>MCLIST.and.REV.and..not.CPSO.and..not.cxwf) then
	      write(KO,*) '*** TOO MANY FORW COUPLINGS ',C1,C2
	      call check(NC,MCLIST,30)
	      endif
           if(NC<=MCLIST)then
             CLIST(C1,C2,NC) =  S*PH12
             NFLIST(C1,C2,NC) = KM
           endif
           NCLIST(C1,C2) =  NC
         if(REV.and..not.CPSO.and..not.cxwf.and.C1.ne.C2) then
           NC = NCLIST(C2,C1)+1
	   if(NC>MCLIST) then
	      write(KO,*) '*** TOO MANY REV COUPLINGS ',C1,C2
	     call check(NC,MCLIST,30)
	   endif
	   CLIST(C2,C1,NC) =  S*conjg(PH12)
	   NFLIST(C2,C1,NC) = KM
	   NCLIST(C2,C1) =  NC
         endif
         IF(C1.NE.C2) ICH = MAX(ICH,C2)
      IF(LISTCC.GE.1) WRITE(KO,43) C1,C2,KNP,KN,KM,IG,
     X IK,NC,LVAL(C1),LVAL(C2), S
43     FORMAT(' Folding:',4i4,2i5,2i3,2i5,F10.5/)
	endif
44    CONTINUE
45    CONTINUE
      call check(NC,MCLIST,30)
      GO TO 200
!***********************************************************************
C
C    TRANSFER COUPLING COEFFICIENTS   ZR: CLIST    FR: KCOEF
C    .......................................................
50    IF(IC1.EQ.IC2) GO TO 200
         IC1 = ICOM(CP,1)
         IC2 = ICOR(CP,1)
         C1FR = ICFROM .EQ. IC1
C    SO NOW NAME(1,IC1) IS LIKE D & NAME(1,IC2) LIKE P IN (D,P) REACTION
      T =      HP(IC1)/HP(IC2)
C                               T > 1
         MAXMV = 1
         DO 60 IN=1,NG(2)
         KN = GPT(2,IN)
         KNP= GPT(1,IN)
         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
         IM = 0
         DO 58 MVP=-LNP,LNP
         DO 58 MV =-LN,LN
           MM = MVP - MV
           IF(MM.LT.0) GO TO 58
           IM = IM + 1
58         CONTINUE
         MAXMV = MAX(MAXMV, IM)
60       ALOSS(IN) = 0.0
         ML1 = 1
         ML2 = 2
         IF(.NOT.LCALL) ML1 = 2
         IF(.NOT.MCALL) ML2 = 1
       DO 85 IL=ML1,ML2
C                        IL=1 => L COUPLINGS,  IL=2 => M-COUPLINGS
        IF(IL.EQ.1.AND.LISTCC.GT.0) WRITE(KO,61) 'L'
        IF(IL.EQ.2) THEN
           READ(14) VFOLD,VCORE,RINS,RINC,XJ
	   if(LISTCC>4) write(KO,601)  'VFOLD:',VFOLD
	   if(LISTCC>4) write(KO,601)  'VCORE:',VCORE
	   if(LISTCC>4) write(KO,*)  'RINS,RINC,XJ',RINS,RINC,XJ
601	   format(1x,a6,' potential:',/(1x,1p,6g14.5))
!           ICV = QQ
	   IC7 = MOD(ABS(QQ),2)
	   ICV = IC7
           NL0 = NL
	   if(IL.eq.2.and.LISTCC>4) write(KO,*)  'ICV,NLO',ICV,NLO
      	   LCLA = .FALSE.
      	   MCLA = .FALSE.
           MXMV1 = 1
           MXMVP1= 1
           MXLN1 = 1
           MXLNP1= 1
           MMX1 = 1
           LDMIN = MAXMUL
           LDMAX = -1
           IF(IL.EQ.2.AND.LISTCC.GT.0.AND.ML2.GE.2) WRITE(KO,61) 'M'
         ENDIF
61       FORMAT('0',A1,'-couplings:')
      REPEAT = .FALSE.
      if(.not.LOCAL.and.IL.eq.2) then
c				do pre-scan to find max size for MCG array:
      MKNL=0
      DO 62 C1=1,NCH
      DO 62 C2=1,NCH
        IF(IC1.NE.PART(C1).OR.IC2.NE.PART(C2)) GO TO 62
	MCL = .false.
	do 605 IN=1,NG(2)
        MCL = MCL .OR. .NOT.LTRANS(IN)
605    CONTINUE
	if(MCL) MKNL=MKNL+1
62	continue
      if(LISTCC>4) write(KO,*)  'ICV,NLO,MAXMV,MKNL',ICV,NLO,MAXMV,MKNL
      allocate (MCG(3,MAXMV,NG(2),MKNL))
      endif

      DO 78 C1=1,NCH
      DO 78 C2=1,NCH
        IF(IC1.NE.PART(C1).OR.IC2.NE.PART(C2)) GO TO 78
           IA1 = EXCIT(C1,1)
           IA2 = EXCIT(C2,1)
           IA1P= IA1
           IA2P= IA2
         IF(LOCAL) GO TO 64
            DO 63 IN=1,NG(2)
            DO 63 ILT=1,KLT
63          KCOEF(IN,ILT) = 0.
64       IF(COPY(1,IC1,IA1P).NE.0) IA1P = COPY(1,IC1,IA1P)
         IF(COPY(2,IC1,IA1 ).NE.0) IA1  = COPY(2,IC1,IA1 )
         IF(COPY(1,IC2,IA2P).NE.0) IA2P = COPY(1,IC2,IA2P)
         IF(COPY(2,IC2,IA2 ).NE.0) IA2  = COPY(2,IC2,IA2 )
      LTMIN = 9999
      LTMAX = -1
      DO 68 IN=1,NG(2)
         KN = GPT(2,IN)
         KNP= GPT(1,IN)
         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
      IF(ABS(AFRAC(ITC(IC2,IA2), ITC(IC1,IA1) ,2,KN)).LT.EPS .OR.
     &   ABS(AFRAC(ITC(IC1,IA1P),ITC(IC2,IA2P),1,KNP)).LT.EPS) GO TO 68
      LTMIN = MIN(LTMIN, MAX(ABS(LVAL(C2)-LN), ABS(LVAL(C1)-LNP) ))
      LTMAX = MAX(LTMAX, MIN(    LVAL(C2)+LN, LVAL(C1)+LNP))
68    CONTINUE
         NUMLT = LTMAX - LTMIN + 1
            IF(LOCAL) NUMLT = 1
       IF(NUMLT.LT.1) GO TO 78
C        CALL CHECK(NUMLT,KLT,10)
       IF(NUMLT.GT.KLT) WRITE(KO,*) 'KLT ERROR:',KLT,NUMLT
C     IF(LISTCC.GE.2) WRITE(KO,*) C1,IA1,IA1P,C2,IA2,IA2P,LTMIN,LTMAX,
C    #        NUMLT,         ,NG(2),LOCAL,CP
      LCL = .FALSE.
      MCL = .FALSE.
      DO 77 IN=1,NG(2)
      IF(IL.EQ.1 .AND. .NOT.LTRANS(IN)) GO TO 77
      IF(IL.EQ.2 .AND.      LTRANS(IN)) GO TO 77
         KNP = GPT(1,IN)
         KN  = GPT(2,IN)
         KM  = IN + LOCF
         LN = QNF(9,KN)
         SN = .5*QNF(10,KN)
         JN = 0.5*QNF(11,KN)
         LNP= QNF(9,KNP)
            IF(KIND.EQ.5) LNP = 0
         JNP= 0.5*QNF(11,KNP)
            IF(KIND.EQ.5) JNP = SN
       R6 = AFRAC(ITC(IC2,IA2),ITC(IC1,IA1), 2,KN)
       R6P= AFRAC(ITC(IC1,IA1P),ITC(IC2,IA2P),1,KNP)
         IF(KIND.EQ.5) R6P = 1.0
C           IF( QNF(12,KNP) .GT. 1) STOP 'DP-TRANSFER'
            IF( QNF(12,KNP) .GT. 1) CALL ABEND(32)
         DO 76 ILT=1,NUMLT
            LTOTAL = LTMIN + ILT - 1
            IF(LOCAL) LTOTAL = LVAL(C1)
      IF(ABS(R6) .LT. EPS.AND.LISTCC.LE.5) GO TO 76
      IF(KIND.GE.6 .AND. (ABS(R6P) .LT. EPS
     &    .OR. ABS(SN-.5*QNF(10,KNP)).GT.EPS).AND.LISTCC.LE.5) GO TO 76
      R3 = RACAH(JTARG(C1),JN,JTOTAL,JVAL(C2),JTARG(C2),JVAL(C1))
     &   * SQRT((TWO*JTARG(C2)+ONE)*(TWO*JVAL(C1)+ONE))
         IF(ABS(R3) .LT. EPS.AND.LISTCC.LE.5) GO TO 76
         FMAX = MIN(SN+JPROJ(C2),LNP+JPROJ(C1),LTOTAL+JVAL(C1))
         FMIN = MAX(ABS(SN-JPROJ(C2)),ABS(LNP-JPROJ(C1)),
     &              ABS(LTOTAL-JVAL(C1)) )
         R5D = 0.0
!      DO 72 F=FMIN,FMAX
      NNF=NINT(FMAX-FMIN)
      DO 72 IIF=0,NNF
      F=FMIN+IIF
      IF(LISTCC.GE.3) WRITE(KO,74) C1,C2,KNP,KN,KM,IN,LN,JN,LVAL(C1),
     &LVAL(C2), JVAL(C1),JVAL(C2),LTOTAL,F
      R50= WIG9J(LVAL(C2)+Z,JPROJ(C2),JVAL(C2),
     &           LN+Z,      SN,       JN,
     &           LTOTAL+Z  ,F        ,JVAL(C1))
     &   * SQRT((TWO*LTOTAL+ONE)*(TWO*F+ONE)*(TWO*JVAL(C2)+ONE)*
     &          (TWO*JN+ONE))
      R51 = RACAH(LNP+Z,SN,JPROJ(C1),JPROJ(C2),JNP,F)
      R52 = RACAH(LVAL(C1)+Z,LNP+Z,JVAL(C1),F,LTOTAL+Z,JPROJ(C1))
      R53 = SQRT((TWO*JNP+ONE)*(TWO*F+ONE)*(TWO*LTOTAL+ONE)*
     X           (TWO*JPROJ(C1)+ONE))
      R54 = PHASE(  NINT(SN + JPROJ(C2) - F))
72    R5D = R5D + R50 * R51 * R52 * R53 * R54
         IF(ABS (R5D) .LT. EPS.AND.LISTCC.LE.5) GO TO 76
      R8 =SQRT((2*LN+1)*(2*LVAL(C1)+1)*(2*LVAL(C2)+1)*ONE)*R4PI
      R9 = WIG3J(LN+Z,LVAL(C1)+Z,LVAL(C2)+Z, Z,Z,Z)
         IF(LOCAL .AND. ABS(R9) .LT. EPS.AND.LISTCC.LE.5) GO TO 76
      R0 = PHASE(LVAL(C2)-LN) / SQRT(2*LVAL(C1)+ONE)
         R8 = R8*R9*R0
         A = R6 * R6P * R3 * R5D
              IF(ABS (A).LT.EPS.AND.LISTCC.LE.5) GO TO 76
         S =     A * R8
         C6 = (0.,1.)**(LVAL(C1)-LVAL(C2))
         IF(LOCAL) THEN
      IF(REV.OR.IC1.EQ.ICTO) THEN
           NC = NCLIST(C1,C2)+1
	   if(NC>MCLIST) call check(NC,MCLIST,30)
	   CLIST(C1,C2,NC) =  T*S*CONJG(C6)
	   NFLIST(C1,C2,NC) = KM
	   NCLIST(C1,C2) =  NC
                 ICH = MAX(ICH,C2)
            ENDIF
      IF(REV.OR.IC2.EQ.ICTO) THEN
           NC = NCLIST(C2,C1)+1
	   if(NC>MCLIST) call check(NC,MCLIST,30)
	   CLIST(C2,C1,NC) =  (T*T)*S*C6
	   NFLIST(C2,C1,NC) = KM
	   NCLIST(C2,C1) =  NC
                 ICH = MAX(ICH,C1)
           ENDIF
         ELSE
            KCOEF(IN,ILT) = A
C
C ..................... THE I**(L2-L1) FACTOR ADDED TO KCOEF IN 'SOURCE'
         ENDIF
!      write(ko,*) QNF(3,KNP),EXCIT(C2,1)
        IF(LISTCC.GT.0) WRITE(KO,74) C1,C2,KNP,KN,KM,IN,LN,JN,LVAL(C1),
     &   LVAL(C2), JVAL(C1),JVAL(C2),LTOTAL,R3,R5D,R6,R6P,A,R8,S
74       FORMAT(/ 1X,7I3,F5.1,2I3,2F5.1,I4,4F9.5,1P,E19.12,0P,
     X            F9.5,2F9.3)
        LCL = LCL .OR. LTRANS(IN)
        MCL = MCL .OR. .NOT.LTRANS(IN)
76       CONTINUE
77    CONTINUE
      IF(LOCAL) GO TO 78
C
C    NON-LOCAL COUPLING FORM FACTORS (2-D ARRAYS) PUT ON FILE 12
C    ...........................................................
      IF(IL.EQ.1.AND.LCL)
     X CALL KERNEL(KCOEF,NLL,ALOSS,ALOSSM,RERR,C1,C2,IC1,IC2,
     &            REPEAT,LVAL,ICTO,ICFROM,REV,PART,DNL,XA,XB,
     &            NG(2),GPT(1,1),FILE,CHNO,CUTOFF,CP,KIND,NIB,
     &            NUMLT,LTMIN,XP,XQ,LTRANS,QNF,Q,HP(ICFROM),HP(ICTO),
     &            VREAL)
       IF(IL.EQ.2.AND.MCL) THEN
       if(.not.allocated(NM)) allocate(NM(MAXQRN,MFNL))
       CALL   KERNEM(VREAL,KCOEF,NLL,KIND.EQ.8,IREM,ICV,NM,LDMIN,LDMAX,
     &               C1,C2,IC1,IC2,REPEAT,  LVAL,ICTO,ICFROM,REV,PART,
     &               NG(1),NG(2),FPT,GPT,CHNO,CP,MCG,MAXMV,MKNL,NL0,QNF,
     &               KNL,NUMLT,LTMIN,LTRANS,IC7,
     X               MXLN1,MXLNP1,MMX1,MXMV1,MXMVP1)
       ENDIF
	LCLA = LCLA .or. LCL 
	MCLA = MCLA .or. MCL 
78    CONTINUE
      IF(IL.EQ.1.and.LCLA) THEN
           IF(.NOT.LOCAL.AND.(JTEST.LE.4.OR.NLPL.GT.0).AND..NOT.DRY.AND.
     X                  ALOSSM*RERR.GT.1E-5)
     X     WRITE(KO,80) CP,LOG10(MAX(ALOSSM,1D0)),
     X             (GPT(1,IN),GPT(2,IN),LOG10(MAX(ALOSS(IN),1D0)),
     X             ALOSS(IN)*RERR*100., IN=1,NG(2))
80    FORMAT(' For coupling #',I2,', acc. loss in kernel =',F5.1,
     X ' digits.',/,(10X,' Between bound states ',I3,' and ',I3,
     X '  loss =',F5.1,' D, so errors =',F10.5,' %'))
C
      ELSE if(IL.EQ.2.and.MCLA) then
      LRANG1 = MAX(LDMAX-LDMIN + 1,3)
      RIN = ONE/RINTP
      HF = HP(ICFROM)
      HT = HP(ICTO)
      HNL = (HF * MLT) / NLT
      RINTO = HT * MR
      NLCN = NLO/2 - NLC
      CENTRE = HF*MLT * NLCN
      write(48,*) 'CPAIR2: IC7 =',IC7
C
      CALL QERNEM(NLN,NLM,NLO,NG(1),NG(2),XA,XB,XP,XQ,FPT,CUTOFF,LTRANS,
     &            MCG,MAXMV,NL0,KNL,NM,MMX1,C1FR,
     &            IC7,ICV,IREM,NNU,RINTO,EPC,HNL,NLT,LRANG1,LDMIN,CP,
     &            VFOLD,VCORE,THM,RINS,RINC,NLL,QNF,NLOC,
     &            FORML,FORML(1,1,2),NLN,RIN,NLC,CENTRE,KIND.EQ.8,HF,HT,
     &            VREAL,XJ,LISTCC,CHNO,RERR,NCH,LVAL,NLPL,
     X            MXMV1,MXMVP1,MXLN1,MXLNP1,FORMC,FORMC(1,1,2),cxwf)
C
      ENDIF
85    CONTINUE
       if(allocated(NM)) deallocate (NM)   
      GO TO 200
C
C     GENERAL SPIN TRANSFERS OF BOTH PROJECTILE AND/OR TARGET
C     .......................................................
C
  90  CPLD(:,:,:) = .false.
      DO 99 C1=1,NCH
      DO 99 C2=1,NCH
      IF(ICTO.NE.PART(C1) .OR. ICFROM.NE.PART(C2)) GO TO 99
      REWIND 2
      L2 = 0
      DO 97 IX=NG(1),NG(2)
       IB = MATRIX(1,IX)
       IA = MATRIX(2,IX)
      IF(EXCIT(C1,1).NE.IB .OR. EXCIT(C2,1).NE.IA) GO TO 97
       LN = MATRIX(3,IX)
       KN = MATRIX(4,IX)
       SN = SPINTR(1,IX)
       JN = SPINTR(2,IX)
      IF(LISTCC.GE.4) WRITE(KO,921) IX,C1,C2,IB,IA,LN,SN,JN,KN
921   format(' IX=',i5,': C1,C2,IB,IA,LN,SN,JN,KN =',4i6,i4,2f5.1,i4)
C
	if(IP3.eq.1) then   ! CHEX2 form ?
      R6 = PHASE(NINT(JVAL(C2)+JN-JTARG(C1)-JTOTAL))
     X     * SQRT((TWO*JVAL(C1)+ONE)*(TWO*JVAL(C2)+ONE) *
     X            (TWO*JPROJ(C1)+ONE)*(TWO*JTARG(C1)+ONE) *
     X            (TWO*SN+ONE)*(TWO*JN+ONE)*(TWO*LVAL(C2)+ONE))
C
       R1 = PHASE(LVAL(C1)-LN) * WIG3J(LVAL(C1)+Z,LVAL(C2)+Z,LN+Z,Z,Z,Z)
     X      * SQRT((TWO*LVAL(C1)+ONE)*(TWO*LN+ONE))
       R2 = RACAH(JVAL(C1),JTARG(C1),JVAL(C2),JTARG(C2),JTOTAL,JN)
       R3 =   WIG9J(LVAL(C1)+Z,  LN+Z,    LVAL(C2)+Z,
     X              JPROJ(C1),   SN,      JPROJ(C2),
     X              JVAL(C1),    JN,      JVAL(C2))
       CH = R6 * R1 * R2 * R3
        IF(LISTCC.GE.1)WRITE(KO,922) R6,R1,R2,R3,real(CH)
922     format('  KIND=1: R6,R1,R2,R3,CH =',4f10.4,f12.5)
	else  ! IP3=0,2,3
C			New derivation of above expressions, 4th March 1998
C			The 1/sqrt(4pi) factor is included in INTER (frxp7.f)
C			Phase corrected, 27th August 1999
       R1 =  CLEB6(LVAL(C2)+Z,Z,LN+Z,Z,LVAL(C1)+Z,Z)
     X      * SQRT((TWO*LVAL(C2)+ONE)*(TWO*LN+ONE))
       R2 = WIG6J(JVAL(C2),JTARG(C2),JTOTAL,JTARG(C1),JVAL(C1),JN)
     X	    * PHASE(NINT(JVAL(C2)+JN+JTARG(C1)+JTOTAL))
!      R3 =   WIG9J(LVAL(C1)+Z,  LN+Z,    LVAL(C2)+Z,
!    X              JPROJ(C1),   SN,      JPROJ(C2),
!    X              JVAL(C1),    JN,      JVAL(C2))
       R3 =   WIG9J(LVAL(C2)+Z,  JPROJ(C2), JVAL(C2),
     X              LN+Z,        SN,        JN,
     X              LVAL(C1)+Z,  JPROJ(C1), JVAL(C1))
     X      * SQRT((TWO*JVAL(C1)+ONE)*(TWO*JVAL(C2)+ONE))
      if(IP3==2) then
!       R4 = sqrt((TWO*JPROJ(C1)+ONE)*(TWO*JTARG(C2)+ONE))   ! scaling by rme when SP,ST=0
       R4 = sqrt((TWO*JPROJ(C1)+ONE))   ! scaling by rme of projectile
       R4 = R4*((TWO*JTARG(C1)+ONE)*(TWO*JTARG(C2)+ONE))**0.25  ! scaling (symmetric!) so monopoles are physical
       else if(IP3==3) then
       R4 = sqrt((TWO*JTARG(C1)+ONE))   ! scaling by rme of projectile
       R4 = R4*((TWO*JTARG(C1)+ONE)*(TWO*JTARG(C2)+ONE))**0.25  ! scaling (symmetric!) so monopoles are physical
       else
       R4 = 1.0
       endif
       CH =  R1 * R2 * R3 * R4
        IF(LISTCC.GE.1) 
     x  WRITE(KO,923) IX,C1,C2,IB,IA,LN,SN,JN,KN,R1,R2,R3,R4,real(CH)
923   format(' IX=',i5,': C1,C2,IB,IA,LT,ST,JT,KN =',5i6,2f5.1,i4,
     X                  ' R1-4,CH =',5f8.4,f12.5)
	endif
            IF(ABS(CH).EQ.0.) GO TO 97
            IF(C1.NE.C2) ICH = MAX(ICH,C2)
      IF(KN.GT.0) THEN
C                     Local form factor:
           NC = NCLIST(C1,C2)+1
	   if(NC>MCLIST) call check(NC,MCLIST,30)
	   CLIST(C1,C2,NC) =  CH*(0.,1.)**(LVAL(C2)-LVAL(C1))
	   NFLIST(C1,C2,NC) = KN
	   NCLIST(C1,C2) =  NC
C	   if(LN<=4) then
C	   if(CPLD(C1,C2,LN)) WRITE(KO,924) C1,C2,IX,LN,SN,JN
C924  	     format(/' Chs ',2i4,' ALREADY coupled by',i3,':',i3,2f5.1)
	   CPLD(C1,C2,LN) = .true.
C	   endif
         IF(REV.and.IA/=IB) THEN
           NC = NCLIST(C2,C1)+1
 	   if(NC>MCLIST) call check(NC,MCLIST,30)
 	   CLIST(C2,C1,NC) =  CH*(0.,1.)**(LVAL(C1)-LVAL(C2))
 	   NFLIST(C2,C1,NC) = KN
 	   NCLIST(C2,C1) =  NC
C	   if(LN<=4) then
C	   if(CPLD(C2,C1,LN)) WRITE(KO,924) C2,C1,IX,LN,SN,JN
	   CPLD(C2,C1,LN) = .true.
C	   endif
 	 ENDIF
      ELSE
C            Non-local form factor:
       LREQ2 = - KN
       IF(LREQ2.NE.L2 .OR. IP3.EQ.-3) THEN
        DO 94 IL=1,LREQ2
 94     READ(66) FNC
       L2 = LREQ2
       ENDIF
        DO 95 I=1,NLL
        DO 95 J=1,NLO
95      FNC(I,J) = FNC(I,J) * CH
       IF(IP3.EQ.-3) THEN
C                         L-dependent additional factor:
      CALL FNLCC(FNC,NLL,NLO,HP(ICTO)*MR,DNL, LN,SN,JN,
     X           C1,C2,LVAL(C1),LVAL(C2),JPROJ(C1),JPROJ(C2),JTARG(C1),
     X           JTARG(C2),JVAL(C1),JVAL(C2),JTOTAL)
       ENDIF
C
       NL = NL + 1
       CALL CHECK(NL,MFNL,6)
       CHNO(NL,1) = C1
       IF(.NOT.REV) CHNO(NL,1) = - C1
       CHNO(NL,2) = C2
       CHNO(NL,3) = 13
       CHNO(NL,4) = NLL
       WRITE(12) FNC
      ENDIF
 97   CONTINUE
 99   CONTINUE
      GO TO 200
C
!@@
      
C
C    PHOTONUCLEAR COUPLING COEFFICIENTS 
C    .......................................................
100    IF(IC1.EQ.IC2) GO TO 200

         IC1 = ICOM(CP,2)
         IC2 = ICOR(CP,2)
      if(LISTCC>0) write(KO,*) ' Photonuclear couplings:',IC1,' <-',IC2
      
         IN2 = 2
         ZC = MASS(2+2,ICOR(CP,2))
         ZP = MASS(2+2,ICOM(CP,2)) - ZC
         XMC = MASS(  2,ICOR(CP,2))
!         XMP = MASS(  2,ICOM(CP,2)) - XMC
         XMP = MASS(  1,ICOR(CP,2)) 
         XRED=XMC*XMP/(XMC+XMP)
         
         if(LISTCC>3) 
     x        write(KO,*) 'ZC,ZP,XMC,XMP,XRED=',ZC,ZP,XMC,XMP,XRED

C    SO NOW NAME(1,IC1) IS LIKE D & NAME(1,IC2) LIKE P IN (D,P) REACTION
!!      T =      HP(IC1)/HP(IC2)
      T =      1.0


      DO 178 C1=1,NCH
      DO 178 C2=1,NCH
        IF(IC1.NE.PART(C1).OR.IC2.NE.PART(C2)) GO TO 178
        
           IA1 = EXCIT(C1,1)
           IA2 = EXCIT(C2,1)
         IF(COPY(2,IC1,IA1 ).NE.0) IA1  = COPY(2,IC1,IA1 )
         IF(COPY(2,IC2,IA2 ).NE.0) IA2  = COPY(2,IC2,IA2 )

!      IF(LISTCC.GE.2) WRITE(KO,*) C1,IA1,C2,IA2
 
      DO 177 IN=1,NK
         KN  = FPT(2,IN)
         IK  = FPT(3,IN)
         KM  = IN + LOCF
!         KMMI(1) = KM ;       KMMF(1) = KMMI(1)
         LN = QNF(9,KN)
         SN = .5*QNF(10,KN)
         JN = 0.5*QNF(11,KN)
       R1 = AFRAC(ITC(IC1,IA1),ITC(IC2,IA2), 2,KN)
!      write(KO,*) 'I2,I1,KN,R1 =',ITC(IC1,IA1),ITC(IC2,IA2),KN,R1
      IF(ABS(R1) .LT. EPS.AND.LISTCC.LE.5) GO TO 176
      IF(ABS(SN-JPROJ(C2)).gt. EPS .AND.LISTCC.LE.5) GO TO 176
      
      R2 = RACAH(JVAL(C2),JTOTAL,JN,JTARG(C1),JTARG(C2),IK+Z)

      R3 = RACAH(LVAL(C2)+Z,JVAL(C2),LN+Z,JN,JPROJ(C2),IK+Z)

      AE = 0.0	
      AM = 0.0
      if(mod(IREM,4)<=1.and.IK==LVAL(C1).and.abs(IK-JVAL(C1))<1e-5) then ! ELECTRIC
      R4 = CLEB6(LN+Z,Z,IK+Z,Z,LVAL(C2)+Z,Z)
      
      R5 = (-1)**nint(JVAL(C2)+JN-LVAL(C2))  * 
     X      sqrt((2*JVAL(C2)+1)*(2*LN+1)*(2*JN+1))  *      
     X      sqrt((2*JTARG(C1)+1)*2*(IK+1)*(2*IK+1)/IK / FINEC)
          
      R6 = XRED**IK * (ZP/XMP**IK + ZC/((-XMC)**IK))
                     
         AE =  R1 * R2 * R3 * R4 * R5 * R6
        IF(LISTCC.GT.0.and.abs(AE)>0.or.LISTCC>3) 
     &   WRITE(KO,174) 'E',IK,C1,C2,KN,KM,IN,LN,JN,
     &   LVAL(C1),LVAL(C2), JVAL(C1),JVAL(C2),IK,R1,R2,R3,R4,R5,R6,AE
      endif
         
!!    CHECK WHETHER THERE SHOULD BE JVAL(C1) dependence in the next block:
      if(mod(IREM,2)==0.and.IK-1==LVAL(C1)) then 			! MAGNETIC
      R4A = RACAH(LVAL(C2)+Z,1+Z,LN+Z,IK-1+Z,LVAL(C2)+Z,IK+Z)
      R4 = CLEB6(LN+Z,Z,IK-1+Z,Z,LVAL(C2)+Z,Z) 
      
      R5 = (-1)**nint(JVAL(C2)+JN-LVAL(C2))  * 
     X     2* sqrt((2*JVAL(C2)+1)*(2*LVAL(C2)+1)*(2*JN+1)*
     X             (2*(IK-1)+1)*IK*(2*IK+1)*LVAL(C2)*(LVAL(C2)+1))
!??  x       / ((IK+1)*(2*IK+1)) *      
     x       / ((IK+1)*sqrt(2.*IK+1.)) *      
     X      sqrt((2*JTARG(C1)+1)*2*(IK+1)*(2*IK+1)/IK / HBC) * mun 
     X       * K(IC1,IA1)
     
      R6 = pmass*(ZP/XMP * (XRED/XMP)**IK - ZC/XMC * ((-XRED)/XMC)**IK)     
                     
         AM = R1 * R2 * R3 * R4A * R4 * R5 * R6
        IF(LISTCC.GT.0.and.abs(AM)>0.or.LISTCC>3) 
     &   WRITE(KO,1741) 'M',IK,C1,C2,KN,KM,IN,LN,JN,
     &   LVAL(C1),LVAL(C2), JVAL(C1),JVAL(C2),IK,
     &   R1,R2,R3,R4A,R4,R5,R6,AM
     
      endif
      A = AE + AM
! 			Above theory is for coef(gamma) = -1
      A = A * HBC/K(IC1,IA1)
! 			Above theory is with k_gamma in cross section
      A = A * sqrt(K(IC1,IA1))

              IF(ABS (A).LT.EPS.AND.LISTCC.LE.5) GO TO 176
          
         C6 = (0.,1.)**(LVAL(C1)-LVAL(C2))

         IF(REV.OR.IC1.EQ.ICTO) THEN
           NC = NCLIST(C1,C2)+1
	   if(NC>MCLIST) call check(NC,MCLIST,30)
	   CLIST(C1,C2,NC) =  T*A*CONJG(C6) 
	   NFLIST(C1,C2,NC) = KM
	   NCLIST(C1,C2) =  NC
                 ICH = MAX(ICH,C2)
         if(LISTCC>0) WRITE(KO,1742) '>',IK,C1,C2,KN,KM,IN,LN,JN,
     &     LVAL(C1),LVAL(C2), JVAL(C1),JVAL(C2),IK,CLIST(C1,C2,NC)
           ENDIF
         IF(REV.OR.IC2.EQ.ICTO) THEN
           NC = NCLIST(C2,C1)+1
	   if(NC>MCLIST) call check(NC,MCLIST,30)
	   CLIST(C2,C1,NC) =  (T*T)*A*C6 * 0.5
	   NFLIST(C2,C1,NC) = KM
	   NCLIST(C2,C1) =  NC
                 ICH = MAX(ICH,C1)
         if(LISTCC>0) WRITE(KO,1742) '<',IK,C2,C1,KN,KM,IN,LN,JN,
     &     LVAL(C2),LVAL(C1), JVAL(C2),JVAL(C1),IK,CLIST(C2,C1,NC)
           ENDIF

174       FORMAT(/ 1X,a1,i1,1x,6I3,F5.1,2I3,2F5.1,I4,6F9.5,1P,E14.6)
1741       FORMAT(/ 1X,a1,i1,1x,6I3,F5.1,2I3,2F5.1,I4,7F9.5,1P,E14.6)
1742       FORMAT(/ 1X,a1,i2,1x,6I3,F5.1,2I3,2F5.1,I4,2G14.5)
1749       FORMAT(/ 1X,a,4i4,2f10.4,i3)
176       CONTINUE
177    CONTINUE


178    CONTINUE      
      
      
      
      
C
200   if(allocated(MCG)) deallocate(MCG)
C This next deallocation caused some x86_64 to crash (wierd)
!      deallocate(KCOEF,ALOSS,THM,FNC)
      RETURN
      END
*****ROTORNC************************************************************
      FUNCTION ROTORNC(Q,I,IP,K,KP)
      IMPLICIT NONE
        real*8,  intent(IN):: I,IP,K,KP
      	integer, intent(IN):: Q
        real*8:: ROTORNC,Z,WIG3J
      Z = 0.0
      ROTORNC = (-1)**NINT(K+IP) * SQRT( (2*I+1)*(2*IP+1) )
     &         * WIG3J(Q+Z,I,IP,Z,-K,KP)
      RETURN
      END
*****ROTOR**************************************************************
      FUNCTION ROTOR(J,I,L,IP,LP,Q,K,ROT)
      IMPLICIT REAL*8 (A-H,O-Z)
	real*8, intent(IN):: I,IP,K,J
      	integer, intent(IN):: L,LP,Q
      	logical, intent(IN):: ROT
      Z = 0.0
      PH = (-1)**NINT( (IP-I + ABS(IP-I))*0.5 )
      ROTOR  = SQRT((2.*L+1.) * (2.*LP+1.) )
     &       * (-1)**NINT(J - IP - L + LP)
     &       * RACAH(L+Z,LP+Z,I,IP,Q+Z,J)
     &       * CLEB6(L+Z,Z,LP+Z,Z,Q+Z,Z)
C
C     If 'ROT', then also include the additional Clebsh-Gordon
C               coefficient for the rotational model,
C               using the projection K quantum numbers.
C     Otherwise assume that this coefficient, or its equivalent,
C          is built into either the M(Ek) number for Coulomb
C              (but - after FRV - NOT the Mn(Ek) number)
C          or the 'reduced deformation length' RDEF(k) for nuclear
C          transitions.
C
       Rin = ROTOR
      IF(ROT) ROTOR = ROTOR * SQRT(2*IP+1) * PH
     X                      * CLEB6(IP,K,Q+Z,Z,I,K)
      RETURN
      END
*****TENSORS************************************************************
      FUNCTION TENSOR(TYPE,L,S1,J,S2,JT,LP,S1P,JP,S2P,LAM,K1,K2,K1P,K2P,
     &                ZP,ZT,ITYP)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8,intent(in):: S1,J,S2,JT,S1P,JP,S2P,K1,K2,K1P,K2P,ZP,ZT
	integer,intent(in):: TYPE,L,LP,LAM,ITYP
      REAL*8 J2,J2MIN,J2MAX
C
C   calculate matrix elements of spin tensors between partial waves:
C
C      < (L,S1)J, S2; JT  /  T(TYPE)  / (LP,S1P)JP, S2P; JT>
C           K1    K2                        K1p     K2p     projections
C
C      i.e. for single-particle wave functions with coupling order  3
C      and  for the main projectile/target channels
C
C      NB: NO i**L factors here, so results always real.
C
C    where  T(0) = ZP.ZT
C           T(1) = T(2) = 1
C           T(3) = r/1 . s/1(S1 S1)   = 2 l.S1  = S1 spin-orbit
C           T(4) = r/1 . s/1(S2 S2)   = 2 l.S2  = S2 spin-orbit
C           T(5) = r/2 . s/2(S1 S1)   = Tr for S1
C           T(6) = r/2 . s/2(S2 S2)   = Tr for S2
C           T(7) = r/2 . s/2(S1 S2)   = S1+S2 spin tensor
C           T(8) = S1.S2              = S1.S2 spin-spin
C           T(9)                         not used here
C           T(10)= T2(S1,R)              deformed projectile
C           T(11)=                    =  deformed target
C           T(12)=                    =  table couplings of projectile
C           T(13)=                    =  table couplings of target
C           T(14)=                    =  2nd-order projectile 
C           T(15)=                    =  2nd-order target 
C           T(16)=                    =  2nd-order projectile and target
C           T(30) = L(L+1)
C           T(others)                         not used here
C
      Z = 0.0
      TENSOR = 0.0
      IF(LAM.EQ.0.AND.ABS(S1-S1P)+ABS(S2-S2P).GT.0.01 .OR.TYPE.LT.0)
     X   RETURN
      GO TO (05,15,15,35,45,55,65,75,85,1,105,115,105,115,145,145,145), 
     X       		TYPE+1
	if(TYPE==30) go to 305
   1  RETURN
C
  05  IF(L.EQ.LP .AND. ABS(J-JP).LT.0.01) TENSOR = ZP*ZT
      RETURN
C
  15  IF(L.EQ.LP .AND. ABS(J-JP).LT.0.01) TENSOR = 1.0
      RETURN
C
  35  IF(.NOT.(L.EQ.LP .AND. ABS(J-JP).LT.0.01)) RETURN
      TENSOR = J * (J+1.) - L * (L+1.) - S1*(S1+1)
      RETURN
C
  45  IF(L.NE.LP) RETURN
      J2MIN = MAX(ABS(L-S2),ABS(S1-JT),ABS(LP-S2))
      J2MAX = MIN(    L+S2 ,    S1+JT ,    LP+S2 )
      T = 0.0
!      DO 50 J2=J2MIN,J2MAX
      NJ2=NINT(J2MAX-J2MIN)
      DO 50 IJ2=0,NJ2
      J2=J2MIN+IJ2
  50  T = T + (2*J2+1) * RACAH(S1,L+Z ,JT,S2,J,J2)
     &                 * RACAH(S1,LP+Z,JT,S2,JP,J2)
     &                 * (J2*(J2+1) - L*(L+1.) - S2*(S2+1))
      TENSOR = T * SQRT((2*JP+1)*(2*J+1)) * (-1)**NINT(J-JP+LP-L)
      RETURN
C
  55  IF(ABS(J-JP).GT.0.01) RETURN
      IF(S1.LT.1.) RETURN
      TENSOR = SQRT((2*L+1.)*(2*S1+1)) * (-1)**NINT(J-L-S1)
     &         * RACAH(L+Z,LP+Z,S1,S1,2+Z,J)
     &         *SQRT(2*(2*LP+1)/(3.*(2*L+1.)))*CLEB6(LP+Z,Z,2+Z,Z,L+Z,Z)
     &         * (3*S1*S1-S1*(S1+1))/(SQRT(6.)*CLEB6(S1,S1,2+Z,Z,S1,S1))
      RETURN
C
  65  J2MIN = MAX(ABS(L-S2),ABS(LP-S2),ABS(S1-JT))
      J2MAX = MIN(    L+S2 ,    LP+S2 ,    S1+JT )
      IF(S2.LT.1.) RETURN
      T = 0.0
!      DO 70 J2=J2MIN,J2MAX
      NJ2=NINT(J2MAX-J2MIN)
      DO 70 IJ2=0,NJ2
      J2=J2MIN+IJ2
  70  T = T + (2*J2+1) * RACAH(S1,L+Z, JT,S2,J, J2)
     &                 * RACAH(S1,LP+Z,JT,S2,JP,J2)
     &    * RACAH(L+Z,LP+Z,S2,S2,2+Z,J2)* (-1)**NINT(J2-L-S2)
      TENSOR = T * SQRT((2*L+1.)*(2*S2+1)*(2*J+1.)*(2*JP+1))
     &           * (-1)**NINT(J-JP+LP-L)
     &         *SQRT(2*(2*LP+1)/(3.*(2*L+1.)))*CLEB6(LP+Z,Z,2+Z,Z,L+Z,Z)
     &         * (3*S2*S2-S2*(S2+1))/(SQRT(6.)*CLEB6(S2,S2,2+Z,Z,S2,S2))
      RETURN
C
  75  SMIN = MAX(ABS(L-JT),ABS(S1-S2))
      SMAX = MIN(    L+JT ,    S1+S2 )
      T = 0.0
!      DO 80 S=SMIN,SMAX
      NS=NINT(SMAX-SMIN)
      DO 80 IS=0,NS
      S=SMIN+IS
      T1 = SQRT((2*J+1)*(2*S+1)) * RACAH(L+Z, S1,JT,S2,J,S)
        SPMIN= MAX(ABS(LP-JT),ABS(S1-S2),ABS(S-2))
        SPMAX= MIN(    LP+JT ,    S1+S2 ,    S+2 )
!      DO 80 SP=SPMIN,SPMAX
      NSP=NINT(SPMAX-SPMIN)
      DO 80 ISP=0,NSP
      SP=SPMIN+ISP
      T2 = SQRT((2*JP+1)*(2*SP+1)) * RACAH(LP+Z,S1,JT,S2, JP,SP)
      T3 = SQRT((2*L+1.)*(2*S+1.)) * (-1)**NINT(JT-L-SP)
     &      * RACAH(L+Z,LP+Z,S,SP,2+Z,JT)
      T4 = SQRT(2*(2*LP+1)/(3.*(2*L+1.))) *CLEB6(LP+Z,Z,2+Z,Z,L+Z,Z)
      T5 = SQRT((2*SP+1)*(2*2+1)*(2*S1+1)*(2*S2+1))
     &      * WIG9J(S,   SP,   2+Z,
     &              S1,  S1,   1+Z,
     &              S2,  S2,   1+Z)
     &      * SQRT(S1*(S1+1) * S2*(S2+1))
      T = T + T1 * T2 * T3 * T4 * T5
 80    CONTINUE
      TENSOR = T
      RETURN
C
C       					S1.S2 spin.spin
  85  IF(L.NE.LP) RETURN
      SMIN = MAX(ABS(L-JT),ABS(S1-S2))
      SMAX = MIN(    L+JT ,    S1+S2 )
      T = 0.0
!      DO 90 S=SMIN,SMAX
      NS=NINT(SMAX-SMIN)
      DO 90 IS=0,NS
      S=SMIN+IS
      T1 = SQRT((2*J+1)*(2*JP+1))*(2*S+1)
     X       	 * RACAH(L+Z, S1,JT,S2,J,S)
     X           * RACAH(LP+Z,S1,JT,S2,JP,S)
      T2 = (S*(S+1) - S1*(S1+1) - S2*(S2+1))*0.5
      T = T + T1 * T2 
 90   CONTINUE
      TENSOR = T
      RETURN
 95   RETURN
C                                                 PROJECTILE ROTOR or TA
 105  CH = ROTOR(J,S1,L,S1P,LP,LAM,K1,TYPE.LE.11)
      IF(ITYP.EQ.0) CH = CH * ZT
      IF(ABS(J-JP).gt.1e-5) CH = 0.0
      TENSOR = CH
      RETURN
C                                                 TARGET ROTOR or TABLE
 115  CH = 0.0
      J2MIN = MAX(ABS(S2-L),ABS(JT-S1),ABS(S2P-LP))
      J2MAX = MIN(    S2+L,     JT+S1 ,   S2P+LP )
      CH = 0.0
      T1 = (-1)**NINT(J-JP-L+LP)* SQRT((2*J+1)*(2*JP+1))
!      DO 120 J2=J2MIN,J2MAX
      NJ2=NINT(J2MAX-J2MIN)
      DO 120 IJ2=0,NJ2
      J2=J2MIN+IJ2
         T2 = (2*J2+1) * RACAH(S1,L+Z,JT,S2,J,J2)
     &                * RACAH(S1,LP+Z,JT,S2P,JP,J2)
         C6 = ROTOR(J2,S2,L,S2P,LP,LAM,K2,TYPE.LE.11)
 120   CH = CH + T1 * T2 * C6
      IF(ITYP.EQ.0) CH = CH * ZP
      TENSOR = CH
      RETURN
C                                            2nd ORDER PROJECTILE AND/OR TARGET
 145  WRITE(KO,*) 'Normal 2nd-order algebra should not be called!'
      WRITE(KO,*) TYPE,L,S1,J,S2,JT,LP,S1P,JP,S2P,LAM,K1,K2,K1P,K2P,ITYP
	stop

 305  IF(L.EQ.LP .AND. ABS(J-JP).LT.0.01) TENSOR = L*(L+1)
      RETURN
C
      END
      FUNCTION TENSLS(TYPE,L,S1,S,S2,JT,LP,S1P,SP,S2P,LAM,K1,K2,K1P,K2P,
     &                ZP,ZT,ITYP)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 JT,K1,K2,K1P,K2P, J2,J2MIN,J2MAX,J1,J1MIN,J1MAX
      INTEGER TYPE
C
C   calculate matrix elements of spin tensors between partial waves:
C
C      < (L,(S1,S2)S; JT  /  T(TYPE)  / (LP,(S1P,S2P)SP; JT>
C            K1 K2                           K1p K2p     projections
C
C      i.e. for single-particle wave functions with coupling order 1
C
C    where  T(0) = ZP*ZT
C           T(1) = T(2) = 1
C           T(3) = r/1 . s/1(S1 S1)   = 2 l.S1  = S1 spin-orbit
C           T(4) = r/1 . s/1(S2 S2)   = 2 l.S2  = S2 spin-orbit
C           T(5) = r/2 . s/2(S1 S1)   = Tr for S1
C           T(6) = r/2 . s/2(S2 S2)   = Tr for S2
C           T(7) = r/2 . s/2(S1 S2)   = S1+S2 spin tensor
C           T(8) = S1.S2              = S1.S2 spin-spin
C           T(10)= T2(S1,R)           =  deformed projectile
C           T(11)=                    =  deformed target
C           T(12)=                    =  table couplings of projectile
C           T(13)=                    =  table couplings of target
C           T(14)=                    =  2nd-order projectile 
C           T(15)=                    =  2nd-order target 
C           T(16)=                    =  2nd-order projectile and target
C           T(20)                         JLS-dependent potentials
C           T(21)                         JLS-dependent potentials
C
      Z = 0.0
      TENSLS = 0.0
      IF(LAM.EQ.0.AND.ABS(S1-S1P)+ABS(S2-S2P).GT.0.01 .OR. TYPE.LT.0)
     X   RETURN
      GO TO (05,15,15,35,45,55,65,75,85,1,105,115,105,115,
     X	     145,145,145,1,1,1,205,205),TYPE+1
   1  RETURN
C
  05  IF(L.EQ.LP .AND. ABS(S-SP).LT.0.01) TENSLS = ZP*ZT
      RETURN
C
  15  IF(L.EQ.LP .AND. ABS(S-SP).LT.0.01) TENSLS = 1.0
      RETURN
C
  35  IF(L.NE.LP) RETURN
      J1MIN = MAX(ABS(L-S1),ABS(S2-JT),ABS(LP-S1))
      J1MAX = MIN(    L+S1 ,    S2+JT ,    LP+S1 )
      T = 0.0
!      DO 40 J1=J1MIN,J1MAX
      NJ1=NINT(J1MAX-J1MIN)
      DO 40 IJ1=0,NJ1
      J1=J1MIN+IJ1
  40  T = T + (2*J1+1) * RACAH(L+Z, S1,JT,S2,J1,S)
     &                 * RACAH(LP+Z,S1,JT,S2,J1,SP)
     &                 * (J1*(J1+1) - L*(L+1.) - S1*(S1+1))
      TENSLS = T * SQRT((2*SP+1)*(2*S+1))
      RETURN
  45  IF(L.NE.LP) RETURN
      J2MIN = MAX(ABS(L-S2),ABS(S1-JT),ABS(LP-S2))
      J2MAX = MIN(    L+S2 ,    S1+JT ,    LP+S2 )
      T = 0.0
!      DO 50 J2=J2MIN,J2MAX
      NJ2=NINT(J2MAX-J2MIN)
      DO 50 IJ2=0,NJ2
      J2=J2MIN+IJ2
  50  T = T + (2*J2+1) * RACAH(L+Z ,S2,JT,S1,J2,S)
     &                 * RACAH(LP+Z,S2,JT,S1,J2,SP)
     &                 * (J2*(J2+1) - L*(L+1.) - S2*(S2+1))
      TENSLS = T * SQRT((2*SP+1)*(2*S+1)) * (-1)**NINT(S-SP)
      RETURN
C
  55  J1MIN = MAX(ABS(L-S1),ABS(S2-JT),ABS(LP-S1))
      J1MAX = MIN(    L+S1 ,    S2+JT ,    LP+S1 )
      T = 0.0
!      DO 60 J1=J1MIN,J1MAX
      NJ1=NINT(J1MAX-J1MIN)
      DO 60 IJ1=0,NJ1
      J1=J1MIN+IJ1
  60  T = T + (2*J1+1) * RACAH(L+Z, S1,JT,S2,J1,S)
     &                 * RACAH(LP+Z,S1,JT,S2,J1,SP)
     &         * RACAH(L+Z,LP+Z,S1,S1,2+Z,J1)* (-1)**NINT(J1-L-S)
      TENSLS = T * SQRT((2*SP+1)*(2*S+1)*(2*L+1.)*(2*S1+1))
     &         *SQRT(2*(2*LP+1)/(3.*(2*L+1.)))*CLEB6(LP+Z,Z,2+Z,Z,L+Z,Z)
     &         * (3*S1*S1-S1*(S1+1))/(SQRT(6.)*CLEB6(S1,S1,2+Z,Z,S1,S1))
      RETURN
C
  65  J2MIN = MAX(ABS(L-S2),ABS(LP-S2),ABS(S1-JT))
      J2MAX = MIN(    L+S2 ,    LP+S2 ,    S1+JT )
      T = 0.0
!      DO 70 J2=J2MIN,J2MAX
      NJ2=NINT(J2MAX-J2MIN)
      DO 70 IJ2=0,NJ2
      J2=J2MIN+IJ2
  70  T = T + (2*J2+1) * RACAH(L+Z, S2,JT,S1,J2,S )
     &                 * RACAH(LP+Z,S1,JT,S1,J2,SP)
     &    * RACAH(L+Z,LP+Z,S2,S2,2+Z,J1)* (-1)**NINT(J2-L-S2)
      TENSLS = T * SQRT((2*L+1.)*(2*S2+1)*(2*S+1)*(2*SP+1))
     &           * (-1)**NINT(S-SP)
     &         *SQRT(2*(2*LP+1)/(3.*(2*L+1.)))*CLEB6(LP+Z,Z,2+Z,Z,L+Z,Z)
     &         * (3*S2*S2-S2*(S2+1))/(SQRT(6.)*CLEB6(S2,S2,2+Z,Z,S2,S2))
      RETURN
C
  75  T = 0.0
      T3 = SQRT((2*L+1.)*(2*S+1.)) * (-1)**NINT(JT-L-SP)
     &      * RACAH(L+Z,LP+Z,S,SP,2+Z,JT)
      T4 = SQRT(2*(2*LP+1)/(3.*(2*L+1.))) *CLEB6(LP+Z,Z,2+Z,Z,L+Z,Z)
      T5 = SQRT((2*SP+1)*(2*2+1)*(2*S1+1)*(2*S2+1))
     &      * WIG9J(S,   SP,   2+Z,
     &              S1,  S1,   1+Z,
     &              S2,  S2,   1+Z)
     &      * SQRT(S1*(S1+1) * S2*(S2+1))
      TENSLS = T3 * T4 * T5
      RETURN
C
C       					S1.S2 spin.spin
 85   TENSLS = (S*(S+1) - S1*(S1+1) - S2*(S2+1)) * 0.5
      RETURN
C                                                 PROJECTILE ROTOR or TA
 105  IF(ABS(S2-S2P)+ABS(K2-K2P)+ABS(K1-K1P).GT.0.01) RETURN
      J1MIN = MAX(ABS(S1-L),ABS(JT-S2),ABS(S1P-LP))
      J1MAX = MIN(    S1+L,     JT+S2 ,    S1P+LP )
      CH = 0.0
      T1=(-1)**NINT(S-SP + S-S1-S2 + SP-S1P-S2P)* SQRT((2*S+1)*(2*SP+1))
!      DO 110 J1=J1MIN,J1MAX
      NJ1=NINT(J1MAX-J1MIN)
      DO 110 IJ1=0,NJ1
      J1=J1MIN+IJ1
         T2 = (2*J1+1) *  RACAH(L+Z, S1, JT,S2,J1,S )
     &                 *  RACAH(LP+Z,S1, JT,S1,J1,SP)
         C6 = ROTOR(J1,S1,L,S1P,LP,LAM,K1,TYPE.LE.11)
 110   CH = CH + T1 * T2 * C6
      IF(ITYP.EQ.0) CH = CH * ZT
      TENSLS = CH
      RETURN
C                                                 TARGET ROTOR or TABLE
 115  IF(ABS(S1-S1P)+ABS(K2-K2P)+ABS(K1-K1P).GT.0.01) RETURN
      J2MIN = MAX(ABS(S2-L),ABS(JT-S1),ABS(S2P-LP))
      J2MAX = MIN(    S2+L,     JT+S1 ,    S2P+LP )
      CH = 0.0
      T1 = (-1)**NINT(S-SP)* SQRT((2*S+1)*(2*SP+1))
!      DO 120 J2=J2MIN,J2MAX
      NJ2=NINT(J2MAX-J2MIN)
      DO 120 IJ2=0,NJ2
      J2=J2MIN+IJ2
         T2 = (2*J2+1) *  RACAH(L+Z, S2, JT,S1,J2,S )
     &                 *  RACAH(LP+Z,S2, JT,S2,J2,SP)
         C6 = ROTOR(J2,S2,L,S2P,LP,LAM,K2,TYPE.LE.11)
 120   CH = CH + T1 * T2 * C6
      IF(ITYP.EQ.0) CH = CH * ZP
      TENSLS = CH
      RETURN
C                                            2nd ORDER PROJECTILE AND/OR TARGET
 145  CH = 0.0
      TENSOR = 1.234567
	WRITE(KO,*) 'TENSLS 2nd-order algebra not yet implemented!'
	stop

C
 205  IF(ABS(S1-S1P)+ABS(S2-S2P)+ABS(K1-K1P)+ABS(K2-K2P) +
     &   ABS(S - SP) .GT.0.01) RETURN
C                                     specific nucleon-nucleon potential
      IF(ABS(S).LT.0.01) THEN
C                                     singlet
      IF(L.NE.LP) RETURN
      IF(L.EQ.0 .AND. LAM.EQ.1 .OR.
     &   L.EQ.1 .AND. LAM.EQ.5 .OR.
     &   L.EQ.2 .AND. LAM.EQ.11)   TENSLS = 1.0
      ELSE
C                                     Triplet
      IJ = NINT(JT)
      IF(L.EQ.LP) THEN
C                                             diagonal
        IF(L.EQ.0 .AND. LAM.EQ.2 .OR.
     &     L.EQ.1 .AND. IJ.EQ.0 .AND. LAM.EQ.6 .OR.
     &     L.EQ.1 .AND. IJ.EQ.1 .AND. LAM.EQ.7 .OR.
     &     L.EQ.1 .AND. IJ.EQ.2 .AND. LAM.EQ.8 .OR.
     &     L.EQ.2 .AND. IJ.EQ.1 .AND. LAM.EQ.4 .OR.
     &     L.EQ.2 .AND. IJ.EQ.2 .AND. LAM.EQ.12.OR.
     &     L.EQ.3 .AND. IJ.EQ.2 .AND. LAM.EQ.10)    TENSLS = 1.0
      ELSE
C                                             off-diagonal
        IF(L.EQ.0 .AND. LP.EQ.2 .AND. IJ.EQ.1 .AND. LAM.EQ.3 .OR.
     &     L.EQ.2 .AND. LP.EQ.0 .AND. IJ.EQ.1 .AND. LAM.EQ.3 .OR.
     &     L.EQ.1 .AND. LP.EQ.3 .AND. IJ.EQ.2 .AND. LAM.EQ.9 .OR.
     &     L.EQ.3 .AND. LP.EQ.1 .AND. IJ.EQ.2 .AND. LAM.EQ.9) TENSLS=1.
        ENDIF
      ENDIF
      RETURN
C
      END
      FUNCTION TENDEF(TYPE,L,S1,J,LP,JP,LAM,K,ZP,ZT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 J,JP,K
      INTEGER TYPE
C
C   calculate matrix elements of deformation between partial waves:
C
C      < (L,S1)J K /  T(TYPE, deformation=LAM  / (LP,S1)JP K >
C
C      i.e. for single-particle wave functions with coupling order 2
C
C    where  T(0) = ZP*ZT
C           T(1) = T(2) = 1
C           T(3) = r/1 . s/1(S1 S1)   = 2 l.S1  = S1 spin-orbit
C
      Z = 0.0
      TENDEF = 0.0
      IF(TYPE.GT.3 .OR. TYPE.LT.0) RETURN
      GO TO (05,05,05,35), TYPE+1
C
  05  IF(J.GT.K .OR. JP.GT.K) RETURN
      T      = (-1)**NINT(J-JP+LAM)
     &         * SQRT((2*J+1)*(2*LAM+1)/(2*JP+1))
     &         * CLEB6(J,S1,LAM+Z,Z,JP,S1)
     &         * CLEB6(J,K, LAM+Z,Z,JP,K )
      IF(TYPE.EQ.0) T = T * ZP*ZT
      TENDEF = T
      RETURN
C
  35  IF(.NOT.(L.EQ.LP .AND. ABS(J-JP).LT.0.01)) RETURN
C              the spin orbit potential is not itself deformed (yet|)
      TENDEF = J * (J+1.) - L * (L+1.) - S1*(S1+1)
      RETURN
C
      END
      FUNCTION TENS0(TYPE,L,S1,J,S2,LP,S1P,JP,S2P,LAM,ZP,ZT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 J,JP
      INTEGER TYPE
C
C   calculate matrix elements of spin tensors between partial waves:
C
C      < (L,S1)J  /  T(TYPE)  / (LP,S1P)JP>
C
C      i.e. for single-particle wave functions with coupling order 0
C
C    where  T(0) = ZP*ZT
C           T(1) = T(2) = 1
C           T(3) = r/1 . s/1(S1 S1)   = 2 l.S1  = S1 spin-orbit
C
      Z = 0.0
      TENS0 = 0.0
      IF(ABS(S1-S1P)+ABS(S2-S2P).GT.0.01 .OR. LAM.NE.0) RETURN
      IF(.NOT.(L.EQ.LP .AND. ABS(J-JP).LT.0.01)) RETURN
      GO TO (05,15,15,35), TYPE+1
         RETURN
C
  05  TENS0 = ZP*ZT
      RETURN
C
  15  TENS0 = 1.0
      RETURN
C
  35  TENS0 = J * (J+1.) - L * (L+1.) - S1*(S1+1)
      RETURN
      END
************************************************************************

      FUNCTION TENSOR2(TYPE,L,S1,J,S2,JT,LP,S1P,JP,S2P,LAM,LL1,LL2,
     &                 K1,K2,K1P,K2P,ZP,ZT,ITYP)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8,intent(in):: S1,J,S2,JT,S1P,JP,S2P,K1,K2,K1P,K2P,ZP,ZT
	integer,intent(in):: TYPE,L,LP,LAM,ITYP,LL1,LL2
      REAL*8 J2,J2MIN,J2MAX,I2,I2MIN,I2MAX
C
C   calculate matrix elements of spin tensors between partial waves:
C
C      < (L,S1)J, S2; JT  /  T(TYPE)  / (LP,S1P)JP, S2P; JT>
C           K1    K2                        K1p     K2p     projections
C
C		of  multipole LAM=LL1*LL2
C
C      NB: NO i**L factors here, so results always real.
C
C    where  T(14)=                    =  2nd-order projectile 
C           T(15)=                    =  2nd-order target 
C           T(16)=                    =  2nd-order projectile and target
C           T(others)                         not used here
C
      Z = 0.0
      TENSOR2 = 0.0
      IF(TYPE<14.or.TYPE>16) RETURN
C
C                                            2nd ORDER PROJECTILE AND/OR TARGET
 145  CH = 0.0
      TENSOR2 = CLEB6(LL1+Z,Z,LL2+Z,Z,LAM+Z,Z) * 1.23456

      END

      SUBROUTINE CCSET(JTOTAL,PARITY,ETOTAL,JTMIN,KINTL,
     X        NEX,NCHAN,GIVEXS,QVAL,ENEX,PEL,EXL,LMAX,JEX,ECM,ECMC,
     X        LVAL,JVAL,PART,EXCIT,COPY,JPROJ,JTARG,CUTL,CUTR,CUTOFF,HP,
     X        RMASS,INCOME,BLOCKD,LUSED,MAL1,LJMAX,
     X        MINTL,INITL,ITC,ITL,IEX,NCH,NCHPART,CHBASE,
     X        K,ETA,CHL,NAME,BAND,N,ISOCEN,RTURN,SMALLS,NSMALL,DROPPED)
	use io
	use parameters
	use trace
	use searchpar, only: final
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CHL(LMAX1,MXPEX,2)
      INTEGER PEL,EXL,NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     X    ITC(MXP,MXX),NCHPART(MXP),CHBASE(MXP),SMALLS(MXPEX)
      INTEGER LVAL(MAXCH),PARITY,C,PART(MAXCH,3),EXCIT(MAXCH,3)
     X       ,INITL(MAXCH),CUTOFF,DROPPED

      LOGICAL BLOCKD(MAXCH),GIVEXS(MXP)
C
C   STATEMENT FUNCTIONS
C   -------------------
      LOGICAL FRAC,FAIL3
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      REAL*8 RMASS(MXP),HP(MXP),JEX(6,MXP,MXX)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1)
C
C    INCOMING ENERGIES
C    -----------------
      REAL*8 JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH),
     X       JTOTAL,JAP,JTMIN,JN,ECM(MAXCH,3),LJMAX,
     X  K(MXP,MXX),ETA(MXP,MXX),ETOTAL,ECMC(MXP,MXX),RMK
C
C
      CHARACTER*8 NAME(2,MXP+1)
      CHARACTER*1 PSIGN(3)
      CHARACTER*1 BLANK,ANI,INCOME(MAXCH)
C
C
      DATA PSIGN / '-','?','+' /, BLANK,ANI / ' ','I' /
C
C     PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-NINT(X)).GT.1D-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      Z = 0.0
      CUTOFF = MAX(abs(CUTL)*INT(JTOTAL), CUTR/HP(1)) + 1
C     if(CUTR.lt.0.) CUTOFF = max(CUTOFF,int((RTURN+CUTR)/HP(1)))
      if(CUTR.lt.0.) CUTOFF = max(CUTOFF,
     *                  int(10.**(CUTR/(JTOTAL+1))*RTURN/HP(1)) )
      C = 0
      MINTL = 0
      DROPPED = 0
      IEX = 0
!      if(MAXPW<LMAX) WRITE(KO,*) ' CCSET: Max partial wave maxpw =',maxpw
      DO 240 IC=1,NCHAN
      NA = NEX(IC)
	NCHPART(IC) = 0
	CHBASE(IC) = C
      DO 240 IA=1,NA
      IF(JTMIN.LT.Z .AND.
     X JTOTAL.LE.-JTMIN-.1 .AND. (IC.NE.PEL .OR. IA.NE.EXL)) GO TO 240

       if(SMALLS(ITC(IC,IA)).gt.NSMALL) then
	 DROPPED = DROPPED+1
	 go to 240
	 endif
      LPAR = PARITY *  SIGN(1, BAND(1,IC,IA)*BAND(2,IC,IA) )
      ITP = 2.*JEX(1,IC,IA) + 1.5
      JN = JEX(1,IC,IA) + JEX(2,IC,IA)
      LAPF1 = MAX(0, NINT(JTOTAL-JN)) + 1
      LAPL1 = MIN(MAXPW,LMAX, NINT(JTOTAL+JN)) + 1
      DO 230 LAP1=LAPF1,LAPL1
         LAP = LAP1 - 1
         IF(LPAR.NE.(-1)**LAP) GO TO 230
C-----------------------------------------------------------
       IF(COPY(2,IC,IA,2).EQ.IC) THEN
C HAVE EXCHANGED NUCLEI.
C                        Omit EVEN L if spin zero nuclei.
          IF(ABS(JEX(1,IC,IA)).LT.1E-5 .AND.
     X           (-1)**LAP .LT. 0) GO TO 230
        ENDIF
C-----------------------------------------------------------
      DO 220 IJAP=1,ITP
         JAP = LAP + IJAP - JEX(1,IC,IA) - 1.
         IF(JAP.LT.0.0 .OR. FRAC(2.*JN) ) GO TO 220
         T = LAP
         IF(FAIL3(JAP,JEX(1,IC,IA),T)) GO TO 220
         IF(FAIL3(JAP,JEX(2,IC,IA),JTOTAL)) GO TO 220
      C = C + 1
      IF(C.GT.MAXCH) GO TO 220
      ECM(C,1) = ECMC(IC,IA)
      LVAL(C)= LAP
      JVAL(C)= JAP
      PART(C,1)= IC
      EXCIT(C,1)=IA
      NCHPART(IC) = NCHPART(IC) + 1
        DO 210 IN=1,2
        PART(C,IN+1)= IC
        EXCIT(C,IN+1)=IA
        IN2 = IN
         IF(COPY(IN,IC,EXCIT(C,IN+1),2).NE.0) THEN
             PART (C,IN+1)=COPY(IN,IC,EXCIT(C,IN+1),2)
             IN2 = 3 - IN
           ENDIF
         IF(COPY(IN2,PART(C,IN+1),EXCIT(C,IN+1),1).NE.0)
     X       EXCIT(C,IN+1)=COPY(IN2,PART(C,IN+1),EXCIT(C,IN+1),1)
  210   CONTINUE
      JPROJ(C)= JEX(1,IC,IA)
      JTARG(C)= JEX(2,IC,IA)
!@@
      IF(RMASS(IC).lt.1e-5) then   ! coefficient of second derivative
          ECM(C,2) = - HBC**2 / ECMC(IC,IA)
      ELSE 
          ECM(C,2) = -1.0 / (FMSCAL * RMASS(IC))
      ENDIF
!@@
      ECM(C,3) = HP(IC)
C     CUTOFF(C) = MAX( ETA(IC,IA)*LAP/(K(IC,IA)*25.)/ECM(C,3), abs(CUTL)*LAP)
C     CUTOFF(C) = MAX(abs(CUTL)*LAP, CUTR/HP(IC)) + 1

      INCOME(C) = BLANK
      BLOCKD(C) = ABS(BAND(1,IC,IA)).EQ.9 .OR. ABS(BAND(2,IC,IA)).EQ.9
      IF(LUSED.LT.LAP) LUSED = LAP
      IF(GIVEXS(IC).AND.LJMAX.LT.ABS(JTOTAL-LAP)) LJMAX=ABS(JTOTAL-LAP)
      IF(.NOT.(IC.EQ.PEL.AND.IA.EQ.EXL)) GO TO 220
       IF(MINTL.LT.KINTL.or.KINTL<0) then
         MINTL = MINTL+1
         INITL(MINTL) = C
         INCOME(C) = ANI
        CALL CHECK(C,MAXICH,25)
	endif
  220 CONTINUE
  230 CONTINUE
  240 IF(ITC(IC,IA).LE.ITL) IEX = C
      CALL CHECK(C,MAXCH,-5)
C         IF(C.LT.MINTL) STOP
          IF(C.LT.MINTL) CALL ABEND(8)
      NCH = MIN(C,MAXCH)
          IEX = MIN(IEX,NCH)
      CALL CHECK(LUSED+1,MIN(LMAX1,MAL1),11)
C
      IF(NCH*MINTL.EQ.0) RETURN
      IF(SMATL.GE.3) WRITE(KO,1290)
      R1 = (CUTOFF-1)*HP(1)
      IF(CHANS+LISTCC+SMATL.GE.3.and.final)  then
           WRITE(KO,1270) JTOTAL,PSIGN(PARITY+2),
     X                              NCH    ,IEX,R1,RTURN
 1270 FORMAT( /,1X,118('#'),/,' #',116X,'#',/
     X           ,' #',' Total SPIN and PARITY =',F8.1,1X,A1,',',I7,
     X' channels,',I5,' in 1st block.  Rmin & Coul turning =',
     X   2F7.1,' fm.  #'/ ' #',116X,'#',/1X,118('#')/)
	   else if(CHANS+LISTCC+SMATL.GE.1.and.final)  then
           WRITE(KO,1271) JTOTAL,PSIGN(PARITY+2),NCH,IEX,R1,RTURN
 1271 FORMAT(/' Total SPIN, PARITY =',F8.1,1X,A1,',',I7,
     X' chs,',I7,' cc.  Rmin & Coul turning =',2f7.1)
     	   endif
      IF(NCH.EQ.0) RETURN
C
      IF(CHANS.GT.0.and.final) WRITE(KO,1280)
 1280 FORMAT(//'    C Projectl Target   #   EX. ',
     X             '  (L  Proj) J  + Targ = Jtotal     E-cm      Re K',
     X '     Re Eta     RM*K            CH')
C
      M = N-1
      DO 260 C=1,NCH
         IC = PART(C,1)
         IA = EXCIT(C,1)
         L1 = LVAL(C) + 1
        IF(ISOCEN.eq.1) L1 = JTOTAL + 0.1 + 1
        IF(ISOCEN.eq.2) L1 = LVAL(INITL(1)) + 1

         RMK = (M-1)*ECM(C,3) * ABS(K(IC,IA))
C        R0 = (CUTOFF(C) -1) * ECM(C,3)
         R1 = K(IC,IA)
         R2 = ETA(IC,IA)
         R3 = RMK
      IF(CHANS.GT.0.and.final)  then
	 if(abs(aimag(CHL(L1,ITC(IC,IA),1)*(0.0,0.5)))<1e3) then
           WRITE(KO,1290) C,NAME(1,IC),NAME(2,IC),IA,
     X      INCOME(C),mod(LVAL(C),1000),JPROJ(C),mod(JVAL(C),1000d0),
     X      JTARG(C),JTOTAL,ECM(C,1),R1,R2,R3,
     X      CHL(L1,ITC(IC,IA),1)*(0.0,0.5),
     X      (PART(C,IN+1),mod(EXCIT(C,IN+1),100),IN=1,2)
	else
           WRITE(KO,1291) C,NAME(1,IC),NAME(2,IC),IA,
     X      INCOME(C),mod(LVAL(C),1000),JPROJ(C),mod(JVAL(C),1000d0),
     X      JTARG(C),JTOTAL,ECM(C,1),R1,R2,R3,
     X      CHL(L1,ITC(IC,IA),1)*(0.0,0.5)
	endif
	endif
  260    CONTINUE
 1290 FORMAT(I5,2(1X,A8),' #',I5,': ',A1,I4,F4.1,F7.1,f5.1,f8.1,
     X  f12.5,2F10.5,F10.4,2F10.5,3I2,I4)
 1291 FORMAT(I5,2(1X,A8),' #',I5,': ',A1,I4,F4.1,F7.1,f5.1,f8.1,
     X  f12.5,2F10.5,F10.4,1p,2e11.3)
C
      IF(MINTL.EQ.0)  WRITE(KO,1310)
 1310 FORMAT
     X(//' ***** INCOMING PROJECTILE STATE NOT INCLUDED HERE *****'//)
      IF(MINTL.EQ.0) RETURN
      IF(CHANS.GT.0.and.CHANS.le.9) CHANS = CHANS - 1
      RETURN
      END
      SUBROUTINE NUMCC(NCH,IEX,MINTL,MEL,JTOTAL,PARITY,JTMIN,KINTL,
     X        NEX,NCHAN,GIVEXS,PEL,EXL,LMAX,JEX,COPY,NFUSCH,
     X        RMASS,MAL1,ITC,ITL,BAND,SMALLS,NSMALL,CHPRES)
	use parameters
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER PEL,EXL,NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     X    ITC(MXP,MXX),PARITY,C,SMALLS(MXPEX),CHPRES(MXPEX),C1
      LOGICAL FRAC,FAIL3,GIVEXS(MXP)
      REAL*8 RMASS(MXP),JEX(6,MXP,MXX)
      REAL*8 JTOTAL,JAP,JTMIN,JN
C
      FRAC(X) = ABS(X-NINT(X)).GT.1D-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      Z = 0.0
      C = 0
      MINTL = 0
      MEL = 0
      IEX = 0
      NFUSCH = 0
      DO 240 IC=1,NCHAN
      NA = NEX(IC)
      DO 240 IA=1,NA
      IF(JTMIN.LT.Z .AND.
     X JTOTAL.LE.-JTMIN-.1 .AND. (IC.NE.PEL .OR. IA.NE.EXL)) GO TO 240
       if(SMALLS(ITC(IC,IA)).gt.NSMALL) go to 240
      LPAR = PARITY *  SIGN(1, BAND(1,IC,IA)*BAND(2,IC,IA) )
      ITP = 2.*JEX(1,IC,IA) + 1.5
      JN = JEX(1,IC,IA) + JEX(2,IC,IA)
      LAPF1 = MAX(0, NINT(JTOTAL-JN)) + 1
      LAPL1 = MIN(LMAX, NINT(JTOTAL+JN)) + 1
      C1 = C
      DO 230 LAP1=LAPF1,LAPL1
         LAP = LAP1 - 1
         IF(LPAR.NE.(-1)**LAP) GO TO 230
C-----------------------------------------------------------
       IF(COPY(2,IC,IA,2).EQ.IC) THEN
C HAVE EXCHANGED NUCLEI.
C                        Omit EVEN L if spin zero nuclei.
          IF(ABS(JEX(1,IC,IA)).LT.1E-5 .AND.
     X           (-1)**LAP .LT. 0) GO TO 230
        ENDIF
C-----------------------------------------------------------
      DO 220 IJAP=1,ITP
         JAP = LAP + IJAP - JEX(1,IC,IA) - 1.
         IF(JAP.LT.0.0 .OR. FRAC(2.*JN) ) GO TO 220
         T = LAP
         IF(FAIL3(JAP,JEX(1,IC,IA),T)) GO TO 220
         IF(FAIL3(JAP,JEX(2,IC,IA),JTOTAL)) GO TO 220
      C = C + 1
      if(ITC(IC,IA)<=NFUS) NFUSCH=C
      if(IC.EQ.PEL.AND.IA.EQ.EXL.and.(MINTL.LT.KINTL.or.KINTL<0))then
         MINTL = MINTL+1
	 MEL = max(MEL,C)
	 endif
  220 CONTINUE
  230 CONTINUE
	CHPRES(ITC(IC,IA)) = C-C1
  240 IF(ITC(IC,IA).LE.ITL) IEX = C
      NCH = C
      RETURN
      END
      SUBROUTINE XCH(PSI,N,NCH,NICH,SMAT,JTOTAL,LVAL,JVAL,JPROJ,JTARG,
     &               PART,EXCIT,COPY,MXP,MXX,MAXN,SMATS,DROP,EXONLY,
     & 		     EXCH,MAXCH,EL)
      use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 PSI(MAXN,NICH),PSIN(N,NICH)
      COMPLEX*16 SMAT(NCH),SNEW(NCH)
      REAL*8 JTOTAL,JVAL(NCH),JPROJ(NCH),JTARG(NCH),EXCH(MAXCH,MAXCH)
      INTEGER LVAL(NCH),PART(NCH),EXCIT(NCH),COPY(2,MXP,MXX,2),SMATS,
     &        C1,C2,EL
      LOGICAL DROP(NCH),EXONLY
C
      Z = 0.0
       DROP(:) = .FALSE.
       EXCH(:,:) = 0.
	SNEW(:) = SMAT(:)
	PSIN(:,:) = PSI(1:N,1:NICH)
      DO 50 C1=1,NCH
       IF(ABS(JPROJ(C1)-JTARG(C1)).GE.1E-5) THEN
C                                               DIFFERENT STATES
         PX = 1.0
      ELSE IF(JPROJ(C1)-AINT(JPROJ(C1)).GT.1E-5) THEN
C                                   FERMIONS: EQUAL HALF-INTEGRAL SPINS
         PX = - 1.0
      ELSE
C                                   BOSONS :  EQUAL      INTEGRAL SPINS
         PX = 1.0
      ENDIF
      DO 50 C2=1,NCH
C
C                      COPY C2 ONTO C1 IF C2 IS AN EXCHANGE COPY OF C1
C
      IF(LVAL(C1)/=LVAL(C2)) GO TO 50
      IF( COPY(1,PART(C2),EXCIT(C2),2)/=PART(C1) .and.  ! different partition
     &   (COPY(2,PART(C2),EXCIT(C2),2)/=PART(C1).or.
     &    PART(C1)/=PART(C2)) )  GO TO 50		! same partition
      IF(EXCIT(C2)/=EXCIT(C1)) GO TO 50
      DIAG = 0.
      if(C2==EL) DIAG = 1.0
C
      T = PX * (-1)**NINT(JTOTAL + LVAL(C2) - JVAL(C1) - JVAL(C2))
     &    * SQRT((2*JVAL(C1)+1.) * (2.*JVAL(C2)+1.)) *
     &    RACAH(JPROJ(C1),LVAL(C2)+Z,JTOTAL,JTARG(C1),JVAL(C1),JVAL(C2))
     &    * (-1)**LVAL(C2)
      EXCH(C1,C2) = T
      if(EXONLY) then
        IF(SMATS.GE.5) WRITE(KO,15) C1,PX,T,C2
15      FORMAT(' Channel',I3,' has',F4.0,' Exchange',
     &         ' by adding',F9.5,' times channel',I3)
      else
      SNEW(C1) = SNEW(C1) + T * (SMAT(C2) - DIAG)
        IF(SMATS.GE.5 .AND. ABS(SMAT(C2)).GT.-1E-6.or..true.)
     &                WRITE(KO,20) C1,SNEW(C1),PX,T,C2,SMAT(C2)-DIAG
20      FORMAT(' S-mat(',I3,') =',2F10.6,' after',F4.0,' Exchange',
     &         ' by adding',F9.5,' times S-mat(',I3,')-D =',2F10.6)
      IF(C1.GT.NICH .OR. C2.GT.NICH) GO TO 40
      PSIN(:,C1) = PSIN(:,C1) + T * PSI(:,C2)
      endif
40    DROP(C2) = .TRUE.
50    CONTINUE
	if(EXONLY) return
C
      SMAT(:) = SNEW(:)
      PSI(1:N,1:NICH) = PSIN(1:N,1:NICH)
      DO 70 C2=1,NCH
       IF(.NOT.DROP(C2).or.COPY(2,PART(C2),EXCIT(C2),2)==PART(C2)) 
     x         GO TO 70
       SMAT(C2) = 0.0
       IF(C2.GT.NICH) GO TO 70
         PSI(:,C2) = 0.0
70      CONTINUE
      RETURN
      END
      
      SUBROUTINE CFRACT(N,G0,G,A,C,D)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C
C*******************************************************************
C
C  dfract converts polynomial g to the corresponding continued
C         fraction, in stieltjes form with coefficients a
C
C   g(z) = g0+g1*z+g2*z**2+g3*z**3+...+gn*z**n
C
C   gn(z)= cn(z)/dn(z)
C        = g0/1+ a1*z/1+a2*z/1+a3*z/1+.../1+an*z
C
C  data:
C   n     order of g, even or odd       input
C   g0    constant term                 input
C   g     vector g(k), k=1,n            input
C   a     vector a(k), k=1,n            output
C   c     output: numerator   polynomials, convergents n and n-1
C   d     output: denominator polynomials, convergents n and n-1
C         caller provides space for a,c,d
C   c and d arrays contain nt=2*((n+1)/2) elements each
C   storage conventions for c and d:
C     c0=g0 and d0=1.0
C     coefficients in sequence, with m=n-1 and k=(n+1)/2 :
C       cn1,cm1,cn2,cm2,...cni,cmi,...  i=1,k
C       dn1,dm1,dn2,dm2,...dni,dmi,...  i=1,k
C     note that dnk=0.0 if n is odd
C   algorithm: rk nesbet, 82.10.27, to be published
C
C*******************************************************************
C
      DIMENSION G(N),A(N),C(N+1),D(N+1)
      REAL*8 ZERO,ONE
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
      IF(N.LE.0) RETURN
      NT=2*((N+1)/2)
      DO 10  I=1,NT
      C(I)=ZERO
  10  D(I)=ZERO
      CT=ZERO
      DT=ONE
      DN= G0
      DO 500 K=1,N
      DD=DN
      DN=G(K)
      DO 50 I=1,K,2
      IF((I+1)/2.GE.K) GOTO 50
      DN=DN+G(K-(I+1)/2)*D(I)
  50  CONTINUE
      A(K)=-DN/DD
C plant ak=0.0 and return if sequence Truncates
      IF(A(K).EQ.0D0) RETURN
      DO 100 I=1,K,2
      CI=C(I)
      C(I)=CI+A(K)*CT
      CT=C(I+1)
      C(I+1)=CI
      DI=D(I)
      D(I)=DI+A(K)*DT
      DT=D(I+1)
      D(I+1)=DI
  100 CONTINUE
      CT=G0
      DT=ONE
  500 CONTINUE
C
      RETURN
      END
      FUNCTION CEVALP(B,X,INUM,K)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION B(1)
C
C     evaluate polynomial 'b' at x
C       b(n*k+1) used up to n = inum,
C     evalp = b(1) + b(k+1)*x + b(2k+1)*x**2 + ... + b(inum.k+1)*x**inum
C
      F = 0.0
      IF(INUM.LT.0) GO TO 20
      DO 10 I=1,INUM
      N = INUM+1 - I
10    F = (B(N*K+1) + F) * X
      F = F + B(1)
20    CEVALP = F
      RETURN
      END
      SUBROUTINE PAD(PADE,SMAT,NCH,ITNL,SPAD,MAXCH,MAXIT2)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER PADE,C
      COMPLEX*16 SMAT(NCH),SPAD(MAXCH,MAXIT2,MAXIT2),C6
C
C               THIS IS THE PADE APPROXIMANT BY THE EPSILON ALGORITHM
C               -----------------------------------------------------
      DO 440 C=1,NCH
      SPAD(C,MIN(ITNL,MAXIT2),1) = 0.0
  440 SPAD(C,MIN(ITNL,MAXIT2),2) = SMAT(C)
      IF(PADE.EQ.0) RETURN
C
      IT1 = ITNL + 1
      DO 460 J=3,IT1
      IT = ITNL + 2 - J
         T = 0.0
         R1 = 0.
         R2 = 0.
      DO 450 C=1,NCH
      T = T + ABS(SPAD(C,IT+1,J-1) - SPAD(C,IT,J-1))**2
      R1= R1+  ABS(SPAD(C,IT+1,J-1))
  450 R2= R2+  ABS(                   SPAD(C,IT,J-1))
      R0 = 0.0
      IF(T.NE.0 .AND. R1*R2.NE.0) R0 = 1.0/T
      DO 460 C=1,NCH
        C6 = 0.
         IF(R0.NE.0.) C6 = CONJG(SPAD(C,IT+1,J-1) - SPAD(C,IT,J-1)) * R0
  460 SPAD(C,IT,J) = SPAD(C,IT+1,J-2) + C6
      DO 470 C=1,NCH
  470 SMAT(C) = SPAD(C,IT1-IT1/2*2+1, IT1/2*2)
      RETURN
      END
      SUBROUTINE FNLCC(FNC,NLL,NLO,H,DNL, LTR,PTR,TTR,
     X                 C1,C2,L1,L2,JP1,JP2,JT1,JT2,J1,J2,JTOTAL)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 FNC(NLL,NLO)
      REAL*8 DNL(NLO),JP1,JP2,JT1,JT2,J1,J2,JTOTAL
      INTEGER C1,C2
C                    Any additional L-dependent factors:
      F = 1.0
      IF(C1.NE.C2) F = 0.0
      DO 20 I=1,NLL
       R = (I-1)*H
      DO 20 J=1,NLO
        RP = R + DNL(J)
        IF(RP.LE.0.0) GO TO 20
        FNC(I,J) = FNC(I,J)* L1 * F
20      CONTINUE
      RETURN
      END
*****CRISS**************************************************************
      SUBROUTINE CRISS(IC,IA,JEX,NSA,NJA,PEL,EXL,JIX,
     &         MAL,K,ETA,RMASS,CSIG,KOORDS,MMXCH,MAXF,ENLAB,LEN,
     &         THMIN,THMAX,THINC,EXCH,LAMPL,NEARFA,KQMAX0,SIGR,XSTRAC,
     &         LJMAX,NLJ,NJ,JBORD,JUMP,ITX,ITEL,IEXCH,PP,PPK,CDCC,IP,
     &         IBIN,PENEX,PI,HEADNG,INFAM,OUTFAM,LCROSS,LFAM,
     & 		LXSEC,BSIGN,ECMI,A1,A2,ECMF,A3,A4,LEG,MAXPLM,
     &          DSPINS,DMULTIES,NMULTIES,DLEVEL,DNAME)
	use io
	use parameters
	use searchpar
	use searchdata
      IMPLICIT REAL*8(A-H,O-Z)
      ! PARAMETER(KQMAX1=7)
      REAL*8 JEX(2),JIX(2),RMASS(MXP),MBP,MJAP,MJBP,JTOTAL,MAA,PENEX,
     & 	     JVAL(MMXCH),MODSQ,JBORD(NJ+1),LJMAX,MA,MAP,MT,MB,MJA,MJB,
     &       JT,JLAST,AI(4),JPREV,XSEC(KQMAX0,KQMAX0),TENS(5),TENHJ(3)
      INTEGER PEL,EXL,PART(MMXCH),EXCIT(MMXCH),LVAL(MMXCH),C,EL,XSTRAC
     &       ,JUMP(7,3),PP,CDCC,BSIGN,INFAM,OUTFAM,KAD(mds),KANG(mds),
     &        IANGLS(mdl,mds),DMULTIES(NMULTIES),DLEVEL(NMULTIES),iSgn
      REAL*8 RR(MMXCH,MAXF),ANGLES(mds*mdl),KYY,CM_LAB(mds*mdl)
      COMPLEX*16 O,C6,C7,C9,C0,FCOUL,C6D,FAM(MAXF),SMAT(MMXCH),FCOUL2,
     &        AMP,FA(KQMAX0,KQMAX0,KQMAX0,2-KQMAX0:KQMAX0),CI,FAR,NEAR
      REAL*8 CSIG(LMAX1,MXPEX),K(MXP,MXX),ETA(MXP,MXX),FUSL(1+NFUS1)
      REAL*8 PL(LMAX1,MAXPLM+1,2),PLEG(KQMAX0),SIGG(KQMAX0),
     &       RCOEF(KQMAX0,NMULTIES),DSPINS(NMULTIES),GXSECS(NMULTIES),
     &       TGXSECS(NMULTIES)
      COMPLEX*16,allocatable:: AMPL(:),YSIG(:,:)
      LOGICAL EXCH,DESCR,IFO,EXTRA
      CHARACTER*12 SIDES(3),PPK*10
      CHARACTER*80 HEADNG,HDGQ
      CHARACTER*8 DNAME
      DATA  SIDES / 'BOTH SIDES','FAR SIDE','NEAR SIDE' /,
     X      DESCR / .false. /
      MODSQ(O) = CONJG(O) * O
!	 write(51,*) 'CRISS:',IC,IA,XSTRAC,LCROSS,SIGR
!	call flush(51)
      Z = 0.0
      CI= (0D0,1D0)
      R1 = 4*PI
      NSB = 2*JEX(1) + 1.1
      NJB = 2*JEX(2) + 1.1
      XCOEF = 10.0/(NSA * NJA)
         MAM = NSA * NJA * NSB * NJB
!	 write(48,*) 'CRISS:',NSA,NJA,NSB,NJB,NSA*NJA*NSB*NJB,MAM
!	 call flush(48)
      KQMAX1 = KQMAX0
      KQMAX2  = 1

      if(mod(PP,4).eq.0) KQMAX1  = MIN( nint(2*JIX(1)+1) , KQMAX1 )
      if(PP.eq.1) KQMAX1  = MIN( nint(2*JIX(2)+1) , KQMAX1 )
      if(PP.eq.2) KQMAX1  = MIN( nint(2*JEX(1)+1) , KQMAX1 )
      if(PP.eq.3) KQMAX1  = MIN( nint(2*JEX(2)+1) , KQMAX1 )
      if(PP.eq.4) KQMAX2  = MIN( nint(2*JEX(1)+1) , KQMAX1 )
      IKYY = 0
      if(PP==4.and.KQMAX2>1) IKYY = 1
      IF(IEXCH.NE.0) WRITE(KO,*)  ' Elastic Projectile-Target EXCHANGE '
     x         ,IEXCH
      
      KA = 0; NANGLE=1
	      IANGLS(:,:) = -1
      do id=1,datasets
      if(IC==data_ic(id).and.IA==data_ia(id)) then
        KA=KA+1
        KAD(KA) = id
           do ip=1,datalen(id)
	    if(data_type(id)<=2.and.
     x         (abs(ENLAB-data_energies(ip,id))<1e-5
     x	       .or. data_energies(ip,id)<0.)) then
             TH = datangles(ip,id)
	     xcm_lab = 1.
	     if(data_lab(id)) then  ! convert lab datangles to cm ANGLES
	       call lab2cm(TH,A1,A2,A3,A4,ECMI,ECMF,xcm_lab,PI)
	     endif
             ANGLES(NANGLE) = TH
             CM_LAB(NANGLE) = xcm_lab
	     IANGLS(ip,id) = NANGLE
             NANGLE=NANGLE+1
	    endif
           enddo
      endif
      enddo
         NANGLE=NANGLE-1
      
      HDGQ = '"'//HEADNG(1:lnbl(HEADNG))//'"'
      KQ1PR = 0
      IF(XSTRAC.NE.0) KQ1PR = MIN(KQMAX1,ABS(XSTRAC))
         MAM = NSA * NJA * NSB * NJB
         MAL1 = MAL + 1
!	 write(48,*) 'CRISS:',NSA,NJA,NSB,NJB,NSA*NJA*NSB*NJB,MAM
!	 write(51,*) 'CRISS:',NSA,NJA,NSB,NJB,NSA*NJA*NSB*NJB,MAM
!	 call flush(51)
         LASTA = MAM * MAL1
           NBJ = 0
           DO 10 I=1,NJ
C         IF(XSTRAC.GE.1) WRITE(KO,*) ' FOR',I,', JUMPS =',(JUMP(I,C)
C    #            ,C=1,3)       ,' FROM',JBORD(I)
10         IF(JUMP(I,1).GT.1) NBJ = MAX(NBJ,JUMP(I,2)+1)
            IF(NBJ.GT.0 .AND. NBJ.LT.5) THEN
               WRITE(KO,*) '0RERUN WITH FINER J-STEPS !!!!!'
               GO TO 301
               ENDIF
            LB = LASTA - MAM * (1+NLJ)
	    allocate (AMPL(LASTA + MAM*NLJ*NBJ))
         CALL CHECK(MAM,MAXF,-20)
         IF(MAM.GT.MAXF) GO TO 300
      IT1=NINT(JIX(1)+JIX(2)+JEX(1)+JEX(2))
      MIT = 2*IT1 + 1
	allocate (YSIG(0:MAL,MIT))
!@@
      IF(RMASS(IC).lt.1e-5) then
         C0 = SQRT(          RMASS(PEL)*amu/ (HBC*K(PEL,EXL)) )
!@@@     C0 = SQRT( K(IC,IA)*RMASS(PEL)*amu/ (HBC*K(PEL,EXL)) )
     &        / K(PEL,EXL)
      ELSE IF(RMASS(PEL).lt.1e-5) then
         C0 = 1./SQRT(          RMASS(IC)*amu/ (HBC*K(IC,IA)) )
      ELSE
         C0 = SQRT(K(IC,IA)/RMASS(IC) / (K(PEL,EXL)/RMASS(PEL)))
     &        / K(PEL,EXL)
      ENDIF
!     see justification in frxx3.f 
!
!     C0 = SQRT(K(IC,IA)/RMASS(IC) / (K(PEL,EXL)/RMASS(PEL)))
!    &          / K(PEL,EXL)
!@@
C
      REWIND 10
       DO 20 IAM=1,MAM
       DO 20 LP=0,MAL
20      AMPL(LP*MAM+IAM) = 0.
C
      NEXT = 0
      JLAST = - 1.0
      JPREV = - 1.0
      DO 23 JBL=1,NJ
C     IF(XSTRAC.GE.2) WRITE(KO,213) JBL,JUMP(JBL,3)
      IF(JUMP(JBL,1).LT.1) GO TO 23
      DO 211 IT=1,MIT
         MP = -IT1 + IT - 1
         MPA= ABS(MP)
      DO 211 LP=MPA,MAL
         T = CSIG(LP+1,ITX) - CSIG(1,ITX)
         IF(JUMP(JBL,1).GT.1) T = -T
         C6 = EXP(CI*T)
         R7 = YLMC(LP,MP)
211   YSIG(LP,IT) = C6 * R7
         DO 19 IBJ=NEXT+1,NBJ
         DO 19 LJ=1,NLJ
         DO 19 IAM=1,MAM
           IB = LB + IAM + LJ*MAM + IBJ*NLJ*MAM
19       AMPL(IB) = 0.0
C
      DO 223 ICC=1,JUMP(JBL,3)
C
      READ(10,END=223) JTOTAL,NCH,MINTL,IPARI,JUMPER
	if(NCH>MMXCH) then
	  WRITE(KO,*) ' Num. channels =',NCH,' > ',MMXCH,'!!!'
	  write(51,*) ' Num. channels =',NCH,' > ',MMXCH,'!!!'
	  call flush(51)
	  stop
	endif
      READ(10) (LVAL(C),JVAL(C),PART(C),EXCIT(C),C=1,NCH)
         IF(ABS(JTOTAL-JLAST).GT.1E-5 .OR.NEXT.EQ.0) NEXT = NEXT + 1
         JLAST = JTOTAL
         IF(JUMP(JBL,1).EQ.1) JPREV = JTOTAL
C
      DO 222 JIN=1,MINTL
      READ(10) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1)
      LA = LVAL(EL)
      R3 = SQRT((2*LA+1)/(4*PI))
!      IF(IEXCH.NE.0) R3 = R3 * (1 - (-1)**(LA + (IEXCH+1)/2))
      T = CSIG(LA+1,ITEL) - CSIG(1,ITEL)
         C7 = C0 * EXP(CI*T) * R3
C     IF(XSTRAC.GE.2) WRITE(KO,213) JBL,ICC,JTOTAL,NCH,EL
C213   FORMAT(/' CC#',2I4,' IS',    F5.1,2I4)
         IAM = 0
      DO 222 IMA = 1,NSA
      MA = -JIX(1) + IMA-1
      MAA = 0
      R2 = CLEB6(LA+Z,MAA,JIX(1),MA,JVAL(EL),MAA+MA)
      DO 222 IJA = 1,NJA
      MJA = -JIX(2) + IJA-1
      MT = MAA + MA + MJA
      R2P= CLEB6(JVAL(EL),MAA+MA,JIX(2),MJA,JTOTAL,MT)
      R22P = R2 * R2P
      DO 222 IMB = 1,NSB
      MB = -JEX(1) + IMB-1
      DO 222 IJB = 1,NJB
      MJB = -JEX(2) + IJB-1
         IAM = IAM + 1
      MP = NINT(MT - MB - MJB)
         IT = MP + IT1 + 1
C     IF(XSTRAC.GE.5) WRITE(KO,218) MA,MJA,MB,MJB,ICC,LA,MP,R0,R1,R2,R2P
C    X                                ,R3
C218   FORMAT(/1X,4F5.1,2I3,      I5,5F9.5)
C
      DO 220 C=1,NCH
      IF(PART(C) .NE. IC) GO TO 220
      IF(EXCIT(C) .NE. IA) GO TO 220
      IF(JVAL(C) .LT. ABS(MT-MJB)) GO TO 220
      IF(JTOTAL .LT. ABS(MT)) GO TO 220
      LP = LVAL(C)
      IF(LP .LT. ABS(MP)) GO TO 220
C
      IF(JIN.NE.1) R145 = RR(C,IAM)
      IF(JIN.NE.1) GO TO 219
      R4 = CLEB6(LP+Z,MP+Z,JEX(1),MB,JVAL(C),MT-MJB)
      R5 = CLEB6(JVAL(C),MT-MJB,JEX(2),MJB,JTOTAL,MT)
      R145 = R1 * R4 * R5
      RR(C,IAM) = R145
219   C6 = YSIG(LP,IT)
      C9 = SMAT(C) 
      IF(C.EQ. EL) C9 = (SMAT(C) - 1) 
      C9 = C9 * (0.0, -0.5)
C
      AMP = (R22P * R145) * (C7*C9*C6)
C
         IF(JUMP(JBL,1).GT.1) THEN
            LJ = NINT(LP - JTOTAL + LJMAX) + 1
            IB = LB + IAM + LJ*MAM + NEXT*NLJ*MAM
         ELSE
            IB = LP*MAM+IAM
         ENDIF
      AMPL(IB) = AMPL(IB) + AMP
C
C     IF(LAMPL.EQ.IC.AND.(XSTRAC.GE.6.OR.XSTRAC.GE.5.AND.LP.EQ.0))
C    #WRITE(KO,368) C,PART(C),EXCIT(C),JVAL(C),LP,R145,C6,C7,C9,
C    X          AMP,AMPL(IB)
C
C 368    FORMAT(' TRACE',3I3,F5.1,I3,7F8.4,2F9.5,2F10.5)
220   CONTINUE
222   CONTINUE
223   CONTINUE
C
      IF(JUMP(JBL,1).EQ.1) NEXT = 0
      IF(JUMP(JBL,1).EQ.1) GO TO 23
            IF(NEXT.GT.NBJ .OR. NEXT.LT.4
     &       .OR.JBL.LT.NJ.AND.ABS(JLAST-JBORD(JBL+1)).GT.0.1) THEN
               WRITE(KO,*) JBL,ICC,JTOTAL,JLAST,NEXT,NBJ,NJ
               WRITE(KO,*) 'JUMP:',JUMP
               WRITE(KO,*) 'JBORD:',JBORD,';',JBORD(JBL+1)
C              STOP 'BLOCKING AMPLITUDES'
               CALL ABEND(8)
            ENDIF
C
!         DO 100 JT=JBORD(JBL),JLAST
         NJT=NINT(JLAST-JBORD(JBL))
         DO 100 IJT=0,NJT
         JT=JBORD(JBL)+IJT
               IF(ABS(JT-JPREV).LT.0.1) GO TO 100
C              IF(JBL.LT.NJ .AND. JT.GE.JLAST-0.1) GO TO 100
            T = (JT - JBORD(JBL)) / JUMP(JBL,1)
c
c   james damper
            damp=1.d0
            if (jbl.eq.nj.and.jbl.gt.2) then
               damp = (JT -JBORD(JBL))/(JLAST - JBORD(JBL))
              damp = 2.d0*damp**3 -3.d0*damp**2 +1.d0
c              WRITE(KO,*) 'damp',jt,damp
            endif
c
            CALL SPLINT(T,NEXT,IJ,NIJ,AI)
C       IF(XSTRAC.GE.4) WRITE(KO,*) ' JT,T,IJ,NIJ,AI =',DBLE(JT),DBLE(T)
C    #           ,IJ, NIJ, (DBLE(AI(I)),I=1,NIJ)
            DO 80 I=1,NIJ
               IBJ = IJ + I - 1
            DO 80 LJ=1,NLJ
               LP = LJ-1 + NINT(JT - LJMAX)
                 IF(LP.LT.0 .OR. LP.GT.MAL) GO TO 80
                 O = 2*(CSIG(LP+1,ITX) - CSIG(1,ITX))
                 C6 = AI(I) * EXP(CI*O) * damp
CDIR$ IVDEP
            DO 70 IAM=1,MAM
70          AMPL(LP*MAM+IAM) = AMPL(LP*MAM+IAM) +
     &          C6    * AMPL(LB + IAM + LJ*MAM + IBJ*NLJ*MAM)
80          CONTINUE
100      CONTINUE
         JPREV = JLAST
C                 SHIFT LAST JT-SET OF VALUES TO START NEXT BLOCK
         DO 120 LJ=1,NLJ
            IBT = LB +       LJ*MAM +  1  *NLJ*MAM
            IBF = LB +       LJ*MAM + NEXT*NLJ*MAM
C           IF(XSTRAC.GE.4) WRITE(KO,*) ' SHIFT',LJ,' FROM',NEXT,IBF,IBT
CDIR$ IVDEP
         DO 120 IAM=1,MAM
 120     AMPL(IBT + IAM) = AMPL(IBF + IAM)
         NEXT = 1
C
23    CONTINUE
C
      IF(LAMPL.NE.IC) GO TO 230
C:    IF(LAMPL.EQ.0 ) GO TO 230
      WRITE(LFAM-1,224)
224   FORMAT('0Scattering Amplitude: Legendre coefficients. ')
      DO 225 LP1=1,MAL1
         LP = LP1-1
225   WRITE(LFAM-1,226) LP,(AMPL(LP*MAM+IAM),IAM=1,MAM)
        if(LAMPL.eq.IC) written(LFAM-1) = .true.
226   FORMAT(' LP=',I5,1X,1P,10E12.4 / (10X,10E12.4))
230   NANGL0 = (abs(THMAX) - THMIN)/abs(THINC) + 1 + 0.5
	NANGL0 = max(0,NANGL0)
	if(.not.gettheoryplot) NANGL0=0
      NANGL = NANGL0
      EXTRA = NANGLE>0
      IF(EXTRA) NANGL = NANGL0 +  NANGLE

      IFOUT=200+ITX
      IFOG = IFOUT+100
      IFO = ITX<=10
      IF(IC.EQ.PEL .AND. IA.EQ.EXL) then
	  if(IFOUT/=201)  IFOUT=200
	  IFO=.true.
	endif
      if(NANGL>0) call openif(16)
      if(IFO.and.NANGL>0) then
	call openif(IFOUT)
        if(NMULTIES>0.and.IA>1) call openif(100+IFOUT)
        endif

      DO 290 NEARF = 1,ABS(mod(NEARFA,10))
         IF(NEARFA.GT.0 .AND. NEARF.GT.1
     X                  .AND. (IC.NE.PEL.OR.IA.NE.EXL)) GO TO 290
      IF(ABS(NEARFA).GT.1) WRITE(KO,2301) SIDES(NEARF)
2301  FORMAT(/ ' Scattering from  ',A12/)
         RX =  2/PI * SIGN(1.0, NEARF-2.5)
      IF(KQ1PR.GT.0.and..not.DESCR) THEN
             WRITE(LCROSS,2302) NANGL,KQ1PR,PP,PPK,HDGQ
             IF(NEARF.EQ.1) WRITE(LCROSS,2304)
		written(LCROSS) = .true.
            DESCR = .true.
          ENDIF
      IF(KQ1PR.GT.0) THEN
             if(IFO) WRITE(IFOUT,2302) NANGL,KQ1PR,PP,PPK,HDGQ
	     if(IFO) written(IFOUT) = .true.
            LS = MOD(LEG,5)+1
            WRITE(LCROSS,2303) '@',LEG,IC,IA,NEARF
	    WRITE(LCROSS,23034) '#',LEN-1,ENLAB
	    if(IFO) WRITE(IFOUT,23034) '@',LEN-1,ENLAB
            if(IFO) WRITE(IFOUT,2303) '#',0,IC,IA,NEARF
          IF(LEG.lt.10) WRITE(LCROSS,23031) LEG,LS
          IF(LEG.ge.10.and.LEG<100) WRITE(LCROSS,23032) LEG,LS
          IF(LEG>=100) WRITE(LCROSS,23033) LEG,LS
          if(DGAM==0) then
            WRITE(LCROSS,2305) PPK
            if(IFO) WRITE(IFOUT,2305) PPK
          else
            WRITE(LCROSS,23051)
            if(IFO) WRITE(IFOUT,23051)
          endif
            LEG = LEG + 1
        ENDIF
2302  FORMAT('#',I4,' angles,',I3,' tensor ranks for ',I1,'=',A10,/,
     X       '@subtitle ',A82,/,'@legend ON',/,
     X       '@legend x1 0.2',/,'@legend y1 0.8',/,'@g0 type LOGY'/
     x       '@yaxis ticklabel format power')
2303  FORMAT(a1,'legend string',I4,' "Partition=',I3,' Excit=',I3,
     X       ' near/far=',I2,'"')
23031 FORMAT('@s',I1,' linestyle ',I3)
23032 FORMAT('@s',I2,' linestyle ',I3)
23033 FORMAT('@s',I3,' linestyle ',I3)
23034 FORMAT(a1,'legend string',I4,' "Lab energy =' ,f10.4,'"')
2304  FORMAT('@xaxis label "Scattering angle (degrees)"',/,
     X       '@yaxis label "Cross section (mb/sr)"')
2305  FORMAT('#  Theta       sigma       iT11        T20',
     X        '         T21         T22         Kyy for ',A10)
23051 FORMAT('#  Theta       sigma       iT10        T20',
     X        '        iT30         T40        iT50')
23021 FORMAT('#',I4,' angles, gamma decay distributions from state ',i2,
     x   ' to state ',I2,','/
     X   '#   using ',I1,' tensor ranks, for nucleus ',I1,'=',A10,/,
     X   '@subtitle ',A82,/,'@legend ON')

      IF(XSTRAC*LXSEC.LT.0) WRITE(LXSEC,*) IC,IA,KQ1PR,NANGL,NEARF
      IF(ABS(LAMPL).eq.IC) WRITE(LFAM,*) JIX,JEX,NANGL,NEARF,ENLAB
      if(ABS(LAMPL).eq.IC) written(LFAM) = .true.
      IF(CDCC/=0.and.IA>1) WRITE(57,'(f4.1,i4,f8.4,i4)') 
     X             JEX(1),IP,PENEX,IBIN
      IF(INFAM/=0) READ(abs(INFAM),*) R,R,R,R,I,I
      IF(OUTFAM/=0) then
	WRITE(abs(OUTFAM),*) JIX,JEX,NANGL,NEARF
	written(abs(OUTFAM)) = .true.
	endif

      SIGG(:) = 0.
      SUMXS = 0.0
      DO 280 ITH=1,NANGL
         IF(ITH.LE.NANGL0) THETA = (ITH-1)*abs(THINC) + THMIN
         IF(ITH.GT.NANGL0) THETA = ANGLES(ITH-NANGL0)
         IF(THETA.LT.0.0 .OR. THETA.GT.180.) GO TO 279
         IF(THETA.LE.0.01.and.IC.EQ.PEL .AND. IA.EQ.EXL) THETA=0.01
         IF(THETA.GE.180-.01.and.IC==PEL .AND. IA==EXL) THETA=180-.01
      IF(EXCH) THETA=180d0-THETA
      IF(ABS(LAMPL).eq.IC) WRITE(LFAM,*) THETA
      IF(OUTFAM/=0) WRITE(abs(OUTFAM),*) THETA
        TH = THETA * PI/180.0
        CTH = COS(TH)
        STH = SIN(TH)
        MP = max(1,MIN(MAL,NINT(JIX(1)+JIX(2)+JEX(1)+JEX(2))))
      CALL PLM(CTH,MAL,MP,LMAX1,PL)
      IF(NEARF.GT.1) CALL QLM(CTH,MAL,MP,LMAX1,PL(1,1,2))
C
      FCOUL = 0.; FCOUL2 = 0.
      IF( IC.EQ.PEL .AND. IA.EQ.EXL) THEN
         O = -2*ETA(PEL,EXL)*LOG(SIN(TH*0.5))
         FCOUL = -(ETA(PEL,EXL)/(2*K(PEL,EXL))) * SIN(TH/2.0)**(-2)
     X           * EXP(CI*O)
         IF(IEXCH.NE.0) THEN
            T = PI - TH
            O = -2*ETA(PEL,EXL)*LOG(SIN(T*0.5))
            FCOUL2 =  IEXCH *
     X      (-ETA(PEL,EXL)/(2*K(PEL,EXL))) * SIN(T/2.0)**(-2)
     X           * EXP(CI*O)
            ENDIF
         ENDIF
2309  MAA = 0
      DO 231 IAM=1,MAM
231   FAM(IAM) = 0.0
      IF(INFAM/=0) then 
      READ(abs(INFAM),*) THIN
       if(abs(THIN-THETA).gt.0.1) 
     X    write(0,*) 'ANGLE MISMATCH in reading file',INFAM,
     X   ': theta =',real(THETA),' expected, but ',real(THIN),' found'
      IF(INFAM>0) READ(abs(INFAM),232) (FAM(IAM),IAM=1,MAM)
      IF(INFAM<0) READ(abs(INFAM),232) C9
232   FORMAT(6E12.4)
 	endif
      IAM = 0
      DO 235 IMA = 1,NSA
      MA = -JIX(1) + IMA-1
      DO 235 IJA = 1,NJA
      MJA = -JIX(2) + IJA-1
      DO 235 IMB = 1,NSB
         MB = -JEX(1) + IMB - 1
      DO 235 IJB = 1,NJB
      MJB = -JEX(2) + IJB-1
         IAM = IAM + 1
      MP = NINT(MAA + MA + MJA - MB - MJB)
      MPA = ABS(MP)
	 if(MPA>MAXPLM) stop 'MAXPLM'
      C6D = 0.0
! NUCLEAR:
      IF(NEARF.EQ.1) THEN
      DO 234 LP=MPA,MAL
234   C6D = C6D + AMPL(LP*MAM + IAM) * PL(LP+1,MPA+1,1)
        ELSE
      DO 2345 LP=MPA,MAL
2345  C6D = C6D + AMPL(LP*MAM + IAM) *
     &             0.5 * CMPLX(PL(LP+1,MPA+1,1),RX*PL(LP+1,MPA+1,2))
      ENDIF

! COULOMB:
      IF(IMA.EQ.IMB.AND.IJA.EQ.IJB) then   !diagonal on m-values

      if(abs(NEARFA)<=11) then   ! all Coulomb is near-side, so for all but far-side
         if(NEARF.NE.2) C6D = C6D + FCOUL

!   split also the Coulomb amplitude if near or far wanted separately:
      else if(NEARF.eq.1) then   ! add all Coulomb
         C6D = C6D + FCOUL
      else if(NEARF.eq.2) then   ! add far-side Coulomb
         iSgn= 1
         call Fuller(ETA(PEL,EXL),TH,iSgn,far)
         C6D = C6D + far*FCOUL
      else if(NEARF.eq.3) then   ! add near-side Coulomb
         iSgn=-1
         call Fuller(ETA(PEL,EXL),TH,iSgn,near)
         C6D = C6D + near*FCOUL
        endif

      endif

      IF(IMA.EQ.IJB.AND.IJA.EQ.IMB .AND. NEARF.NE.2)  ! exchange (if any)
     &                              C6D = C6D + FCOUL2
235   FAM(IAM) = FAM(IAM) + C6D
C
      DO 248 KQ2 = 1,KQMAX2
      DO 248 LQ2 = 2-KQ2,KQ2
         KQE = KQ2 - 1
         LQE = LQ2 - 1
!	if(ITH==1) WRITE(KO,*) 'KQE,LQE =',KQE,LQE

      DO 248 KQ1 = 1,KQMAX1
      DO 248 LQ1 = 1,KQ1
         KQ = KQ1 - 1
         LQ = LQ1 - 1
        FA(KQ1,LQ1,KQ2,LQ2) = 0
!	if(ITH==1) WRITE(KO,*) '        KQ,LQ =',KQ,LQ
      IAM = 0
      DO 245 IMA = 1,NSA
         MA = -JIX(1) + IMA-1
      DO 245 IJA = 1,NJA
         MJA = -JIX(2) + IJA-1
      DO 242 IMB=1,NSB
         MB = -JEX(1) + IMB - 1
      DO 242 IJB = 1,NJB
      MJB = -JEX(2) + IJB-1
        IAM = IAM + 1
       if(PP.eq.0) then
            MAP = LQ + MA
            IAMP= IAM + LQ * NJA*NSB*NJB
            IF(ABS(MAP).GT.JIX(1)) GO TO 242
                 R = CLEB6(JIX(1),MA,KQ+Z,LQ+Z,JIX(1),MAP)
       else if(PP.eq.1) then
            MJAP = LQ + MJA
            IAMP= IAM + LQ *     NSB*NJB
            IF(ABS(MJAP).GT.JIX(2)) GO TO 242
                 R = CLEB6(JIX(2),MJA,KQ+Z,LQ+Z,JIX(2),MJAP)
       else if(PP.eq.2) then
            MBP = LQ + MB
            IAMP= IAM + LQ *         NJB
            IF(ABS(MBP).GT.JEX(1)) GO TO 242
                 R = CLEB6(JEX(1),MB,KQ+Z,LQ+Z,JEX(1),MBP)
       else if(PP.eq.3) then
            MJBP = LQ + MJB
            IAMP= IAM + LQ 
            IF(ABS(MJBP).GT.JEX(2)) GO TO 242
                 R = CLEB6(JEX(2),MJB,KQ+Z,LQ+Z,JEX(2),MJBP)
       else if(PP.eq.4) then
            MAP = LQ + MA
            IAMP= IAM + LQ * NJA*NSB*NJB
            IF(ABS(MAP).GT.JIX(1)) GO TO 242
                 R = CLEB6(JIX(1),MA,KQ+Z,LQ+Z,JIX(1),MAP)
            MBP = LQE + MB
            IAMP= IAMP + LQE *         NJB
            IF(ABS(MBP).GT.JEX(1)) GO TO 242
                 R = R * CLEB6(JEX(1),MB,KQE+Z,LQE+Z,JEX(1),MBP)
       endif
      C6 = FAM(IAM)
      C7 = FAM(IAMP)
      FA(KQ1,LQ1,KQ2,LQ2) = FA(KQ1,LQ1,KQ2,LQ2) 
     x             + C7 * CONJG(C6)  *  R * SQRT((2*KQ+1.)*(2*KQE+1))
C       IF(XSTRAC.GE.1.AND.ITH.EQ.10) WRITE(KO,240) KQ,LQ,MA,MAP,MJA,
C    #   MB,MJB,     FCOUL,C6,C7,R,   FA(KQ1,LQ1,KQ2,LQ2)
C240    FORMAT(' ',2I3,5F5.1,    7F10.4,2F13.5)
242   CONTINUE
245     CONTINUE
248     CONTINUE
      SNZ = MAX(ABS(FA(1,1,1,1)),1E-9+Z)
      DO 250 KQ1=2,KQMAX1
      DO 250 LQ1=1,KQ1
         KQ = KQ1-1
         C9 = 1
         IF(KQ.NE.(KQ/2)*2) C9 = (0.,1.)
250    XSEC(KQ1,LQ1)=FA(KQ1,LQ1,1,1) *C9/SNZ
      XSEC(1,1) = FA(1,1,1,1) * XCOEF
      KYY = 0.
      if(IKYY==1) KYY = sqrt(2./3.)*(FA(2,2,2,2)+FA(2,2,2,0))/SNZ

      T = 100.
      IF(XSEC(1,1).GT.1E-20) T = LOG10(XSEC(1,1))
      if(DGAM>0) then
       IF(ABS(T).GT.4) THEN
        WRITE(KO,256) THETA,(XSEC(KQ1,1),KQ1=1,KQMAX1)
       ELSE
        WRITE(KO,255) THETA,(XSEC(KQ1,1),KQ1=1,KQMAX1)
       ENDIF
      else if(IKYY==0) then
       IF(ABS(T).GT.4) THEN
        WRITE(KO,256) THETA,((XSEC(KQ1,LQ1),LQ1=1,KQ1),KQ1=1,KQMAX1)
       ELSE
        WRITE(KO,255) THETA,((XSEC(KQ1,LQ1),LQ1=1,KQ1),KQ1=1,KQMAX1)
       ENDIF
      else
       IF(ABS(T).GT.4) THEN
        WRITE(KO,2561)THETA,
     x            ((XSEC(KQ1,LQ1),LQ1=1,KQ1),KQ1=1,min(2,KQMAX1)),KYY,
     x            ((XSEC(KQ1,LQ1),LQ1=1,KQ1),KQ1=3,KQMAX1)
       ELSE
        WRITE(KO,2551)THETA,
     x            ((XSEC(KQ1,LQ1),LQ1=1,KQ1),KQ1=1,min(2,KQMAX1)),KYY,
     x            ((XSEC(KQ1,LQ1),LQ1=1,KQ1),KQ1=3,KQMAX1)
       ENDIF
      endif
255   FORMAT(1X,F8.2,' deg.: X-S =', F13.6,' mb/sr,',:,22X,' & pols =',
     & F3.0,F9.5,3X,3F9.5,:,/,78X,4F9.5,:,/,78X,5F9.5,
     &                    :,/,72X,6F9.5,:,/,60X,7F9.5)
256   FORMAT(1X,F8.2,' deg.: X-S =',1P,E13.4,' mb/sr,',:,22X,' & pols ='
     &,0P,F3.0,F9.5,3X,3F9.5,:,/,78X,4F9.5,:,/,78X,5F9.5,
     &                       :,/,60X,6F9.5,:,/,42X,7F9.5)
2551  FORMAT(1X,F8.2,' deg.: X-S =', F13.6,' mb/sr,',:,22X,' & pols =',
     & F3.0,2F9.5,3X,3F9.5,:,/,78X,4F9.5,:,/,78X,5F9.5,
     &                     :,/,72X,6F9.5,:,/,60X,7F9.5)
2561  FORMAT(1X,F8.2,' deg.: X-S =',1P,E13.4,' mb/sr,',:,22X,' & pols ='
     &,0P,F3.0,2F9.5,3X,3F9.5,:,/,78X,4F9.5,:,/,78X,5F9.5,
     &                        :,/,60X,6F9.5,:,/,42X,7F9.5)
      IF(KQMAX1.GE.3.AND.KOORDS.GE.1) THEN
C      Calculate Transverse Tensors TT(k,0) for k=1,2,3
         IF(KQMAX1.GE.3) TENS(2) = XSEC(2,2) * SQRT(2.)
         IF(KQMAX1.GE.3) TENS(3) =
     &         - 0.5 * (XSEC(3,1) + SQRT(6.) * XSEC(3,3))
         IF(KQMAX1.GE.4) TENS(4) =
     &         - 0.5 * (SQRT(3.)*XSEC(4,2) + SQRT(5.)*XSEC(4,4))
         WRITE(KO,259) THETA,(TENS(KQ1),KQ1=2,MIN(KQMAX1,4))
259      FORMAT(1X,F8.2,' deg.:',57X,'TTk0 =',4G12.4)
       ENDIF
      IF(KQMAX1.GE.3.AND.KOORDS.GE.2) THEN
C      Calculate Recoil Tensors T(2,q) for q=0,1,2  (Second rank only)
C      Start by finding Hooton-Johnson second-rank tensors:
       TENHJ(0+1) = SQRT(3./8.) * (1.-CTH) * XSEC(3,3)
     X              - SQRT(1.5) * STH      * XSEC(3,2)
     X              + 0.25*(1.+3.*CTH)     * XSEC(3,1)
       TENHJ(1+1) = - 0.5 * STH            * XSEC(3,3)
     X              + CTH                  * XSEC(3,2)
     X              + SQRT(3./8.)*STH      * XSEC(3,1)
       TENHJ(2+1) = 0.25*(3. + CTH)        * XSEC(3,3)
     X              + 0.5 * STH            * XSEC(3,2)
     X              + SQRT(3./32.)*(1.-CTH)* XSEC(3,1)
C     Now calculate recoil tensors:
       TENS(0+1) = -0.5 * (TENHJ(1) - SQRT(6.) * TENHJ(3))
       TENS(1+1) = - TENHJ(1+1)
       TENS(2+1) = SQRT(3./8.) * TENHJ(1) + 0.5 * TENHJ(3)
       IF(KQMAX1.GE.5) THEN
            CHTH = COS(TH*0.5)
            SHTH = SIN(TH*0.5)
            TENS(4+1) = ((35.*SHTH**4 - 30.*SHTH**2 + 3.)*XSEC(5,1)
     X                  +(7.*SHTH**3 - SHTH)*CHTH * XSEC(5,2)
     X                  + (7.*SHTH**2-1.)*CHTH**2 /6. * XSEC(5,3)
     X                  + CHTH**3*SHTH/6. * XSEC(5,4)
     X                  + CHTH**4/48. * XSEC(5,5)   ) / 8.0
        ENDIF
       WRITE(KO,2592) THETA,(TENS(LQ1),LQ1=1,3),(TENS(5),LQ1=5,KQMAX1)
C2592   FORMAT(65X,'Recoil: T2q =',3F9.5,:,' T40 =',F9.5)
2592  FORMAT(1X,F8.2,' deg.:',50X,'Recoil: T2q =',3F9.5,:,' T40 =',F9.5)
      ENDIF
      IF(KQMAX1.GE.3.AND.KOORDS.GE.3) THEN
C      Print Hooton-Johnson second-rank tensors from above:
       WRITE(KO,2593) (TENHJ(LQ1),LQ1=1,3)
2593   FORMAT(69X,'HJ: T2q =',3F9.5)
      ENDIF
      SUMXS = SUMXS  + XSEC(1,1) * STH
      RUTHR = -1.
      RUTH = 1.
      XEL = XSEC(1,1)
      IF(MODSQ(FCOUL).GT.1D-10) THEN
        RUTH = MODSQ(FCOUL) 
        RUTHR = FA(1,1,1,1) / (RUTH * NSA * NJA)
        WRITE(KO,257) RUTHR
        XEL = RUTHR
257     FORMAT('+',40X, ' /R =',1P,E13.4)
        RUTH = RUTH*10  ! to mb
       ENDIF
       if(THMAX.lt.0.) XEL = XSEC(1,1)
C
	if(NEARF==1) then
		do 2575 IK=1,KA   ! Calculate chi-sq for data
		 id = KAD(IK)
		 ib = data_ib(id)
		 if(ib>0) go to 2575  ! only do particles here
		 kq1 = data_rank_k(id)+1
		 lq1 = data_rank_q(id)+1
		 idir = data_idir(id)
	!		Adjust any datanorm search parameter!
		   datanorm=1.0
		   do ip=1,nvars
		   if(srch_kind(ip)==5.and.srch_datanorm(ip)==id) 
     x	     datanorm = datanorm * srch_value(ip)
		   enddo
!           WRITE(KO,*) ' dataset ',id,kq1,lq1,real(datanorm)
		   if(lq1<=kq1)  theory = XSEC(kq1,lq1)
		   if(kq1==2.and.lq1==3) theory = KYY
!           call flush(6)
		   if(idir>=1) theory = theory/RUTH
			theorycm = theory
		   if(idir==0.and.kq1==1.and.data_lab(id)) 
     x         theory = theory / CM_LAB(ITH-NANGL0) ! convert to lab

		sfactor = exp(2d0*PI*ETA(PEL,EXL)) * ECMI
		do ip=1,datalen(id)
		  if(IANGLS(ip,id)>=0) then       
		  if(ITH == NANGL0+IANGLS(ip,id)) then       
	!           WRITE(KO,*) ' XS',datangles(ip,id),THETA,ITH,IANGLS(ip,id)
           if(idir==2) then  !  convert to ratio to Rutherford, the first time
              datavals(ip,id) = datavals(ip,id)/RUTH
              dataerr(ip,id) = dataerr(ip,id)/RUTH
             endif
           if(idir==-1) then  !  convert to absolute, the first time
              datavals(ip,id) = datavals(ip,id)/sfactor
              dataerr(ip,id) = dataerr(ip,id)/sfactor
             endif

            chi = (theory/datanorm-datavals(ip,id))/dataerr(ip,id)
            data_chisq(id) = data_chisq(id) + chi**2
	    theoryvals(ip,id) = theory
             call flush(6)
           else if(ITH<=NANGL0.and.data_type(id)<=2) then
             theoryplot(ITH,id) = theorycm
!	     write(190,281) ith,id,ip,theory,enlab
!281	     format(' theory ',3i3,'=',f8.5,' at E =',f8.3)
           endif
           endif
           
          enddo
2575	continue
	endif  ! NEARF==1
C---       Print Rutherford-ratio or sigma, in simple forms for plotting
c     IF(THETA.LT.0.01 .AND. IC.EQ.PEL .AND. IA.EQ.EXL) GO TO 260
      IF(KQ1PR.GT.0) WRITE(LCROSS,258)
     1     THETA,XEL,((XSEC(KQ1,LQ1),LQ1=2-MOD(KQ1,2),KQ1),KQ1=2,KQ1PR),
     2     (KYY,I=1,IKYY)
      IF(KQ1PR.GT.0.and.IFO)  WRITE(IFOUT,258)
     1     THETA,XEL,((XSEC(KQ1,LQ1),LQ1=2-MOD(KQ1,2),KQ1),KQ1=2,KQ1PR),
     2     (KYY,I=1,IKYY)
c258   FORMAT(1P,6E12.4)
c258   FORMAT(6G12.5)
258   FORMAT(10G12.4)
      if(CDCC/=0.and.IA>1) WRITE(57,265) (FAM(IAM)*BSIGN,IAM=1,MAM)
      IF(XSTRAC*LXSEC.LT.0) WRITE(LXSEC,*) THETA,XEL
C---
260   IF(ITH.EQ.1 .AND. RUTHR.LT.0..AND.NEARF.EQ.1) WRITE(KO,261) SIGR
261   FORMAT('+',40X, ' REAC ', 1P,E12.4,' mb' )
      IF(KQ1PR.EQ.0 .AND. SIGR .EQ.0. .AND.INFAM==0.and.
     x   .not.(IC.EQ.PEL .AND. IA.EQ.EXL) ) GO TO 290
      IF(KQMAX1.GE.4 .OR. KOORDS.GE.1) WRITE(KO,256)
      if(OUTFAM/=0) WRITE(abs(OUTFAM),265) (FAM(IAM),IAM=1,MAM)
      IF(ABS(LAMPL).NE.IC) GO TO 279
      WRITE(LFAM,265) (FAM(IAM),IAM=1,MAM)
265   FORMAT(1P,6E12.4)
279   if(ITH==NANGL0) then
         IF(KQ1PR.GT.0) WRITE(LCROSS,*) 'END'
         if(KQ1PR.GT.0.and.IFO) WRITE(IFOUT,*) 'END'
      endif

       if(NMULTIES>0.and.ITH.LE.NANGL0) then
         do KQ1=1,KQMAX1
          TKQ=1.0
          if(KQ1>1) TKQ = XSEC(KQ1,1)
          SIGG(KQ1) = SIGG(KQ1) + XSEC(1,1)*TKQ*STH
         enddo
         endif
280   CONTINUE

      IF(ITH.GT.1 .AND. SUMXS.NE.0.)
     #    WRITE(KO,293) 2*PI * SUMXS * (abs(THINC)*PI/180.) ,
     #    THMIN,abs(THMAX),ENLAB
 293   FORMAT(41X, '  Integrated',f10.4,' mb, over [',
     x        f8.3,',',f8.3,'] at ',f10.4,' MeV ')
         IF(KQ1PR.GT.0.and.EXTRA) then
	    WRITE(LCROSS,*) 'END'
            if(IFO) WRITE(IFOUT,*) 'END'
	    endif

         if(NMULTIES>0.and.IA>1) then
        write(KO,1471) DNAME,IA,(DSPINS(I),DMULTIES(I),I=1,NMULTIES)
1471    format(/' GAMMA DECAY CROSS SECTIONS mb/sr of ',A8,' level',i3,
     X     ' as if 100% to states/multipoles:', 10(f5.1,'/',i2,';'))
         IN = PP-1
         do KQ1=1,KQMAX1
          KQ = KQ1-1
          SIGG(KQ1) = SIGG(KQ1) * abs(THINC)*PI/180.  *2d0*PI
          do IFF=1,NMULTIES
           LGAM=DMULTIES(IFF)
           RCOEF(KQ1,IFF)= sqrt(2d0*JEX(IN)+1d0) * (2*LGAM+1) 
     x      * (-1)**nint(JEX(IN)-DSPINS(IFF)+KQ+1) 
     x      * cleb6(LGAM+Z,1d0,LGAM+Z,-1d0,KQ+Z,Z)
     x      * racah(LGAM+Z,LGAM+Z,JEX(IN),JEX(IN),KQ+Z,DSPINS(IFF))
          enddo
         enddo
	  write(KO,1472) SIGG(1:KQMAX1)
!1472      format(/'     Multipole integrals:',10(f10.5,f4.1))
 1472      format(/'     Multipole integrals:',20f10.5)
          do IFF=1,NMULTIES
	  write(KO,1473) IFF,RCOEF(1:KQMAX1,IFF)
1473      format( '     R coefficients to',i2,':',20f10.5)
	  enddo
	  write(KO,*)
           if(IFO) WRITE(IFOG,23021) NANGL0,IA,IFF,KQMAX1,PP,PPK,HDGQ
           if(IFO) written(IFOG) = .true.
            LS = MOD(LEG,5)+1
            WRITE(LCROSS,2303) '@',LEG,IC,IA,NEARF
            WRITE(LCROSS,23034) '#',LEN-1,ENLAB
            if(IFO) WRITE(IFOG,23034) '#',LEN-1,ENLAB
            if(IFO) WRITE(IFOG,2303) '@',LEG,IC,IA,NEARF
          IF(LEG.lt.10) WRITE(LCROSS,23031) LEG,LS
          IF(LEG.ge.10.and.LEG<100) WRITE(LCROSS,23032) LEG,LS
          IF(LEG>=100) WRITE(LCROSS,23033) LEG,LS
            WRITE(LCROSS,23052) 
            if(IFO) WRITE(IFOG,23052) 
23052       FORMAT('#  Theta       sigma (decay gamma)')
            LEG = LEG + 1

            TGXSECS(:) = 0.
         DO 288 ITH=1,NANGL
            IF(ITH.LE.NANGL0) THETA = (ITH-1)*abs(THINC) + THMIN
            IF(ITH.GT.NANGL0) THETA = ANGLES(ITH-NANGL0)
            IF(THETA.LT.0.0 .OR. THETA.GT.180.) GO TO 288
             TH = THETA * PI/180.0
             CTH = COS(TH)
             STH = SIN(TH)
             CALL PLM(CTH,KQMAX1-1,0,KQMAX1,PLEG)
            do IFF=1,NMULTIES
            GXSECS(IFF) = 0.
             do KQ1=1,KQMAX1
              GXSECS(IFF) = GXSECS(IFF) + 
     x            RCOEF(KQ1,IFF) * PLEG(KQ1)*SIGG(KQ1)/(4d0*PI)
             enddo
            enddo
           write(KO,285) THETA,GXSECS
285       FORMAT(1X,F8.2,' deg. GX-S =',10F13.6)
   
                    WRITE(LCROSS,258) THETA,GXSECS
           IF(IFO)  WRITE(IFOG,258) THETA,GXSECS
            TGXSECS(:) = TGXSECS(:) + GXSECS(:)*STH

C
		if(NEARF==1) then
		do 287 IK=1,KA   ! Calculate chi-sq for data
		 id = KAD(IK)
		 kq1 = data_rank_k(id)+1
		 lq1 = data_rank_q(id)+1
		 ib = data_ib(id)
		  iff=0
		  do I=1,NMULTIES
		    if(DLEVEL(I)==ib) IFF=I
			enddo
		 if(iff==0) go to 287   ! only do gammas here
	!		Adjust any datanorm search parameter!
		   idir = data_idir(id)
		   datanorm=1.0
		   do ip=1,nvars
		   if(srch_kind(ip)==5.and.srch_datanorm(ip)==id) 
     x	     datanorm = datanorm * srch_value(ip)
		   enddo
	!           WRITE(KO,*) ' dataset ',id,kq1,lq1,real(datanorm)
			   theory = GXSECS(IFF)
			   if(idir>=1) theory = theory/RUTH
		   theorycm = theory
!			   if(idir==0.and.kq1==1.and.data_lab(id)) 
!    x         theory = theory / CM_LAB(ITH-NANGL0) ! convert to lab

		sfactor = exp(2d0*PI*ETA(PEL,EXL)) * ECMI
			do ip=1,datalen(id)
		  if(IANGLS(ip,id)>=0) then       
		  if(ITH == NANGL0+IANGLS(ip,id)) then       
	!           WRITE(KO,*) ' XS',datangles(ip,id),THETA,ITH,IANGLS(ip,id)
			   if(idir==2) then  !  convert to ratio to Rutherford, the first time
				  datavals(ip,id) = datavals(ip,id)/RUTH
				  dataerr(ip,id) = dataerr(ip,id)/RUTH
				 endif
			   if(idir==-1) then  !  convert to absolute, the first time
				  datavals(ip,id) = datavals(ip,id)/sfactor
				  dataerr(ip,id) = dataerr(ip,id)/sfactor
				 endif

				chi = (theory/datanorm-datavals(ip,id))/dataerr(ip,id)
				data_chisq(id) = data_chisq(id) + chi**2
			theoryvals(ip,id) = theory
				 call flush(6)
			   else if(ITH<=NANGL0.and.data_type(id)<=2) then
				 theoryplot(ITH,id) = theorycm
!	   			 write(190,281) ith,id,ip,theory,enlab
!281		     format(' theory ',3i3,'=',f8.5,' at E =',f8.3)
			   endif
			   endif
		   
			  enddo
287		continue
		endif  ! NEARF==1

288	  continue  ! THETA loop for gamma angle

      IF(ITH.GT.1)WRITE(KO,289) 2*PI*TGXSECS(:)*(abs(THINC)*PI/180) 
289   FORMAT(/'    Angle integrated:',10F13.6)
      
         WRITE(KO,*)
         WRITE(LCROSS,*) 'END'
         WRITE(IFOG,*) 'END'
 	endif

290    continue
300   deallocate (AMPL)
301   deallocate (YSIG)

      return
      END
	subroutine lab2cm(TH,A1,A2,A3,A4,ECMI,ECMF,XCM_LAB,pi)
	implicit real*8 (a-h,o-z)
		
	x=sqrt(A1*A3*ECMI/(A2*A4*ECMF))
	
	thetalab = TH*pi/180.d0
	coslab=cos(thetalab)
        c2=coslab**2
        s2=1d0-c2
*
c As 0 < TH < 180 deg (at most), we have always sinlab > 0, so it is 
c ok to take only the positive square root of s2.
*
        sinlab=sqrt(s2)
*
        thetacm=thetalab+asin(x*sinlab)
*
        coscm=cos(thetacm)
        thetacm=thetacm*180.d0/pi !switching to degrees
*
        a=1.+x*coscm
        b=abs(a)
	XCM_LAB=b/(1.+x*x+2.*x*coscm)**1.5
!        sigCM=sigLAB*XCM_LAB

	TH = thetacm   ! return cm degrees
	return
	end

C **** SUBROUTINE FULLER START *****************************
C **** PURPOSE                                             *
C Calculates the ratio (CuRat) between the Near (iSgn=-1)  *
C or Far (iSgn=1) Coulomb Amplitude and the Rutherford     *
C Amplitude at Ang (in Rad)                                *
C **** ARGS                                                *
C (eta,Ang,iSgn,CuRat)                                     *
C **** INPUTS                                              *
C    eta  Sommerfeld Parameter                             *
C    Ang  Angle in Rad                                     *
C   iSgn  Side Index                                       *
C    -1   Near-Side                                        *
C    +1   Far-Side                                         *
C **** OUTPUT                                              *
C cuRat   f_Ruth({N, F})/f_Ruth                            *
C **** SUBPROGRAM USED                                     *
C FUNCTION cPSI                                            *
C
      SUBROUTINE FULLER(eta,Ang,iSgn,CuRat)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8 (A-B,D-H,O-z)
      deta = eta
      dAng=Ang
      ceta=DCMPLX(0.D0,deta)
      ci=(0.D0,1.D0)
      dpai=2.D0*ASIN(1.D0)
      dC2=COS(dAng/2.D0)**2
      dS2=SIN(dAng/2.D0)**2
C calculate the Fuller S(theta) function
        IF (dC2.LT.dS2) THEN
        k=0
        c1=cPsi(0.D0)-cPsi(deta)-LOG(dC2)
        c2=DCMPLX(1.D0,0.D0)
        cAdd=c1*c2
        cSum=cAdd
        kont=0
        kGo=1
        dMax=0.D0
          DO WHILE (kGo.EQ.1)
          dk=DFLOAT(k+1)
          c1=c1+1.D0/dk-1.D0/(dk+ceta)
          c2=(dk+ceta)/dk*dC2*c2
          cAdd=c1*c2
          IF (ABS(cAdd).GT.dMax) dMax=ABS(cAdd)
          cSum=cSum+cAdd
            IF (ABS(cAdd/cSum).LT.1.D-14) THEN
            kont=kont+1
            IF (kont.EQ.5) kGo=0
            ELSE
            kont=0
            ENDIF
          k=k+1
          ENDDO
        ELSE
C Direct Formula
        k=0
        c2=DCMPLX(1.D0,0.D0)
        cAdd=DCMPLX(1.D0,0.D0)
        cSum=cAdd
        kont=0
        kGo=1
        dMax=0.D0
          DO WHILE (kGo.EQ.1)
          dk=DFLOAT(k+1)
          c2=(dk+ceta)/dk*dS2*c2
          cAdd=(dk+ceta)/(dk+1.D0+ceta)*dS2*cAdd
          IF (ABS(cAdd).GT.dMax) dMax=ABS(cAdd)
          cSum=cSum+cAdd
            IF (ABS(cAdd/cSum).LT.1.D-14) THEN
            kont=kont+1
            IF (kont.EQ.5) kGo=0
            ELSE
            kont=0
            ENDIF
          k=k+1
          ENDDO
        cSum=cSum/(1.D0+ceta)
        ENDIF
      cSum=ci/(2.D0*dpai)*EXP((1.D0+ceta)*LOG(dS2))*cSum
      dTmp=EXP(-2*dpai*deta)
        IF (iSgn.EQ.-1) THEN
C Near 
        CuRat=1.D0/(1.D0-dTmp)-cSum
        ELSEIF (iSgn.EQ.1) THEN
C Far 
        CuRat=-dTmp/(1.D0-dTmp)+cSum
        ELSE
C Set Zero
        CuRat=(0.D0,0.D0)
        ENDIF
       RETURN
       END
C **** SUBROUTINE FULLER END *****************************
C>>>
C<<<
C **** FUNCTION cPsi START *******************************
C **** PURPOSE                                           *
C Calculate the function psi(1+iy)                       *
C **** ARGS                                              *
C (y)                                                    *
C **** INPUT                                             *
C  y                                                     *
C **** OUTPUT                                            *
C cPsi psi(1+iy)                                         *
      FUNCTION cPsi(y)
      IMPLICIT REAL*8 (y)
      IMPLICIT COMPLEX*16(C)
      cZ=DCMPLX(21.D0,y)
      cZ1=1.D0/cZ
      cZ2=cZ1*cZ1
      cZ4=cZ2*cZ2
      cZ6=cZ4*cZ2
      cF=LOG(cZ)-cZ1/2.D0-cZ2/12.D0+cZ4/120.D0-cZ6/252.D0
        DO i=1,20
        cZ=cZ-1.D0
        cF=cF-1.D0/cZ
        ENDDO
      cPsi=cF
      RETURN
      END
C **** FUNCTION cPsi END *********************************


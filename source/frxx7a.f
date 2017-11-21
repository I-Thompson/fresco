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

!***********************************************************************
      SUBROUTINE MULTIP(FORMF,NF,N,LOCF,NK,FPT,NK1,FPT2,A,B,P,Q,RSP,
     &    FORML,FORMC,NLN,QNF,RIN,NNU,NC,QQ,KQ1,
     &    H,KFRAG,KCORE,KOPT,NLL,HPOT,NF0,PTYPE,
     &    ZF,ZC,ZT,MR,QSCALE,KPCORE,STREN,LAMBDA,CPSO,QCMAX,
     &    NBINS,JEX,BAND,PSIGN,EMID,DELTAE,BE,PCOUP)
C
	use parameters
	use trace
	use io
	use fresco1, only: cxwf,ccbins,sumkql,ompform
!$      use omp_lib
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16,ALLOCATABLE:: QERN(:,:,:)
      REAL*8 FORML(MAXNLR,MSP),PCOUP(NK1)
      REAL*8 BE(MSP,1),EMID(MSP),JEX(6,MXP,MXX)
      COMPLEX*16 FORMC(MAXNLC,MSP),STRENC(NK),MEK!,STRENK(NK1)
      COMPLEX*16 FORMF(MAXM,NF),VC,VF,VO,V,VCQ
      COMPLEX*16 VS(0:47),VCORE(NLN,10),VFRAG(NLN,2),VOPT(NLN,2)
      COMPLEX*16, EXTERNAL:: FFCI
      COMPLEX*16 MEKM(MSP,MSP,KQ1),MIA(MSP,MSP,KQ1)
      COMPLEX QSCALE(0:11)
      REAL*8 JPT,JPF,JCT,JCF,JPIOLD
      REAL*8 DBDE(KQ1),BEI(MSP,KQ1),DELTAE(MSP)
      INTEGER FPT(6,NK),FPT2(8,NK1),QNF(19,MSP),QQ,KQ,KQ1,KQ1M,QC1(8)
      INTEGER PTYPE(8,NF),PKIND,LAMBDA(NF),KSO,LSO
      INTEGER LA,QC,LAA,QCMAX,QLM,QL,NQC,IQC,PKINDD
      INTEGER BAND(2,MXP,MXX),PHASEFAC,PHASE,QQ1,PF,PT
      LOGICAL CPSO,NEWSET,JOPEN,DEFCORE
      CHARACTER*1 PSIGN(3)
      REAL*8 PLEG(NNU,KQ1),UK(NNU),WT(NNU),R2(NNU),RCR(NNU),HPOT(NF)
     &       ,WG(6),XG(6),AI(4),JS,JT,STREN(NF),RK,MEKR
      DATA NWW,XG(4),XG(5),XG(6),WG(4),WG(5),WG(6)/3,
     1   .2386191861D0,.6612093865D0,.9324695142D0,
     2   .4679139346D0,.3607615730D0,.1713244924D0/, PI /3.14159D0/
C
      MSO = 0
      if(CPSO) MSO=3
      QLM=QCMAX*(3+QCMAX)/2 ! NUMBER OF ARRAY ELEMENTS FOR Q MULTIPOLES
      MSOQL=MSO+QLM         ! ADD ON THE 3 NEEDED FOR SPIN-ORBIT
! note set qcmax to have max of 8, ie only can have core spin max of 4
       ALLOCATE(QERN(MAXN,KQ1,0:MSOQL))
      NW=2*NWW
      DO 1 J=1,NWW
      JJ=NW-J+1
      XG(J)=-XG(JJ)
  1   WG(J)=WG(JJ)
      RINTP = 1.0/RIN
      CUT = 1E-5
      KQ = KQ1 - 1
      KQ1M = 1
      IF(QQ.LT.0) KQ1M = -QQ + 1
      DO 10 J=1,NK
      STREN(LOCF+J) = 0.0
      STRENC(J) = 0.0
      HPOT(J+LOCF) = H
      DO 10 I=1,N
10    FORMF(I,J+LOCF) = 0.0
      AB2 = A*B*2
      PQ2 = P*Q*2
         THMAX = PI
           IC3 = NNU/NW
           NNT = IC3 * NW
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
               if(KQ>0) PLEG(K,2) = UK(K)
13             CONTINUE
14           C2 = C2 + CI
               DO 15 L=2,KQ
               DO 15 K=1,NNT
               WD = UK(K)
15         PLEG(K,L+1) = ((2*L-1)*WD*PLEG(K,L-1+1) -(L-1)*PLEG(K,L-2+1))
     &                       / DBLE(L)
      FMAX = 0.0
      DO 20 I=1,NLN
      VFRAG(I,:) = 0.0
      VCORE(I,:) = 0.0
      VOPT(I,:)  = 0.0
      IF(I.EQ.1) GO TO 20
      DO 19 IK=1,NK
        if(.not.cxwf) RP = ABS(FORML(I,FPT(1,IK)) * FORML(I,FPT(2,IK)))
        if(cxwf) RP = ABS(FORMC(I,FPT(1,IK)) * FORMC(I,FPT(2,IK)))
        RP = RP * ((I-1)*RINTP)**(QNF(9,FPT(1,IK))+QNF(9,FPT(2,IK))+2)
19    FMAX = MAX(FMAX, RP )
20    CONTINUE
      RR = FMAX * CUT
      IMIN = MIN(10,NLN)
      DO 21 I=NLN,IMIN,-1
      DO 21 IK=1,NK
        if(.not.cxwf) RP = ABS(FORML(I,FPT(1,IK)) * FORML(I,FPT(2,IK)))
        if(     cxwf) RP = ABS(FORMC(I,FPT(1,IK)) * FORMC(I,FPT(2,IK)))
        RP = RP * ((I-1)*RINTP)**(QNF(9,FPT(1,IK))+QNF(9,FPT(2,IK))+2)
        IF(RP.GT.RR) GO TO 22
21    CONTINUE
22    NLIM = I
!      IF(LISTCC.GE.1) WRITE(KO,225) CUT,FMAX,NLIM,(NLIM-1)*RINTP
!225   FORMAT(' Using cutoff of',1P,E9.1,' times maximum of',0P,F9.2,
!     &       ' gives NLIM =',I4,' at',F6.1,' fm.')
	NLIM = RSP/RINTP
        WRITE(KO,225) RSP,NLIM,(NLIM-1)*RINTP
225   FORMAT(' Using cutoff at',F9.1,' => NLIM =',I4,' at',F6.1,' fm.')
      HFRAG = 1.0/(MR*RIN)
      HCORE = HFRAG
      HOPT  = HFRAG
C
       QC1(:)=0
       NQC=0
      DO 24 JF=1,NF0
!      IF(LISTCC.GE.3) WRITE(KO,'(3x,7i9)') JF,PTYPE(1:6,JF)
      KP = PTYPE(1,JF)
      PKIND = PTYPE(2,JF)
      DEFCORE = KP==KCORE .AND. (PKIND==10 .OR. PKIND==11)
      IF(PTYPE(3,JF).NE.0 .AND. .not.DEFCORE)GOTO 24
      IF(PTYPE(4,JF)<0) GO TO 24  !SKIP UNDEFORMED POTENTIAL TO AVOID DOUBLE COUNTING
C                      (ONLY CENTRAL and PROJECTILE VSO FORCES CONSIDERED HERE)
C              (NOW ALSO INCLUDES DEFORMATION OF CORE)
      QC=0
      PKINDD=-1
      IF(DEFCORE)THEN
       PKINDD = PTYPE(5,JF) ! PKIND OF UNDEFORMED PARENT POT
       IF(PTYPE(3,JF)>0)THEN
        IF(QC<=QCMAX.AND.QCMAX>0)THEN
         QC = PTYPE(3,JF)
         NQC = NQC+1
         IF(NQC<=QCMAX)QC1(NQC) = QC
        ENDIF
       ENDIF
      ENDIF
	             LSO = 1
	IF(PKIND==3) LSO = 2
        IF(QC>0)     LSO = 2+QC

      IF(LISTCC.GE.3)WRITE(KO,'(/4i4,l4,4i4,6l4)',advance='no')
     &JF,KP,PKIND,PKINDD,DEFCORE,NC,KFRAG,KCORE,KOPT,
     &PKIND>3,(PKIND==0.or.PKIND>3),PKIND.NE.0,
     &PKINDD>3,(PKINDD==0.OR.PKINDD>3),PKINDD.NE.0
C
      IF(.NOT.DEFCORE)THEN
        IF(NC.EQ.0 .AND. PKIND>3 .OR.
     &   NC.EQ.1 .AND. (PKIND==0.or.PKIND>3) .OR.
C                         ONLY NUCLEAR, NO COUL
     &   NC.EQ.2 .AND. PKIND.NE.0 
C                         ONLY COULOMB, NO NUCLEAR
     &  ) GOTO 24
      ELSE !IF(DEFCORE)THEN
        IF(NC.EQ.0 .AND. PKINDD>3 .OR.
     &   NC.EQ.1 .AND. (PKINDD==0.OR.PKINDD>3) .OR.
C                         ONLY NUCLEAR, NO COUL
     &   NC.EQ.2 .AND. PKINDD.NE.0 
C                         ONLY COULOMB, NO NUCLEAR
     &  ) GOTO 24
      ENDIF
C
      IF(LISTCC.GE.3.and.((KP==KFRAG).or.(KP==KCORE).or.(KP==KOPT)))
     &WRITE(KO,'(a)',advance='no') '  breakup potential used'
C
      T1 = 1.0
      T2 = 1.0
      T3 = 1.0
      IF(PKIND.EQ.0) THEN
         T1 = ZF*ZT
         T2 = ZC*ZT
         T3 = T1 + T2
         ENDIF
      IF(PKIND.ne.0.and.KP.eq.KOPT.and.PTYPE(6,JF)>=2) THEN
	write(6,226) PKIND,KOPT,JF
226	format(//' ******* NUCLEAR PART',i4,' OF POTENTIAL  ',I3,
     X		 ' AT',i3,' IS **NOT** SUBTRACTED'/)
	T3 = 0.
       endif
      IF(KP.EQ.KFRAG) HFRAG = HPOT(JF)
      IF(KP.EQ.KCORE) HCORE = HPOT(JF)
      IF(KP.EQ.KOPT ) HOPT  = HPOT(JF)
      DO 23 I=1,NLN
         II = (I-1)*MR + 1
      IF(KP.EQ.KFRAG)  VFRAG(I,LSO) = VFRAG(I,LSO) + FORMF(II,JF) * T1
      IF(KP.EQ.KCORE)  VCORE(I,LSO) = VCORE(I,LSO) + FORMF(II,JF) * T2
23    IF(KP.EQ.KOPT)   VOPT(I,LSO)  = VOPT(I,LSO)  + FORMF(II,JF) * T3
24    CONTINUE
      IF(LISTCC.GE.3) WRITE(KO,*)
      IF(LISTCC.GE.3 .OR. HFRAG*HOPT*HCORE.LT.1E-20) THEN
         WRITE(KO,241) 'FRAGMENT',HFRAG,VFRAG(:,1)
         WRITE(KO,241) 'CORE'    ,HCORE,VCORE(:,1)
         WRITE(KO,241) 'OPTICAL' ,HOPT,VOPT(:,1)
	if(CPSO) then
         WRITE(KO,241) 'FRAG Vso',HFRAG,VFRAG(:,2)
         WRITE(KO,241) 'CORE Vso',HCORE,VCORE(:,2)
         WRITE(KO,241) 'OPT. Vso',HOPT,VOPT(:,2)
	endif
        DO IQC=1,NQC
         QC=QC1(IQC)         
         WRITE(KO,242) 'CORE Q=',QC,HCORE,VCORE(:,2+QC)
        ENDDO
241    FORMAT(' ',A8,' potential (step size ',F6.4,' fm) is',/,
     X  (1X,12F10.4) )
242    FORMAT(' ',A8,I2,' potential (step size ',F6.4,' fm) is',/,
     X  (1X,12F10.4) )
         WRITE(KO,*) A,B,P,Q
      ENDIF
C
      SINRF = 1.0/(MR * HFRAG)
      SINRC = 1.0/(MR * HCORE)
      SINRO = 1.0/(MR * HOPT)
C
      DO 30 IK=1,NK
30    PTYPE(3,LOCF+IK) = FPT(3,IK)
C
      INLIM=NLIM/10
      NN = (NLL-1)*MR+1
	write(6,*) ' FF to max ',(NN-1)*H
!***********************************************************************

      DO 270 J=1,NLIM
         R1 = (J-1) / RIN

      QERN = 0.0
        DO 250 I=1,NN,MR
            RR = (I-1) * H
            VO = FFCI(RR*SINRO,VOPT,NLN,0)
            if(CPSO) VO2 = FFCI(RR*SINRO,VOPT(1,2),NLN,0)
C          IF(RR.GE.(NLN-1)*RINTP - R1 * MAX(ABS(B),ABS(Q)))  GO TO 250
            ABR = (A*RR)**2 + (B*R1)**2
            PQR = (P*RR)**2 + (Q*R1)**2
            RRR  = R1*RR
            DO 231 K=1,NNT
               RCR2 =     PQR + PQ2*RRR*UK(K)
               R22  =     ABR + AB2*RRR*UK(K)
               RCR(K) = SQRT(ABS(RCR2)) * SINRC
231            R2(K) = SQRT(ABS(R22)) * SINRF
               DO 240 K=1,NNT
                   VF = FFCI(R2(K),VFRAG,NLN,0)
                   VC = FFCI(RCR(K),VCORE,NLN,0)
		  if(CPSO) then
                   VF2= FFCI(R2(K),VFRAG(1,2),NLN,0)
                   VC2= FFCI(RCR(K),VCORE(1,2),NLN,0)
		  endif
                   VS(0) = (VF + VC - VO) * WT(K)
	 	  if(CPSO) then
                   VS(1) =  - VO2 * WT(K)
                   VS(2) =  VF2 * WT(K)
                   VS(3) =  VC2 * WT(K)
		  endif
                  VS(4:MSOQL) = (0d0,0d0)
                  do IQC=1,NQC
                  QC=QC1(IQC)
                  VCQ = FFCI(RCR(K),VCORE(1,2+QC),NLN,0)
                  do LAA=0,QC
                   QL=QC*(QC+1)/2+LAA
                   if(RCR(K)/=0.) VS(MSO+QL) = VCQ/RCR(K)**QC * WT(K)
     &                             * (Q*R1)**(QC-LAA) * RR**LAA
                  enddo
                  enddo
!               IF(LISTCC.EQ.9) VS(0) = WT(K)
	    DO 236 KSO=0,MSOQL
            DO 236 L1=KQ1M,KQ1
             QERN(I,L1,KSO)=QERN(I,L1,KSO) + VS(KSO) * PLEG(K,L1)
236         CONTINUE
240        CONTINUE
250     CONTINUE

      if(.not.sumkql)then

       DO 260 IK=1,NK
           KNT = FPT(1,IK)
           KNF = FPT(2,IK)
           L1  = FPT(3,IK) + 1
           KSO = FPT(4,IK)
           QC  = FPT(5,IK)
           LA  = FPT(6,IK)
           QL=QC*(QC+1)/2+LA
           if(QL>0) QL=QL+MSO
              IF(L1.LT.KQ1M .OR. L1.GT.KQ1 .OR. QL*KSO/=0) GO TO 260
        RP = 0.5 * R1 ** (QNF(9,KNT)+1 + QNF(9,KNF)+1)  * RINTP
        if(.not.cxwf) V = RP * FORML(J,KNT) * FORML(J,KNF)
        if(     cxwf) V = RP * conjg(FORMC(J,KNT)) * FORMC(J,KNF)
!         STRENC(IK) = STRENC(IK) + 2.0*V
         V = V * SQRT(2.0 * (L1-1) + 1.0)
        DO 255 I=1,NN,MR
         FORMF(I,LOCF+IK) = FORMF(I,LOCF+IK) +   V * QERN(I,L1,KSO+QL)
255     CONTINUE
260    CONTINUE

      else !if(sumkql)then

       DO 261 IK1=1,NK1
           KNF = FPT2(2,IK1)   ! FROM (PRIMED   } IN frxx4.f ) ?????
           KNT = FPT2(1,IK1)   ! TO   (UNPRIMED }            ) ?????
           LA  = FPT2(3,IK1)   ! new multipole of summed formfactors
           L1  = FPT2(4,IK1)+1 ! K+1 multipole for legendre polys
           QL  = FPT2(5,IK1)   !  QL=QC*(QC+1)/2+LAA
           IK  = FPT2(6,IK1)   ! coupling no. in FPT to be summed into
           KSO = FPT2(7,IK1)   ! spin-orbit force
           if(QL>0) QL=QL+MSO
              IF(L1.LT.KQ1M .OR. L1.GT.KQ1 .OR. QL*KSO/=0) GO TO 261
        RP = 0.5 * R1 ** (QNF(9,KNT)+1 + QNF(9,KNF)+1)  * RINTP
        if(.not.cxwf) V = RP * FORML(J,KNT) * FORML(J,KNF)
        if(     cxwf) V = RP * conjg(FORMC(J,KNT)) * FORMC(J,KNF)
!         IF(I==1)STRENC(IK) = STRENC(IK) + 2.0*V
         V = V * SQRT(2.0 * (L1-1) + 1.0)
        DO 256 I=1,NN,MR
        FORMF(I,LOCF+IK)=FORMF(I,LOCF+IK)+V*QERN(I,L1,KSO+QL)*PCOUP(IK1)
256     CONTINUE ! I
261    CONTINUE ! IK1

      endif !sumkql

270   CONTINUE ! J
	timeformf = 0.
!$     timeformf = OMP_GET_WTIME()-timeformf0
       if(timeformf>0.) WRITE(KO,*) 'FORMFACTORS ALL DONE at ',timeformf
       CALL FLUSH(2)

!***********************************************************************
      DEALLOCATE(QERN)
!***********************************************************************

      DO 90 I=1,NN
         IF(MOD(I-1,MR).EQ.0) GO TO 90
         RR = (I-1)
         CALL SPLINT(RR/MR,NLL,II,M,AI)
         DO 80 J=1,M
         DO 80 IK=1,NK
           L1 = LOCF + IK
80       FORMF(I,L1) = FORMF(I,L1) + AI(J) * FORMF((II+J-2)*MR+1,L1)
90    CONTINUE
      IF(KPCORE.ge.10) THEN
         DO 95 IK=1,NK
            L1 = FPT(3,IK) + 1
            DO 95 I=1,N
95          FORMF(I,LOCF+IK) = FORMF(I,LOCF+IK) * QSCALE(L1-1)
!      ELSE IF(KPCORE.EQ.4.OR.KPCORE.EQ.5) THEN
!C                   ADD QSCALE * MULTIPOLE TO MONOPOLE, SET MULTIPOLE=0
!        DO 104 IK=1,NK
!        L1 = FPT(3,IK)+1
!        IF(L1.EQ.1) GO TO 104
!        DO 100 IK1=1,NK
!        IF(FPT(3,IK1).EQ.0 .AND. FPT(1,IK1).EQ.FPT(1,IK)
!     $                     .AND. FPT(2,IK1).EQ.FPT(2,IK)) GO TO 101
!100     CONTINUE
!        GO TO 104
!101     CONTINUE
!       WRITE(KO,*) 'ADDING  ',QSCALE(L1-1),
!     $   ' of multipole ',L1-1,' in column ',IK,' to column ',IK1
!        DO 102 I=1,N
!        FORMF(I,LOCF+IK1) = FORMF(I,LOCF+IK1) +
!     $           QSCALE(L1-1)*FORMF(I,LOCF+IK)
!102     FORMF(I,LOCF+IK)  = 0.0
!104     CONTINUE
       ENDIF
C
      STRENC(1:NK) = 0.0

      IF(NN.ge.N-MR-4) then
       DO 160 IK=1,NK
        LAMBDA(LOCF+IK) = FPT(3,IK)
        DO 160 I=N-MR-23,N-MR-4
         R = (I-1)*H
         STRENC(IK) = STRENC(IK) +
     X                   FORMF(I,LOCF+IK) * R ** (LAMBDA(LOCF+IK)+1)
160    CONTINUE       
       DO 170 IK=1,NK
        STRENC(IK) = STRENC(IK) * 0.05
        STREN(LOCF+IK) = STRENC(IK)    ! real part only!
170    CONTINUE
       DO 175 I=N-MR-4,N
        R = (I-1)*H
        DO 175 IK=1,NK
         FORMF(I,LOCF+IK) = STRENC(IK)/ R**(LAMBDA(LOCF+IK)+1)
175    CONTINUE
      ENDIF
!***********************************************************************
      IF(LISTCC.LE.0) go to 999
      K89=89
      K = MIN(NK,6)
      WRITE(KO,110) (FPT(1,IK),FPT(3:4,IK),FPT(2,IK),IK=1,K)
110   FORMAT(//' Single-Particle multipole form factors for overlaps'//
     &  ' Radius',6(4X,:,'<',I3,' /',2I2,' /',I3,'>'))
      L1=MIN(NN+MR,N)
      DO 120 I=1,L1,MR
120   WRITE(KO,130) (I-1)*H,(FORMF(I,LOCF+IK),IK=1,K)
      DO 121 I=L1+1,N
121   WRITE(KO,130) (I-1)*H,(FORMF(I,LOCF+IK),IK=1,K)
130   FORMAT(1X,F6.2,6(F11.5,F9.3))
      WRITE(KO,140) 'NORMS',(STREN(LOCF+IK),IK=1,K)
140   FORMAT('0',A5, F12.4,5(9X,F12.4))
134   continue
	if(TRENEG>=1) then
	call openif(89)
135	do IK=1,NK
        WRITE(K89,'(''#'',4i5)') IK+LOCF,N,MR,2
        WRITE(K89,141) IK,IK+LOCF,FPT(1,IK),FPT(3:4,IK),FPT(2,IK)
141     FORMAT('# Single-particle multipole form factors ',i4,' at',I5,
     &    ': <',I3,' /',2I2,' /',I3,'>')
        DO 143 I=1,N,MR
143     WRITE(K89,144) (I-1)*H,FORMF(I,LOCF+IK)
144   FORMAT(1X,F8.3,1p,2g13.5)
	WRITE(K89,*) '&'
	enddo
	endif
      call flush(K89)
      IF(sumkql)go to 999
      IF(NC.EQ.1.OR.NN.LT.N-MR-4) go to 999

	if(abs(ZT).lt.1e-3) then
	  write(KO,*) ' No Coulomb multipole strengths as ZT =',ZT
	  go to 999
	  endif
      WRITE(KO,140) 'M(EK)',(STREN(LOCF+IK) * (2*FPT(3,IK)+1)
     &            / (SQRT(4*PI) * COULCN * ZT ),IK=1,K)
         Z = 0.0
      WRITE(KO,*)
      WRITE(KO,*) ' Multipole transitions using pure single-particle',
     &       		' wave functions:'
	call rewop(47)
      MEKM(:,:,:)=(0d0,0d0)
      DO 180 IK=1,NK
      RK = FPT(3,IK)
      IF(FPT(3,IK).EQ.0.or.FPT(4,IK)>0) GO TO 180
           KNF = FPT(2,IK)
           KNT = FPT(1,IK)
             JS = QNF(11,KNF)*.5
             LS = QNF(9,KNF)
             JT = QNF(11,KNT)*.5
             LT = QNF(9,KNT)
             SS = QNF(10,KNF)*.5
C  The R factor is to transform from <L'||Ek||L> to <L'SJ'||Ek||LSJ> 
C       (more correctly, to  <LT SS JT || Ek || LS SS JS> )
             R = sqrt((2*JS+1)*(2*JT+1))
     x       	  * RACAH(SS,JS,LT+Z,RK,LS+Z,JT)
              RR = (2*RK+1) / (SQRT(4*PI) * COULCN * ZT)
             RRR=sqrt((2.*LS+1.)*(2.*LT+1.))*wig3j(LT+Z,RK,LS+Z,Z,Z,Z)
     x       		* (-1)**abs(LS-LT)
      MEK = STRENC(IK) * R  * RR * RRR
      BEK = abs(MEK)**2 / (2 * JS + 1.)
      WRITE(KO,190) IK,FPT(3,IK),JS,KNF,JT,KNT,BEK,FPT(3,IK),MEK
      MEKR = MEK  ! real part only
      if(MEKR.gt.-9.9.and.MEKR.lt.99) then
       WRITE(47,1302) KNT,KNF,FPT(3,IK),MEKR
        else
        WRITE(47,1303) KNT,KNF,FPT(3,IK),MEKR
       endif
180   CONTINUE
C
1302     FORMAT(4X,3I4,F8.4)
1303     FORMAT(4X,3I4,F8.1)
190   FORMAT(' ',I5,': B(E',I1,',',F4.1,I4,' -> ',F4.1,I4,') =',F13.6,
     X       ' from lsj M(E',I1,') =',2F9.4)
C

999   RETURN
      END
!***********************************************************************
      SUBROUTINE INTER(ICTO,ICFROM,IC,KIND,QQ,IP2,IP3,IP4,IP5,
     &                 BETAR,BETAI,KLT,
     &                 NIB,HPOT, REV,NLL,LOCF,COUPLE,ICOR,ICOM,CP,ITER,
     &                 XA,XB,XP,XQ,RSP,FFREAL,FILE,
     &        MASS,HP,JEX,ENEX,QVAL,NEX,BAND,COPY,ITC,CPOT,
     &        PTYPE,NAME,FORMF,FORML,FORMC,BE,D0,BPROJ,POTCAP,
     &        QNF,AFRAC,FPT,NKP,GPT,NFI,USED,NPRIOR,NPOST,OPN,
     &        M,HCM,EPS,NAXQRN,NF,NF0,NSP,NLCN,RIN,WID,CENTR,INFILE,
     &        MATRIX,MEK,NIX,SPINTR, STREN,LAMBDA,HASO,CPSO,
     X        LTRANS,REW14,MTMIN,LCALL,MCALL,MCALLS,
     X        NBINS,PSIGN,EMID,DELTAE)
	use io
	use parameters
	use kcom
	use drier
	use trace
	use factorials
	use searchpar
	use fresco1, only: cxwf,ccbins,sumccbins,sumkql,nrbases
      IMPLICIT REAL*8(A-H,O-Z)
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      REAL*8 MASS(4,MXP+1),HP(MXP),JEX(6,MXP,MXX),HPOT(MLOC)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1),EMID(MSP),DELTAE(MSP)
      INTEGER NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     &  PKIND,ITC(MXP,MXX),CPOT(MXP,MXX),PTYPE(8,MLOC)
      CHARACTER*8 NAME(2,MXP+1),COMM(4)
      CHARACTER*1 PSIGN(3)
C
C    FORM FACTORS AND THEIR PARAMETERS
C    ---------------------------------
      PARAMETER(MMX=2)
      REAL*8 FORML(MAXNLR,MSP,2),FORMR(MAXM),FFR4
      COMPLEX*16 FORMC(MAXNLC,MSP,2),FFC4
      COMPLEX*16 FORMF(MAXM,MLOC),FFCI,FFC,F8,FRMC(MMX*MAXM),VSP
      COMPLEX QSCALE(0:11)
      REAL*8 BE(MSP,3),D0(MSP,2),
     &       FORM(MMX*MAXM,2),BPROJ(MSP,MXP),STREN(MLOC)
      INTEGER QNF(19,MSP),LAMBDA(MLOC),POTCAP(MLOC,2)
C
C    COUPLINGS
C    ---------
      REAL*8 AFRAC(MXPEX,MXPEX,2,MSP),BEN(4,MXPEX,MXPEX),MEK(MPAIR)
     X       ,SPINTR(2,MPAIR),ICORE,ICOREPP,KICORE,KICOREP,JCOM,JCOMP
      INTEGER CP,ICTO,ICFROM,KIND,QQ,MATRIX(6,MPAIR),fails,
     &        QC,QCMIN,QCMAX,QCM,LA,INFILE(MAXCPL),
     &        LOCF,FPT(6,MAXQRN),NKP(2),FIL,USEDLOW(MSP),USEDHIGH(MSP),
     &        ICOR(MAXCPL,2),ICOM(MAXCPL,2),GPT(2,MAXQRN),FILE,NFI(3)
      INTEGER FPT2(8,MAXQRN2)
      REAL*8 P(MAXQRN2)
      LOGICAL REV,USED(2,MSP),COUPLE,NPRIOR,NPOST,OPN(2),FFREAL,
     & 		HASO(kpmax),CPSO,INCOUL
      CHARACTER*12 REM(6),REOR(16),KSCAL(4),GRID(3)
      CHARACTER*9 LNL(2)
      CHARACTER*8 POTL(5),POPR(4),REALCM(2),RIC(3)
C
C    PARTIAL WAVES AND THEIR PARAMETERS
C    ----------------------------------
      COMPLEX*16 VOPT(MAXNLN),VPOT(MAXM,2),FNC(NLL,NLO),ASCALE,S
C
C    ARRAYS FOR NON-LOCAL COUPLING FORM FACTORS
C    ------------------------------------------
      REAL*8 DNL(MAXNLO),JN,JNP
      INTEGER,SAVE:: JFCF
      LOGICAL LTRANS(MAXQRN),LCALL,MCALL,REW14,keep,MCALLS
      CHARACTER*70 FILC*2
   	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN
	namelist/inel/ ib,ia,k,no,kp,a
	namelist/cfp/ in,ib,ia,kn,a,keep
	namelist/scale/ qscale
C
      LOGICAL FRAC,FAIL3
C
      DATA POPR /'POST','PRIOR','PRI-POST','VCORE'/,
     &  LNL /'  LOCAL  ','NON-LOCAL' /,
     &  RIC / '  REAL  ','  IMAG. ',' COMPLEX' /,
     &  KSCAL / 'as given.','as CHEX2','jlmP targ','jlmP proj' /,
     X  GRID / 'printing.','calling FFNL','FNLSET call' /,
     & POTL / 'READ CMP','READ RL:','COUL+NUC','NUCLEAR ','COULOMB '/,
     &  REALCM /'  REAL  ','COMPLEX '/
     &  REM /'CLUSTER FOLD','NONO+COMPLEX','COMPLEX',
     &       'NONE','FULL REAL','NONO + REAL'/
     & ,REOR / '       FULL ','NO DIAG','ONLY DIAG   ','ONLY FROM GS',
     &         'DIAG+To/FrGS','DIAG+To/FrBS','No bin J/=J''',
     &       3*'            ',
     &         'SCALE FULL','NONE,SC','ONLY DIAG,SC','ONLY F.GS,SC',
     &         'DIAG+GS,SC'  ,'DG+ToFrBS,SC' /
!    &         'Q=0 SUMMED &',     'Q=0 SUMMD,NO' /
C
      FRAC(X) = ABS(X-NINT(X)).GT.1D-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      Z = 0D0
      INCOUL=.true.
!	write(KO,*) 'HASO(1:10) =',HASO(1:10)
!	write(KO,*) 'NF =',NF
	 DO 1 IN=1,4
	 DO 1 ITC1=1,MXPEX
	 DO 1 ITCO=1,MXPEX
1        BEN(IN,ITC1,ITCO) = 0d0
232   FORMAT(6E12.4)
      LOCF = NF
      fails = 0
      FPT(:,:)=0
      FPT2(:,:)=0
      POTCAP(:,:) = 0
C
      IF(KIND.NE.1) GO TO 100
C
C  General projectile/target spin transfer, external form factors.
C
	KF = INFILE(CP)
      IF(IP3.GE.0)  WRITE(KO,5) KF,LNL(QQ+1),RIC(IP2+1),
     x  KSCAL(IP3+1),BETAR,BETAI
5     FORMAT(' Therefore from file ',i3,' read ',A9,' form factor of ',
     x A8,' elements, and use ',A12,' with Re,Im scalings of',2f9.4)
      IF(IP3.LT.0) WRITE(KO,6) GRID(-IP3)
6     FORMAT(' Generate radial grids for ',A12)
C
       LOCF = NF
       NKP(1) = NIX + 1
10    if(IP3==1) then
	READ(KF,119,END=50) FNP,HNP,FSCALE,LTR,PTR,TTR,IB,IA,COMM
	RFS = HNP
	NP = nint(FNP)
       else
        READ(KF,12,END=50) NP,HNP,RFS,FSCALE,LTR,PTR,TTR,IB,IA,COMM
       endif
119   FORMAT(3F8.4,I4,2F4.0,2I4,4A8)
12    FORMAT(I4,3F8.4,I4,2F4.0,2I4,4A8)
      IF(ABS(FSCALE).LT.1E-20) GO TO 50
      CALL CHECK(NP,MMX*MAXM,7)
      IF(IB.EQ.0) IB = 1
      IF(IA.EQ.0) IA = 1
        if(IA>NEX(ICFROM)) then
                write(KO,*) 'ERROR: THERE IS NO STATE ',IA,'
     x          : only ',NEX(ICFROM),' in partition ',ICFROM
                stop
                endif
        if(IB>NEX(ICTO)) then
                write(KO,*) 'ERROR: THERE IS NO STATE ',IB,
     x          ': only ',NEX(ICTO),' in partition ',ICTO
                stop
                endif

      SN = JEX(1,ICTO,ABS(IB))
      SNP= JEX(1,ICFROM,IA)
C
      IF(QQ.NE.0) THEN
C                             Non-local form factor
       ASCALE = FSCALE
       NX = NLO
       CALL CHECK(NLO,MAXNLO,9)
       HTARG = HP(ICTO) * MR
       NTARG = NLL
       HNP = HTARG
       NP = NLL
      ELSE
C                         Local form factor
       ASCALE = FSCALE
      IF(IP3.EQ.1) ASCALE = FSCALE * R4PI
     X       * SQRT((2.*SNP+1.)/((2.*PTR+1.)*(2.*SN+1.)))
      IF(IP3.EQ.0) ASCALE = FSCALE * R4PI
        NX = 1
        HTARG = HP(ICTO)
        NTARG = M + 1
        NF = NF + 1
        CALL CHECK(NF,MLOC,24)
        HPOT(NF) = HP(ICTO)
       ENDIF
C
      IF(IP3.GE.0) WRITE(KO,20) NP,HNP,RFS,COMM
20    FORMAT(/' Read ',I4,' point form factor at h =',F6.3,
     x             ' fm from',f6.3,':',4A8)
      WRITE(KO,21) FSCALE,LTR,PTR,TTR,IB,IA
21    FORMAT( ' Scaled by',F9.4,' for L-transfer',I3,
     X',  projectile transfer',F4.1, ', and target transfer =',F5.1,
     X   ' to excited pair ',I3,' from pair',I3)
      SN = JEX(1,ICTO,ABS(IB))
      SNP= JEX(1,ICFROM,IA)
	if(FAIL3(SN,SNP,PTR).or.FRAC(SN+SNP+PTR)) then
	write(ko,211) 'PROJ',SN,SNP,PTR,.not.FAIL3(SN,SNP,PTR)
211	format(/' ****ERROR: NO ',a4,' COUPLING OF ',2f5.1,' with',f5.1,
     x    ' as Delta =',L2,' or integer-spin error')
	endif
      JN = JEX(2,ICTO,ABS(IB))
      JNP= JEX(2,ICFROM,IA)
	if(FAIL3(JN,JNP,TTR).or.FRAC(JN+JNP+TTR)) then
	write(ko,211) 'TARG',JN,JNP,TTR,.not.FAIL3(JN,JNP,TTR)
	endif
C
      IF(HNP.LT.1E-5) HNP = HTARG
      DO 30 JX=1,NX
       IF(IP3.GE.0) THEN
         FORM(:,1:2) = 0.
         IF(IP2.eq.0) READ(KF,*) (FORM(I,1),I=1,NP)
         IF(IP2.eq.1) READ(KF,*) (FORM(I,2),I=1,NP)
         IF(IP2.eq.2) READ(KF,*) (FORM(I,1),FORM(I,2),I=1,NP)
         DO 23 I=1,NP
          R1 = 0.0
          R2 = 0.0
          IF(mod(IP2,2).eq.0) R1 = FORM(I,1)
          IF(IP2       .ge.1) R2 = FORM(I,2)
          FRMC(I) = CMPLX(BETAR*R1,BETAI*R2)
23       CONTINUE
      ELSE
C            NON-LOCAL:
        DNL(JX) = (JX - NLC - 1) * MLT * HP(ICFROM)
      ENDIF
         DO 25 I=1,NTARG
          R = (I-1) * HTARG
          T2 = (R-RFS) / HNP
          F8 = FFCI(T2,FRMC,NP,LTR)
          IF(QQ.EQ.0) FORMF(I,NF) = F8 * ASCALE
          IF(QQ.EQ.1.AND.IP3.GE.0) FNC(I,JX) = F8 * ASCALE
          IF(QQ.EQ.1.AND.IP3.LT.0) FNC(I,JX) = R + (0.,1.)*(R + DNL(JX))
25       CONTINUE
30      CONTINUE
      IF(IP3.EQ.-1) THEN
         WRITE(KO,35) NLL,NLO
35       FORMAT(2I4)
         DO 36 I=1,NLL
36       WRITE(KO,232) (FNC(I,J),J=1,NLO)
         GO TO 999
       ENDIF
      IF(IP3.EQ.-2) CALL FFNL(FNC,NLL,NLO,ASCALE,LTR,PTR,TTR,
     X                             ICTO,ICFROM,IB,IA)
      IF(IP3.EQ.-3) CALL FNLSET(FNC,NLL,NLO,HTARG,DNL,ASCALE,
     X                             LTR,PTR,TTR,ICTO,ICFROM,IB,IA)
C
      NIX = NIX + 1
      CALL CHECK(NIX,MPAIR,16)
      MATRIX(1,NIX) = ABS(IB)
      MATRIX(2,NIX) = IA
      MATRIX(3,NIX) = LTR
      MATRIX(4,NIX) = NF
      MEK(NIX) = 0.0
      SPINTR(1,NIX) = PTR
      SPINTR(2,NIX) = TTR
C
      IF(QQ.EQ.1) THEN
         LREC66 = LREC66 + 1
         MATRIX(4,NIX) = -LREC66
         WRITE(66) FNC
        ENDIF
C
       if(TRENEG>=1) then
        call openif(89)
        WRITE(89,'(''#'',4i5)') NF,N,MR,3
        WRITE(89,141) NIX,NF,IA,LTR,PTR,TTR,IB
141     FORMAT('# General multipole form factor ',i4,' at',I5,
     &    ': <',I3,' /',I2,2f4.0,' /',I3,'>')
        DO 143 I=1,N,MR
143     WRITE(89,144) (I-1)*HTARG,FORMF(I,NF)
144   FORMAT(1X,F8.3,1p,2g13.5)
        WRITE(89,*) '&'
        endif

      IF(LISTCC.GE.2) THEN
       IF(QQ.EQ.0) THEN
		write(KO,*)
              WRITE(KO,*) ' Form factor at ',NF,' is (after scaling) '
              WRITE(KO,132) ((I-1)*HP(ICTO),FORMF(I,NF),I=1,N,MR)
132           FORMAT(4(1X,0P,F6.2,':',1P,2E11.3))
        ELSE
           IF(LISTCC.GE.3) THEN
            DO 45 I=1,NLL
             WRITE(KO,40) (I-1)*HTARG
40           FORMAT(' At R-to =',F8.4,', the non-local range is')
             WRITE(KO,44) (FNC(I,J),J=1,NLO)
44           FORMAT(5(1X,2F9.4))
45          CONTINUE
          ENDIF
         CALL DISPLR(FNC,NLL,NLO,NLL,SCALR)
           SCALI = SCALR * 1E-6
         CALL DISPLI(FNC,NLL,NLO,NLL,SCALI)
          DO 47 I=1,NLL
           VOPT(I) = 0.
           DO 47 J=1,NLO
 47        VOPT(I) = VOPT(I) + FNC(I,J) * MLT * HP(ICFROM)
          WRITE(KO,44)
          WRITE(KO,*) ' Low-energy local equivalent potential every ',
     X                   HTARG
          WRITE(KO,44) (VOPT(I),I=1,NLL)
        ENDIF
       ENDIF
C
      IF(IB.GT.0) GO TO 10
C
 50   NKP(2) = NIX
      GO TO 999
C
100   IF(KIND.EQ.8.AND.IP2.EQ.0) GO TO 876
!       IF(LISTCC.GE.1) WRITE(KO,*) ' LOCF = ',LOCF
C
C    READ IN COEFFICIENTS OF FRACTIONAL PARENTAGE FOR TRANSFER STATES
C    ................................................................
C                             (DO THIS FOR KINDS 2, 3,4, 5,6,7,8)
8735  continue
!      READ(KI,8739) IN,IB,IA,KN,A
      IN=0; IB=0; IA=0; KN=0; A=0
      READ(KI,nml=cfp) 
8739  FORMAT(4X,4I4,F8.4)
      IF(IN.EQ.0) GO TO 876
      I = IN
      IN = ABS(IN)
	ipa = 0
!		Adjust any search parameter!
	numafrac = numafrac + 1
	do ip=1,nvars
	if(srch_kind(ip)==2.and.numafrac==srch_nafrac(ip)) then
	    if(abs(srch_value(ip)-nul)>.001) then
              A = srch_value(ip)
	    else
	     srch_value(ip) =  A
	    endif
	  ipa = ip
	endif
	enddo
      IF(KN.GT.NSP.OR.QNF(8,KN).EQ.0) GO TO 8755
      IF(IB.GT.NEX(ICOM(CP,IN)).OR.IA.GT.NEX(ICOR(CP,IN))) GO TO 8755
      ITC1= ITC(ICOM(CP,IN),IB)
      ITCO= ITC(ICOR(CP,IN),IA)
       AFRAC(ITC1,ITCO,IN,KN) = A
      DO 874 KN2=KN,NSP-1
      IF(QNF(1,KN2+1) .NE. KN) GO TO 8745
       IF(.not.sumccbins) AFRAC(ITC1,ITCO,IN,KN2+1) = A
874   CONTINUE
8745  WRITE(KO,875) I,IB,IA,KN,A,NAME(IN,ICOM(CP,IN)),IB,
     &                       NAME(IN,ICOR(CP,IN)),IA,KN,A,BE(KN,1),KN2
875   FORMAT(/' data =',4I3,F8.4, ',  so  < [',A8,' #',I3,'] / [',
     &A8,' #',I3,']  * Eigen-form #',I3,'>  =',F9.4,12X,'BE =',F7.3,I6)
	if(ipa>0) 
     &WRITE(srch_afrac_overlap(ipa),87501) NAME(IN,ICOM(CP,IN)),IB,
     &                       NAME(IN,ICOR(CP,IN)),IA,KN
87501 FORMAT('< [',A8,' #',I2,'] / [',A8,' #',I2,']  * KN#',I3,'>')
         BEN(IN,ITC1,ITCO) = BEN(IN,ITC1,ITCO) + BE(KN,1) * A*A
         BEN(2+IN,ITC1,ITCO) = BEN(2+IN,ITC1,ITCO) +          A*A
      IF(KIND.EQ.5.OR.KIND.EQ.6) THEN
      DO 8751 KN2=KN,NSP
      IF(.not.sumccbins) A = AFRAC(ITC1,ITCO,IN,KN2)
8751   IF(ABS(A).GT.1E-5.AND.QNF(1,KN2).EQ.KN) BPROJ(KN2,ICOM(CP,2))
     &       =   BE(KN2,1) - QVAL(ICOR(CP,1)) + QVAL(ICOM(CP,1))
     &                  + ENEX(2,ICOR(CP,1),IB) - ENEX(2,ICOM(CP,1),IA)
      ENDIF
         J = 1
         IF(QNF(16,KN).GE.ICHAR('N')) J = -1
         IF(QNF(7,KN).ge.6) J = (-1)**QNF(13,KN)
      IF((QNF(7,KN).EQ.0 .OR. QNF(7,KN).EQ.6)  .AND. (
     &        FAIL3(JEX(IN,ICOR(CP,IN),IA),
     X         QNF(11,KN)*0.5D0,JEX(IN,ICOM(CP,IN),IB))
     &.OR.BAND(IN,ICOR(CP,IN),IA)*BAND(IN,ICOM(CP,IN),IB)
     &    *(-1)**QNF(9,KN)*J .LT. 0))   then
	WRITE(KO,8753) KN,
     &      JEX(IN,ICOR(CP,IN),IA),QNF(11,KN)*.5,JEX(IN,ICOM(CP,IN),IB)
8753  FORMAT(///' **** WARNING *** -- FORM FACTOR',I3,' GIVES  NO COUP'
     &,'LING AS PARITY ERROR, OR D(',F4.1,2(',',F4.1),') ',
     &  ' = (JCORE,JV,JCOMP)  FAILS  ****'//)
       fails = fails+1
       endif
      J = QNF(6,KN)
      IF((IN.NE.QNF(6,KN).OR.ICOR(CP,IN).NE.QNF(2,KN).OR.ICOM(CP,IN).NE.
     & QNF(4,KN)).and.(NAME(IN,ICOM(CP,IN)).ne.NAME(J ,QNF(4,KN)).or.
     &                 NAME(IN,ICOR(CP,IN)).ne.NAME(J,QNF(2,KN))) ) 
     &  WRITE(KO,8754) KN,NAME(J ,QNF(4,KN)),NAME(J,QNF(2,KN))
     &                            ,IN,J,QNF(4,KN) ,        QNF(2,KN)
8754  FORMAT(///' **** WARNING *** -- FORM FACTOR',I3,' WAS ORIGINALLY'
     &,' DEFINED FOR A  ',A8,' / ',A8,' OVERLAP ****',4I4//)
	call flush(KO)
8755  IF(I.GT.0) GO TO 8735
C
C                     HAVE NOW TO CALCULATE TRANSFER OR INELASTIC
C                     FORM FACTORS FOR THE VARIOUS MULTIPOLES.
876    IF(ABS(KIND).LE.4) then
         IN1 = MOD(ABS(KIND)-1,2)+1
         if(KIND==2) IN1=2
         ZC = MASS(2+IN1,ICOR(CP,IN1))
         ZF = MASS(2+IN1,ICOM(CP,IN1)) - ZC
         ZT = MASS(2+(3-IN1),ICOM(CP,IN1))
	endif
	USEDLOW(:) = MXX+1
	USEDHIGH(:) = 0
      if(.not.sumccbins)then
       
       DO 8763 IN=1,2
        NEXA = NEX(ICOR(CP,IN))
        NEXB = NEX(ICOM(CP,IN))
	IF(ABS(KIND)>4) IN1=IN
       DO 8763 KN=1,NSP
          USED(IN,KN) = .FALSE.
       DO 8763 IA=1,NEXA
       DO 8763 IB=1,NEXB
        A = AFRAC(ITC(ICOM(CP,IN),IB),ITC(ICOR(CP,IN),IA),IN,KN)
         IF(LISTCC.GE.5) then
          WRITE(KO,8761) IN,KN,IA,IB,A
8761            format('IN,KN,IA,IB,A = ',4i4,f10.5)
          WRITE(KO,*) 'IM,IR=',ITC(ICOM(CP,IN),IB),ITC(ICOR(CP,IN),IA)
          endif
       USED(IN,KN) = ABS(A).GT.EPS.AND.QNF(8,KN).NE.0 .OR. USED(IN,KN)
       IF(KIND.EQ.5.or.KIND.eq.2) USED(1,KN) = KN.EQ.1
       if(IN==IN1.and.ABS(A)>EPS.and.QNF(8,KN)/=0) then
	 USEDLOW(KN) = min(IB,USEDLOW(KN))
	 USEDHIGH(KN) = max(IB,USEDHIGH(KN))
       endif
8763   CONTINUE
       
      endif
        IF(LISTCC.GE.1) then
           WRITE(KO,*) ' LOCF = ',LOCF
!           WRITE(KO,*) ' cxwf,ccbins,sumkql =',cxwf,ccbins,sumkql
!           WRITE(KO,*) ' USEDHIGH(*) = ',USEDHIGH(:)
!           WRITE(KO,*) ' USEDLOW(*) = ',USEDLOW(:)
!            WRITE(KO,*) ' USED(1,*) = ',USED(1,1:NSP)
!            WRITE(KO,*) ' USED(2,*) = ',USED(2,1:NSP)
         endif
!           WRITE(KO,*) ' USEDHIGH(*) = ',USEDHIGH(:)
!           WRITE(KO,*) ' USEDLOW(*)  = ',USEDLOW(:)

      NK = 0
      NK1= 0
      NG = 0
      QCM= 0
      LAM= 0
      IF(IP4>=0) QCM=IP4
      IF(IP5>=0) LAM=IP5
      LAM1= 0
         IN1 = 1
         IN2 = 2
         IF(ABS(KIND).LE.4.and.KIND/=2) then
             IN1 = MOD(ABS(KIND)-1,2)+1
             IN2 = IN1
 		endif
         LCALL = .FALSE.
         MCALL = .FALSE.
      DO 878 KNP=1,NSP
      DO 878 KNT=1,NSP
      IF(.NOT.(USED(IN1,KNP).AND.USED(IN2,KNT))) GOTO 878
      IF(KIND.NE.5.and.KIND.ne.2 .AND. (QNF(16,KNT).NE.QNF(16,KNP)
     &       .OR.QNF(12,KNT) .NE. QNF(12,KNP)
     &       .OR.ABS(QNF(10,KNT)-QNF(10,KNP)).GT.1E-5) ) GO TO 878
      IF(QNF(12,KNP).NE.0 .AND. (QNF(13,KNP).NE.QNF(13,KNT) .OR.
     &                           QNF(14,KNP).NE.QNF(14,KNT))) GOTO 878
         LN = QNF(9,KNT)
         LNP= QNF(9,KNP)
         JN = QNF(11,KNT) * 0.5
         JNP= QNF(11,KNP) * 0.5
         SN = QNF(10,KNT) * 0.5
!***********************************************************************
      IF(KIND.GE.5.and.KIND.le.8) THEN
C                           TRANSFER FORM FACTORS SUMMED INTO 'IG'
      !if(KIND==7.and.ICFROM==2.AND.ICTO==1.and.QNF(5,KNP)/=1)goto 878 ! temp removal of back couplings to continuum states
      IF(QNF(12,KNP).LE.1) THEN
         NG = NG + 1
         IF(NG.GT.MAXQRN) GO TO 878
         GPT(1,NG) = KNP
         GPT(2,NG) = KNT
         LTRANS(NG) = LN + LNP .LT. MTMIN .or. KIND.le.6.or.KIND==2
         LCALL = LCALL .OR. LTRANS(NG)
         MCALL = MCALL .OR. .NOT.LTRANS(NG)
         KLT = MAX(KLT, 2*LNP+1)
       ENDIF
         NK = NK + 1
         IF(NK.GT.MAXQRN) GO TO 878
         FPT(1,NK) = KNP
         FPT(2,NK) = KNT
             J = MAX(QNF(12,KNP)-1,0)
            DO 8781 IG=1,NG
8781        IF(GPT(1,IG).EQ.KNP-J .AND. GPT(2,IG).EQ.KNT-J)
     &         FPT(3,NK) = IG
          IF(LISTCC.GE.0)
     &    WRITE(KO,87811) KNT,KNP,NK,NG,FPT(3,NK),
     X            FPT(3,NK)+LOCF,LTRANS(FPT(3,NK)),LCALL,MCALL
87811    FORMAT(' KNT =',i4,', KNP =',i4,': NK,NG,loc=',3i4,i4,3L3)

!***********************************************************************
      ELSE IF(KIND.GE.3.and.KIND.le.4.and..not.sumkql) THEN
C        KIND = 3 & 4 :    S.P. INELASTIC FORM FACTORS OF MULTIPOLE 'IK'
      KFRAG = NINT(BETAR)
      KCORE = NINT(BETAI)
         IF(KCORE.EQ.0) KCORE = CPOT(ICOR(CP,IN1),1)
      KOPT = CPOT(ICOM(CP,IN1),1)
	CPSO = .false.
	if(KFRAG>0.and.KCORE>0.and.KOPT>0.and.KIND==3) 
     x	CPSO = HASO(KFRAG) .or. HASO(KCORE) !   .or. HASO(KOPT) 
!	write(KO,*) HASO(KFRAG),HASO(KCORE),HASO(KOPT),' so CPSO=',CPSO
!       ! if no core excited states or no def core 3-body coup.
      QCMIN=0
      QCMAX=0

      DO 8782 IK=ABS(LN-LNP), MIN(LNP+LN,ABS(QQ))
         IF(QQ.LT.0 .AND. IK.NE.ABS(QQ)) GOTO 8782

                        IF(MOD(IK+LNP+LN,2).NE.0) GOTO 8789
                        IF(FAIL3(IK+Z,JN,JNP)) GOTO 8789
	if(USEDHIGH(KNP)>USEDLOW(KNT))
     x	WRITE(KO,'(4i7,'' tried '')')NK+1,KNP,KNT,IK
!                       Impose order to use symmetry when wfs real:
!        IF(.not.cxwf  .and.USEDHIGH(KNP)>USEDLOW(KNT)
!    x       .and..not.CPSO) GO TO 8782
         NK = NK + 1
         WRITE(KO,'(4i7)') NK,KNP,KNT,IK
         IF(NK<=MAXQRN) then
            FPT(1,NK)   = KNP
            FPT(2,NK)   = KNT
            FPT(3,NK)   = IK
            FPT(4:6,NK) = 0
	 endif
	 if(CPSO)then
	  do KSO=1,3
	   if(KSO==0 .or. KSO==1.and.HASO(KOPT) .or.
     X        KSO==2.and.HASO(KFRAG) .or. KSO==3.and.HASO(KCORE)) then
            NK = NK + 1
         if(.not.ccbins.or.LISTCC>=1)WRITE(KO,'(5i7)')NK,KNP,KNT,IK,KSO
            IF(NK<=MAXQRN) then
              FPT(2,NK)   = KNP
              FPT(1,NK)   = KNT
              FPT(3,NK)   = IK
              FPT(4,NK)   = KSO
              FPT(5:6,NK) = 0
	     endif
	   endif
	  enddo
         endif
8787     continue !LA
8789     continue !QC
8782   CONTINUE !IK
	if(CPSO) WRITE(KO,*) ' FOLDED SPIN-ORBIT FORCES'
!***********************************************************************
   ! ccbins treatment removed in this version
!***********************************************************************
      ELSE IF(KIND==2) THEN
!@@      
!      IF(QNF(12,KNP).LE.1) THEN
!         NG = NG + 1
!         IF(NG.GT.MAXQRN) GO TO 878
!         GPT(1,NG) = KNP
!         GPT(2,NG) = KNT
!       ENDIF


      DO 87821 IK=1, ABS(QQ)
         IF(QQ.LT.0 .AND. IK.NE.ABS(QQ)) GO TO 87821
         NK = NK + 1
	 if(LISTCC>0) 
     &    WRITE(KO,*) 'Need:',NK,' uses bs',KNP,KNT,' q=',IK
         IF(NK<=MAXQRN) then
            FPT(1,NK) = KNP
            FPT(2,NK) = KNT
            FPT(3,NK) = IK
            FPT(4,NK) = 0
	 endif


          IF(LISTCC.GE.3)
     &    WRITE(KO,87811) KNT,KNP,NK,NG,FPT(3,NK),FPT(3,NK)+LOCF
87821   CONTINUE
!@@      
      ENDIF
!***********************************************************************
878   CONTINUE !KNP,KNT
         CALL FLUSH(KO)
         CALL CHECK(MAX(NK,NG),MAXQRN,13)
         CALL CHECK(NK1,MAXQRN2,13)
         CALL CHECK(MAX(NK,NG),MAXQRN,13)
         NKP(1) = NK
         NKP(2) = NG
         NAXQRN = MAX(NAXQRN,NK)


      IF(NK.GT.0)                   GO TO 8764
         WRITE(KO,8783)
8783     FORMAT('0NO COUPLED PAIRS OF FORM FACTORS ARE FOUND!!!'/)
      DO 8784 IN=1,2
8784  WRITE(KO,8785) NAME(IN,ICFROM),NSP,(USED(IN,KN),KN=1,NSP)
8785  FORMAT('0Use with ',A8,' of the',I3,' form factors are'/(1X,60L2))
         COUPLE = .FALSE.
	 fails = fails+1
         go to 999
8764  IF(KIND.GE.7)  GO TO 877
      IF(KIND.EQ.5) WRITE(KO,8765) BETAR,BETAI
      IF(KIND.EQ.6) WRITE(KO,87651) BETAI
8765  FORMAT(/5X,'So ZERO-RANGE TRANSFER, with D0 =',F10.3,' and',
     &  ' finite-range effective radius',F9.4, ' fm.'/)
87651 FORMAT(/5X,'So ZERO-RANGE TRANSFER, with D0 derived from project',
     &'ile states.'/,5X,'If bound, derive finite-range effective radii',
     &'from the "D" values of those states;',
     & 'If unbound or IP3=1, use E.R. =',F9.4)
     
     
      IF(KIND.eq.2) then      
!@@
      WRITE(KO,88400) QQ, IP2, IP3
88400 FORMAT(/5X,'Electromagnetic coupling : IP1 =',I2,', IP2='
     &,I2,', IP3=',I2)
      if(mod(IP2,4)/=1) WRITE(KO,88401) BETAR, BETAI
88401 format(5x,' Particle g factor =',F10.3,' and'
     &,' target g factor',F9.4)
!
      IF(QQ.GT.0 .AND. mod(IP2,4).EQ.0 ) WRITE(KO,88402) QQ
88402 FORMAT(5X,' All Electric and Magnetic multipoles up to ', I2)
      IF(QQ.GT.0 .AND. mod(IP2,4).EQ.1 ) WRITE(KO,88403) QQ
88403 FORMAT(5X,' All Electric multipoles up to ', I2)
      IF(QQ.GT.0 .AND. mod(IP2,4).EQ.2 ) WRITE(KO,88404) QQ
88404 FORMAT(5X,' All Magnetic multipoles up to ', I2)     
!      
      IF(QQ.LT.0 .AND. mod(IP2,4).EQ.0 ) WRITE(KO,88405) abs(QQ),abs(QQ)
88405 FORMAT(5X,' Calculates E',I2,' and M',I2, ' multipoles')
      IF(QQ.LT.2 .AND. mod(IP2,4).EQ.1 ) WRITE(KO,88406) abs(QQ)
88406 FORMAT(5X,' Calculates E',I2,' multipole only')
      IF(QQ.LT.2 .AND. mod(IP2,4).EQ.2 ) WRITE(KO,88407) abs(QQ)
88407 FORMAT(5X,' Calculates M',I2,' multipole only')   
	if(IP2>=4) write(KO,*) '     INCLUDING the Siegert remnant'
! 
!!      IF(IP3.EQ.0 .AND. mod(IP2,4).EQ.0 ) WRITE(KO,88408) 
!!88408 FORMAT(5X,' NO Intrinsic convection couplings')
!!      IF(IP3.EQ.1 .AND. mod(IP2,4).EQ.1 ) WRITE(KO,88409) 
!!88409 FORMAT(5X,' Intrinsic convection coupling for projectile only')
!!      IF(IP3.EQ.2 .AND. mod(IP2,4).EQ.2 ) WRITE(KO,88410) 
!!88410 FORMAT(5X,' Intrinsic convection coupling for target only')
!!      IF(IP3.EQ.3 .AND. mod(IP2,4).EQ.0 ) WRITE(KO,88411) 
!!88411 FORMAT(5X,' Intrinsic convection coupling for projectile 
!!     Xand target')
!!      IF(IP3.EQ.4 .AND. mod(IP2,4).EQ.1 ) WRITE(KO,88412) 
!!88412 FORMAT(5X,' Magnetization term for electric components ON')     
!@@ 

C
C    ALLOCATE BOUND STATE FORM FACTORS FOR PHOTONUCLEAR COUPLINGS
C    ................................................................
      LOCF = NF
      NF = NF + NK	! for bs
      if(IP2>=4) then
	 NF = NF + NK   ! for bs vertex function
	 !  add in number of potential shapes for scattering in particle partition
	 npots=0
	 do JF=1,NF0
	  I = 0
	  do IA=1,NEX(ICFROM)
	  if(PTYPE(1,JF)==CPOT(ICFROM,IA)) then
	     I = NK
!	     write(KO,*) 'Form ',JF,' copied to ',NF,'+',NK
	     endif
	  enddo
          if(I>0) then
            if(NF+I<=MLOC) POTCAP(NF+1:NF+NK,1) = JF
!            if(NF+I<=MLOC) POTCAP(NF+1:NF+NK,2) = (/ (IN,IN=1,NK) /)
            NF = NF + NK
            npots = npots+1
          endif
         enddo
         endif

      CALL CHECK(NF,MLOC,24)
      FORMF(1:N,LOCF+1:NF) = 0.0
       if(IP2<4) npots=-1
      POTCAP(LOCF+1:NF,2) = npots  ! in all places!
        
      DO IN=1,NK
         KN = FPT(2,IN)
         KNP= FPT(1,IN)
         IK = FPT(3,IN)
         IG = IN
         IF(QNF(6,KN).EQ.INH) HDN = HP(QNF(2,KN))
         IF(QNF(6,KN).NE.INH) HDN = HP(QNF(4,KN))
C         HDN = step size required for channel form factor,
C                NOT that already used to store bound states (HCM).

C     IF(LISTCC>1)
      WRITE(KO,8767) KN,KNP,IN,IG+LOCF,HDN
!      READ(8,REC=KN) (FORMR(I),I=1,N)
!      DO  I=1,N
!         R = (I-1)*HDN
!         T2 = R / HCM
!       FORMF(I,IG+LOCF) = FORMF(I,IG+LOCF) +  FFR4(T2,FORMR,N) 
!      ENDDO
      DO  I=1,N
         R = (I-1)*HDN
         T1 = R / RINTP
          RP =  R ** QNF(9,KN)
         if(.not.cxwf) S = FFR4(T1,FORML(1,KN,1),NLN)
         if(     cxwf) S = FFC4(T1,FORMC(1,KN,1),NLN)
        FORMF(I,IG+LOCF) = FORMF(I,IG+LOCF) +  S * RP 	     ! bound state wf
       if(IP2>=4) then
         VSP = 0.0
         IF(.not.cxwf.and.S.NE.(0.,0.)) VSP=FFR4(T1,FORML(1,KN,2),NLN)
         IF(     cxwf.and.S.NE.(0.,0.)) VSP=FFC4(T1,FORMC(1,KN,2),NLN)
        FORMF(I,IG+NK+LOCF) = FORMF(I,IG+NK+LOCF) +  VSP * RP   ! vertex function for bs
	do IIF=1,npots
         JF = LOCF+(IIF+1)*NK + IN
        FORMF(I,JF) = FORMF(I,JF) + S*RP*FORMF(I,POTCAP(JF,1))   ! bs wf * scattering potl 
	enddo
       endif
       ENDDO

        IF(LISTCC.ge.1) then
         do IIF=-1,npots
            JF = LOCF+(IIF+1)*NK + IN
              WRITE(KO,*) ' Form factor at ',JF,' from ',POTCAP(JF,1),
     x                 ' is :'
              WRITE(KO, 132) ((I-1)*HDN,FORMF(I,JF),I=1,N,MR)
	   if(TRENEG>=1) then
	    call openif(89)
              WRITE(89,'(''#'',4i5)') JF,N,MR,2
              WRITE(89,*) ' Gamma ff at ',JF,' from ',POTCAP(JF,1)
              DO I=1,N,MR
               WRITE(89,144) (I-1)*HDN,FORMF(I,JF)
	       ENDDO
	      WRITE(89,*) '&'
	   endif
          enddo
            endif
      ENDDO
	write(KO,*)

      go to 999
      endif
      IF(KIND.NE.5.AND.KIND.NE.6) GO TO 901
C
C    CALCULATE LOCAL FORM FACTORS FOR ZR TRANSFERS (PERHAPS WITH LEA)
C    ................................................................
      LOCF = NF
      NF = NF + NG
      CALL CHECK(NF,MLOC,24)
       DO 8766 IG=LOCF+1,NF
       DO 8766 I=1,N
8766   FORMF(I,IG) = 0.0
      DO 8769 IN=1,NK
         KN = FPT(2,IN)
         KNP= FPT(1,IN)
         IG = FPT(3,IN)
         IF(QNF(6,KN).EQ.INH) HDN = HP(QNF(2,KN))
         IF(QNF(6,KN).NE.INH) HDN = HP(QNF(4,KN))
C         HDN = step size required for channel form factor,
C                NOT that already used to store bound states (HCM).
      KOPT = CPOT(ICOM(CP,1),1)
      KCORE= CPOT(ICOR(CP,1),1)
      DO 87661 I=1,N
      VPOT(I,1) = 0.0
87661  VPOT(I,2) = 0.0
      NC = 0
      DO 87668 JF=1,NF0
      KP = PTYPE(1,JF)
      IF(PTYPE(3,JF).NE.0) GO TO 87668
C                                   (ONLY SCALAR FORCES CONSIDERED HERE)
      PKIND = PTYPE(2,JF)
      IF(NC.EQ.0 .AND. PKIND.GT.2 .OR.
     &   NC.EQ.1 .AND. (PKIND.NE.1.AND.PKIND.NE.2) .OR.
C                         ONLY NUCLEAR, NO COUL
     &   NC.EQ.2 .AND. PKIND.NE.0 )  GOTO 87668
C                         ONLY COULOMB, NO NUCLEAR
      T1 = 1.0
      T2 = 1.0
      IF(PKIND.EQ.0) THEN
         T1  = MASS(2+1,ICOM(CP,1)) * MASS(2+2,ICOM(CP,1))
         T2  = MASS(2+1,ICOR(CP,1)) * MASS(2+2,ICOR(CP,1))
         ENDIF
      DO 87665 I=1,N
      IF(KP.EQ.KOPT)   VPOT(I,1) = VPOT(I,1) + FORMF(I,JF) * T1
87665 IF(KP.EQ.KCORE)  VPOT(I,2) = VPOT(I,2) + FORMF(I,JF) * T2
87668  CONTINUE
      DM =  MASS(1,ICOM(CP,1))-MASS(1,ICOR(CP,1))
      RM = DM * MASS(1,ICOR(CP,1)) / MASS(1,ICOM(CP,1))
      IF(KIND.EQ.5) THEN
             R1 = BETAR
             R0 = BPROJ(KN,ICOM(CP,2))
             T = BETAI**2
             R2 = R1/(1.0 - T * R0*RM*FMSCAL)
      ELSE
             R1 = D0(KNP,1)
             R2 = D0(KNP,2)
             R0 = BE(KNP,1)
             T = MAX(1-R1/R2,Z) /(R0*RM*FMSCAL)
             IF(BE(KNP,1).LT.0.0.OR.IP3.EQ.1) T = BETAI**2
      ENDIF
C     IF(T.GT.1E-5)
      WRITE(KO,8767) KN,KNP,IN,IG+LOCF,R0,R1,R2,DM,RM,SQRT(ABS(T)),HDN
8767  FORMAT(1X,4I3,10F10.3)
      if(.not.cxwf) READ(8,REC=KN) (FORMR(I),I=1,N)
      if(     cxwf) READ(8,REC=KN) (FRMC(I),I=1,N)
      DO 8768 I=1,N
         R = (I-1)*HDN
         T2 = R / HCM
         if(.not.cxwf) F8 = FFR4(T2,FORMR,N)
         if(     cxwf) F8 = FFC4(T2,FRMC,N)
         S = 1.0
      IF(T.GT.1E-5) THEN
         T1 = R / RINTP
         if(.not.cxwf) S = FFR4(T1,FORML(1,KN,1),NLN)
         if(     cxwf) S = FFC4(T1,FORMC(1,KN,1),NLN)
         VSP = 0.0
         IF(.not.cxwf.and.S.NE.(0.,0.)) VSP=FFR4(T1,FORML(1,KN,2),NLN)/S
         IF(     cxwf.and.S.NE.(0.,0.)) VSP=FFC4(T1,FORMC(1,KN,2),NLN)/S
         S = 1.0 + T * RM * FMSCAL * (VPOT(I,2)+VSP-VPOT(I,1)+R0)
         ENDIF
       FORMF(I,IG+LOCF) = FORMF(I,IG+LOCF) + R1 * F8 * S
8768   CONTINUE
        IF(LISTCC.LT.1) GO TO 8769
              WRITE(KO,*) ' Form factor at ',IG+LOCF,' is :'
              WRITE(KO, 132) ((I-1)*HDN,FORMF(I,IG+LOCF),I=1,N,MR)
8769  CONTINUE
      IF(KIND.EQ.6) GO TO 8853
      GO TO 999
C
C               SO NOW HAVE NON-LOCAL TRANSFER FORM FACTOR.
C               ..........................................
877   NC = nint(BETAR)
      BETAR = 0.0
      IC7 = MOD(ABS(QQ),2)
	if(QQ==2) IC7=2 ! for prior-post  March 2012
      IF(QQ.LE.-3 .AND. IP2.EQ.0 ) IP2 = 1
      NONO = 0
      IF(KIND.NE.8) GO TO 8788
      IP2 = 0
      NONO = 1
      NPRIOR = NPRIOR .OR. IC7.GT.0 .OR. REV
      NPOST  = NPOST  .OR. IC7.EQ.0 .OR. REV
8788  ICSUB = ICTO
      IF(IC7   .GT.0) ICSUB = ICFROM
      KSUB = CPOT(ICSUB,1)
      XMUT = MASS(2,ICTO) / MASS(2,ICFROM)
      XMUP = MASS(1,ICTO) / MASS(1,ICFROM)
      IF(XMUP.LT.1.) THEN
C                                STRIPPING FROM PROJECTILE : (D,P)
      ICOREP = ICTO
      ICORET = ICFROM
      ICV = IC7
      XMUT = 1.0/XMUT
      XLAM = 1. - XMUT * XMUP
      XA = 1./XLAM
      XB = - XMUP/XLAM
      XP = + XMUT/XLAM
      XQ= - XA
C--------------------    SIGN OF JACOBIAN CORRECTED 28TH APRIL, 1988.
      XJAC = +1.0/XLAM**3
      ELSE
      ICOREP = ICFROM
C                     HERE,      STRIPPING TO PROJECTILE : (P,D)
      ICORET = ICTO
      ICV = 1 - IC7
      XMUP = 1.0/XMUP
      XLAM = 1. - XMUT * XMUP
      XB = 1./XLAM
      XA = - XMUP/XLAM
      XQ= + XMUT/XLAM
      XP = - XB
      XJAC = +1.0/XLAM**3
      ENDIF
      IF(IP3.EQ.0) IP3 = CPOT(ICOREP,1)
      DO 8831 II=1,NLN
      VPOT(II,1) = 0.0
8831  VPOT(II,2) = 0.0
      HSUB = 0.
      HCOR = 0.
      IF(IP2.NE.0) THEN
      NC = 0
      if(NC.lt.0.or.NC.gt.2) NC = 0
      if(IP2==-3)then
! IP2=-3 for valid cluster folding potential label and PRIOR only
       if((JFCF<=1.or.JFCF>=MAXQRN).and.IC7==0) IP2=-1
      endif
      if(LISTCC>=3.and.IP2==-3)write(KO,*)'JFCF=',JFCF
      DO 8838 JF=1,NF0
      KP = PTYPE(1,JF)
      IF(PTYPE(3,JF).NE.0) GO TO 8838
C                                   (ONLY SCALAR FORCES CONSIDERED HERE)
      PKIND = PTYPE(2,JF)
      IF(NC.EQ.0 .AND. PKIND.GT.2 .OR.
     &   NC.EQ.1 .AND. (PKIND.NE.1.and.PKIND.NE.2) .OR.
C                         ONLY NUCLEAR, NO COUL
     &   NC.EQ.2 .AND. PKIND.NE.0 )  GOTO 8838
C                         ONLY COULOMB, NO NUCLEAR
      T1 = 1.0
      T2 = 1.0
      IF(PKIND.EQ.0) THEN
         T1  = MASS(2+1,ICSUB) * MASS(2+2,ICSUB)
         T2  = MASS(2+1,ICOREP) * MASS(2+2,ICORET)
	 if(.not.INCOUL) then
	  if(abs(T1)+abs(T2)>0.) then
	   write(KO,*) ' ***** COULOMB PART OMITTED IN REMNANT ***** '
	  endif
	  T1=0; T2=0
	  endif
         ENDIF
      IF(KP.EQ.KSUB) HSUB = HPOT(JF)
      IF(KP.EQ.IP3) HCOR = HPOT(JF)
      DO 8835 II=1,NLN
        I = (II-1) * MR + 1
      IF(KP.EQ.KSUB)    VPOT(II,1) = VPOT(II,1) + FORMF(I,JF) * T1
! IP2=-3 => use Cluster Folding core+valence potential
!           from single particle excitations
! note have to add to VPOT as VOPT subtracted in CF potential
! only valid for PRIOR ( PKIND=0 so only add once )
! set back to IP2=-1 later to treat as normal complex remnant
      IF(KP.EQ.KSUB.and.IP2==-3.and.PKIND==0)
     >                  VPOT(II,1)= VPOT(II,1) + FORMF(I,JFCF)
8835  IF(KP.EQ.IP3)  VPOT(II,2) = VPOT(II,2) + FORMF(I,JF) * T2
8838  CONTINUE
      ENDIF
      FFREAL = IP2.GE.0 .and..not.cxwf
c     FFREAL = .false.       				TEST ONLY!
      if(cxwf.and.REV) then
	WRITE(KO,8850)
8850	format(//'  Complex particle states AND reversible transfer ',
     x           'couplings requested!',/'  This combination is not ',
     x           'implemented correctly: please enter one or both ',
     x           'of the one-way couplings'/)
	if(ccbins) stop
      endif
        RINTO = HP(ICTO) * MR
        RRC = HCOR * MR
        RRS = HSUB * MR
        IF(RRC.LT.1E-5) RRC = RINTO
        IF(RRS.LT.1E-5) RRS = RINTO
        HNL   =(HP(ICFROM) * MLT)/ NLT
        CENTRE= HP(ICFROM) * MLT * NLCN
        HFHT = HP(ICFROM) / HP(ICTO)
        IC71 = IC7 + 1
        IF(QQ.LE.-3) IC71 = 4
      IF(NONO.EQ.0) then
	 WRITE(KO,884) QQ,POPR(IC71),IP2,REM(IP2+4), NK, NG
884   FORMAT(/5X,'So FINITE-RANGE TRANSFER : IP1 =',I2,' (',A8,'), IP2='
     &,I2,' =',A12,',  for',I4,' pairs of bound states summed into',I4,
     & ' pairs')
      IF(QQ.GE.-2.and.QQ<2) 
     &   WRITE(KO,8841) NAME(1+ICV,ICV*ICORET+(1-ICV)*ICOREP)
8841  FORMAT(5X,' The main coupling potential is that binding the tra',
     &'nsferred particle to the ',A8,' core.'/)
      IF(QQ.LE.-3) WRITE(KO,8842) NAME(1,ICOREP),NAME(2,ICORET),IP3
8842  FORMAT(5X,' Using ONLY the core-core (',A8,' - ',A8,') potential',
     X',  No.',I4/)
      IF(IP2/=0.and.NC/=0) WRITE(KO,8843) NC,POTL(NC+2)
8843  FORMAT(5X,' Using NC =',I2,':',A8,' parts of remnant potentials'/)
	endif
      IF(QQ.LE.-1) WRITE(KO,8844)
8844  FORMAT(5X,' With angular quadrature down from theta = 180 deg.'/)
      IF(NONO.EQ.1) WRITE(KO,8845) IC7,POPR(IC7+1)
8845  FORMAT(/5X,'So NON-ORTHOGONALITY SUPPLEMENT for an  IC7=',I1,' (',
     &       A4,')  coupling'/)
      IF(NONO.EQ.1) GO TO 8857
      IF(IP2.NE.0.AND.QQ.GE.-2.and.QQ.le.1)THEN
      if(IP2/=-3)then
      WRITE(KO,8851)NAME(1,ICOREP),NAME(2,ICORET),IP3,KSUB
8851  FORMAT(5X,'     Using between ',A8,' & ',A8,' the potential',
     & ' No.',I4,', and subtracting optical potential No.',I3/)
      else
      WRITE(KO,8849)NAME(1,ICOREP),NAME(2,ICORET),IP3,JFCF
8849  FORMAT(5X,'     Using between ',A8,' & ',A8,' the potential',
     & ' No.',I4,', and subtracting cluster folding potential No.',I3/)
      endif
      ENDIF
      IF(LCALL.AND.MCALL)
     X        WRITE(KO,8852) (GPT(1,I),GPT(2,I),LTRANS(I),I=1,NG)
8852  FORMAT(3(' Between forms',I3,' &',I3,', L-trans =',L2,:,'  '))
	MCALLS = MCALLS .or. MCALL
C
      IF(IP2==-3)IP2=-1 ! set back to complex remnant now cf pot added
C
      IF(LISTCC.GE.4) THEN
          WRITE(KO,1001) XMUT,XMUP,XLAM,XA,XB,XP,XQ,XJAC
1001  	FORMAT(1X,8F10.5)
        WRITE(KO,241) 'CORE'    ,HCOR*MR,VPOT(1:NLN,2)
        WRITE(KO,241) 'OPTICAL' ,HSUB*MR,VPOT(1:NLN,1)
241    FORMAT(' ',A8,' potential (step size ',F6.4,' fm) is',/,
     X  (1X,12F10.4) )
      ENDIF
!	write(48,*) 'COUPLE,ITER,MINL,MAXL =',COUPLE,ITER,MINL,MAXL
8857    continue
!	call flush(48)
      IF(.NOT.COUPLE.OR.(ITER.EQ.0.and.nrbases==0) .OR. MINL.GT.MAXL)
     x              GO TO 999
      IF(     FFREAL) NQ = 1
      IF(.NOT.FFREAL) NQ = 2
           CALL CHECK(NLO,MAXNLO,9)
           CALL CHECK(MAXL+1,MAXMUL,11)
           CALL CHECK(NNU,MAXNNU,12)
      IF(MCALL) THEN
C                    SOME M-TRANSFERS
        IF(.NOT.REW14) THEN
          IF(MACH==8.or.MACH==7.or.MACH==6.or.MACH==3.or.MACH==4) THEN
               OPEN(14,FILE=TMPN(1:lnbl(TMPN))//'14',
     X                 ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
             ELSE
               OPEN(14,ACCESS='SEQUENTIAL',STATUS='SCRATCH',
     X                 FORM='UNFORMATTED')
             ENDIF
               REWIND 14
          ENDIF
        WRITE(14) (VPOT(I,1),I=1,NLN),(VPOT(I,2),I=1,NLN),
     X             RRS,RRC,XJAC
        REW14 = .TRUE.
        ENDIF
      IF(.NOT.LCALL) GO TO 8853
C                    SOME L-TRANSFERS
             FIL = 13 - 2*NQ
           IF(.NOT.OPN(NQ)) THEN
	         inquire(iolength=lreal8) T1
		 write(KO,*) '  A real*8 is ',lreal8,' bytes'
		 write(KO,*) '  NLL,NLN,NLO =',NLL,NLN,NLO
                 I = NFI(FIL-8)+NG*(MAXL-MINL+1)
                 WRITE(KO,8858) FIL,NLL*NLO*NQ,'REAL*8',
     X	                 I,real(I)*NLL*NLO*NQ*lreal8*1E-6
 8858           FORMAT('0File ',I3,' needs RL =',I8,' ',A9,' numbers',
     X		       ', and NR =',I6,'. i.e.',F9.2,' Megabytes')
                 IF(MACH.eq.8.or.MACH==3.or.MACH==4) THEN
                 FILC = '0'//char(ichar('0') + FIL)
C-SHARED           OPEN(FIL,FILE=TMP(1:lnbl(TMP))//'fort.'//FILC,
                   OPEN(FIL,FILE=TMPN(1:lnbl(TMPN))//FILC,
     X                ACCESS='DIRECT',RECL=NQ*lreal8*NLO*(NLN-1))
C            -----------------------------SHARED FILE
C-SHARED          SETIOMODE(FIL,3)
                  ELSE if(MACH.eq.6.or.MACH.eq.7) THEN
               OPEN(FIL,ACCESS='DIRECT',RECL=NQ*lreal8*NLO*(NLN-1),
     X                FILE='fort.nq'//char(nq+ichar('0')))
	write(KO,*) '  File ',FIL,' RECL =',NQ*lreal8*NLO*(NLN-1)
                  ELSE
               OPEN(FIL,ACCESS='DIRECT',RECL=NQ*lreal8*NLO*(NLN-1),
     X                STATUS='SCRATCH')
	write(KO,*) '  File ',FIL,' RECL =',NQ*lreal8*NLO*(NLN-1)
                  ENDIF
                 OPN(NQ) = .TRUE.
               ENDIF
           FILE = NFI(FIL-8)
		NIBL=MAXL-MINL+1
                NIB = MAXL+1
      IF(LISTCC.GE.4) THEN
      DO IN=1,NK
         KN = FPT(2,IN)
	 write(KO,*) ' KN=',KN
       if(.not.cxwf) then
        WRITE(KO,241) 'wftn'  ,RINTP,FORML(1:10,KN,1)
        WRITE(KO,241) 'bind'  ,RINTP,FORML(1:10,KN,2)
	else
        WRITE(KO,241) 'cwftn'  ,RINTP,real(FORMC(1:10,KN,1))
        WRITE(KO,241) 'cbind'  ,RINTP,real(FORMC(1:10,KN,2))
	endif
	enddo
	endif
	 call flush(6)
      CALL QERNEL(NFI(FIL-8),NLN,NLM,NLO,MAXL+1,NK,XA,XB,
     & XP,XQ,FPT(1,1),NG,IC7,ICV,QQ,IP2,NNU,MINL+1,RINTO,
     & EPC,HNL,NLT,VPOT(1,1),VPOT(1,2),RRS,RRC,NLL,NIBL,cxwf,
C                  VFOLD     VCORE
     & FORML,FORML(1,1,2),FORMC,FORMC(1,1,2),NLN,RIN,NLC,WID,CENTR,
     & CENTRE,NONO,HFHT,NQ.EQ.1,XJAC,LTRANS)
C
      IF(NIBL.LT.MAXL-MINL+1)
     x WRITE(KO,*) 'Multipoles calculated in blocks of ',NIBL,
     x'  To speed up, need to increase MAXV '
C
8853  ICOREP = ICOR(CP,1)
      ICORET = ICOR(CP,2)
      NA = NEX(ICOREP)
      NAT= NEX(ICORET)
      DO 8855 IACP=1,NA
      DO 8855 IACT=1,NAT
      DO 8855 IAP=1,NAT
           T=BEN(3,ITC(ICORET,IAP),ITC(ICOREP,IACP))
           IF(ABS(T).LE.EPS) GO TO 8855
        BENP=BEN(1,ITC(ICORET,IAP),ITC(ICOREP,IACP))/T
      DO 8854 IAT=1,NA
         IF(COPY(1,ICORET,IAP,1)+COPY(1,ICOREP,IACP,1) +
     &      COPY(2,ICOREP,IAT,1)+COPY(2,ICORET,IACT,1).GT.0) GO TO 8854
           R=BEN(4,ITC(ICOREP,IAT),ITC(ICORET,IACT))
           IF(ABS(R).LE.EPS) GO TO 8854
        BENT=BEN(2,ITC(ICOREP,IAT),ITC(ICORET,IACT))/R
      E =-QVAL(ICOREP) + BENT + ENEX(1,ICOREP,IACP)+ENEX(2,ICOREP,IAT)
     & -(-QVAL(ICORET) + BENP + ENEX(1,ICORET,IAP)+ENEX(2,ICORET,IACT))
      IF(ABS(E).LT.0.001) GO TO 8854
      WRITE(KO,8856) NAME(1,ICORET),IAP,NAME(1,ICOREP),IACP,BENP,
     &           NAME(2,ICOREP),IAT,NAME(2,ICORET),IACT,BENT,E
8854  CONTINUE
8855  CONTINUE
8856  FORMAT('0Binding Energies: ',
     & 2('(',A8,' #',I2,' / ',A8,' #',I2,')(BE =',F7.3,')-'),
     & ' DISCREPANCY =',F8.4)
      GO TO 999
C
901   XB = MASS(IN1,ICOR(CP,IN1)) / MASS(IN1,ICOM(CP,IN1))
        XP = +1.0
      IF(IN1.EQ.1) THEN
        XA = 1.0
        XQ = XB - 1.
      ELSE
        XA = - 1.0
        XQ = +1 - XB
      ENDIF
C                           THESE ARE TO CALCULATE THE MATRIX ELEMENTS
C                           OF THE MULTIPOLES OF THE SINGLE-PARTICLE
C                           EXCITING POTENTIAL   V(SP)+V(CORE)-V(OPT)
C                           ..........................................
      IF(IP3.GE.10) THEN
C                          read in scaling factors for multipoles
!         READ(KI,232) (QSCALE(I),I=MAX(0,-QQ),ABS(QQ))
	  qscale(:)=0.
         READ(KI,nml=scale) 
         DO 950 I=MAX(0,-QQ),ABS(QQ)
950      WRITE(KO,955) I,QSCALE(I)
955      FORMAT('0 The MULTIPOLE ',I2,' forms are all MULTIPLIED by ',
     X          'the complex number',F10.5,' +',F10.5,' i')
        ENDIF
      LOCF = NF
      NF = NF + NK
      CALL CHECK(NF,MLOC,24)
           CALL CHECK(NNU,MAXNNU,12)
           CALL CHECK(ABS(QQ)+1,MAXMUL,11)
      KFRAG = NINT(BETAR)
         IF(KFRAG.EQ.0) WRITE(KO,*) 'NO S.P. EXCITING POTENTIAL KNOWN!!'
C        IF(KFRAG.EQ.0) STOP 'KFRAG'
         IF(KFRAG.EQ.0) CALL ABEND(16)
      KCORE = NINT(BETAI)
         IF(KCORE.EQ.0) KCORE = CPOT(ICOR(CP,IN1),1)
      KOPT = CPOT(ICOM(CP,IN1),1)
      call flush(KO)
      IF(COUPLE)THEN
       CALL MULTIP(FORMF,NF,N,LOCF,NK,FPT(1,1),NK1,FPT2(1,1),XA,XB,
     &            XP,XQ,RSP,FORML,FORMC,NLN,QNF,RIN,NNU,
     &   IP2,QQ,ABS(QQ)+1, HP(ICOM(CP,IN1)),KFRAG,KCORE,KOPT,NLL,HPOT,
     &   NF0,PTYPE,ZF,ZC,ZT,MR,QSCALE,IP3,STREN,LAMBDA,CPSO,QCM,
     &   NBINS,JEX,BAND,PSIGN,EMID,DELTAE,BE,P)
       JFCF=LOCF+1 ! define label for cluster folding potential
       IF(sumkql)THEN
        IF(QQ>=0)THEN
          QQ=MIN(LAM,LAM1)
        ELSE
         WRITE(KO,*)'ERROR: SUMMED COUPLED CHANNELS FORMFACTORS'
         WRITE(KO,*)' AND INDIVIDUAL MULTIPOLES NOT IMPLEMENTED'
         fails=fails+100
        ENDIF
       ENDIF !sumkql
      ENDIF !COUPLE
C
      IF(.not.sumkql.or.IP4==0)THEN
      WRITE(KO,910) NAME(IN1,ICOM(CP,IN1)),NAME(IN1,ICOR(CP,IN1)),
     &      QQ,POTL(IP2+3),IP2,KFRAG,KCORE,KOPT,
     &      REOR(IP3+1),IP3
      ELSE
      WRITE(KO,911) NAME(IN1,ICOM(CP,IN1)),NAME(IN1,ICOR(CP,IN1)),
     &      QQ,POTL(IP2+3),IP2,KFRAG,KCORE,KOPT,
     &      REOR(IP3+1),IP3
      ENDIF
910   FORMAT(/5X,'So INELASTIC EXCITATIONS OF SINGLE PARTICLES in the',
     & 1X,A8,' nucleus outside the ',A8,' core'//
     & 5X,'by Q =',I3,' multipoles of the ',A8,' parts (',I1,') of ',
     & 'V-sp(',I2,') + V-CC(',I2,') - V-opt(',I2,') with ',A12,
     & ' REORIENTATION (',I2,')'/)
911   FORMAT(/5X,'So INELASTIC EXCITATIONS OF SINGLE PARTICLES in the',
     & 1X,A8,' nucleus outside the ',A8,' core'//
     & 5X,'by LA=',I3,' multipoles of the ',A8,' parts (',I1,') of ',
     & 'V-sp(',I2,') + V-CC(',I2,') - V-opt(',I2,') with ',A12,
     & ' REORIENTATION (',I2,')'/)
C
999   if(fails>0) then
	write(KO,*) 'STOPPING, as some coupling failures'
	write(KO,*) '   See WARNING messages above'
!        call ABEND(8)
	endif
      RETURN
      END
      SUBROUTINE PARTEX(NAME,MASS,NEX,QVAL,RMASS,HCM,HP,COPY,EXTRA,
     X     GIVEXS,JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X     NCHPMAX,PWFLAG,FCWFN,IOFAM)
	use io
	use parameters
	use searchpar, only: final
      IMPLICIT REAL*8(A-H,O-Z)
C
      INTEGER PEL,NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),ICOP(2),
     X  ITC(MXP,MXX),BHEAD(2,2,0:9 ),CPOT(MXP,MXX),NCHPMAX(MXP),
     X  IOFAM(2,MXP,MXX)
      LOGICAL GIVEXS(MXP),EXTRA(2,MXP,MXX),FLIP,PWFLAG(MXP),FCWFN,EX(2)
      REAL*8 MASS(4,MXP+1),RMASS(MXP),HP(MXP),JEX(6,MXP,MXX)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1)
      CHARACTER PSIGN(3)
      LOGICAL PWF
      CHARACTER*8 NAME(2,MXP+1)
C
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      IT = 0
      PEL = 1
      FLIP = .FALSE.
      MXX10 = MXX + 10
      DO 20 IN=1,2
      DO 20 IC=1,MXP
      DO 20 IA=1,MXX
      DO 20 I=1,2
   20 COPY(IN,IC,IA,I) = 0
      IC = 1
   30 continue
!	READ(KI,1130) NAME(1,IC),MASS(1,IC),MASS(2+1,IC),NEX(IC),
!     X   PWF,NAME(2,IC),MASS(2,IC),MASS(2+2,IC),QVAL(IC)
	call readpt(KI,NAME(1,IC),MASS(1,IC),MASS(2+1,IC),NEX(IC),
     X   PWF,NAME(2,IC),MASS(2,IC),MASS(2+2,IC),QVAL(IC))
      IF(MASS(1,IC)+MASS(2,IC).LT.1E-5) GO TO 80
 1130 FORMAT(A8,2F8.4,I4,A1,1X,A8,2F8.4,F8.4)
      CALL CHECK(IC,MXP,2)
      if(final) then
      WRITE(KO,1140) IC,NAME(1,IC),MASS(1,IC),MASS(2+1,IC),NEX(IC),
     X   PWF,NAME(2,IC),MASS(2,IC),MASS(2+2,IC),QVAL(IC)
      else
      WRITE(KO,1141) IC,NAME(1,IC),MASS(1,IC),MASS(2+1,IC),NEX(IC),
     X   PWF,NAME(2,IC),MASS(2,IC),MASS(2+2,IC),QVAL(IC)
      endif
 1140 FORMAT(//' *********** PARTITION NUMBER ',I2,1X,90('*'),
     X //'   PROJ=',A8,'   MASS=',F8.4,' Z=',F6.1,',        # STATES=',
     X I5,L1, ', TARG=',A8,' MASS=',F8.4,' Z=',F6.1,
     X  ',     Q-VALUE =',F9.4,' MEV.')
 1141 FORMAT(//' Partition ',I2,' ***',
     X /'   PROJ=',A8,'   MASS=',F8.4,' Z=',F6.1,',        # STATES=',
     X I5,L1, ', TARG=',A8,' MASS=',F8.4,' Z=',F6.1,
     X  ',     Q-VALUE =',F9.4,' MEV.')
C                                         Piecewise flags
      PWFLAG(IC)= PWF .and. FCWFN
C
      RMASS(IC) = MASS(1,IC)*MASS(2,IC) / (MASS(1,IC) + MASS(2,IC))
      HP(IC) = HCM
C     IF(INH.GT.0) HP(IC) = HCM * MASS(INH,1)/MASS(INH,IC)
      IF(INH.GT.0) HP(IC) = HCM *
     X    (MASS(MOD(INH-1,2)+1,1)/MASS(MOD(INH-1,2)+1,IC))**((INH+1)/2)
      GIVEXS(IC) = NEX(IC).GT.0
      NEX(IC) = ABS(NEX(IC))
      NCHPMAX(IC) = 0
      I = NEX(IC)
      IF(I.EQ.0) GO TO 75
      CALL CHECK(I,MXX,1)
      DO 40 IA=0,9
      DO 40 J=1,2
      DO 40 IN=1,2
   40 BHEAD(IN,J,IA) = 0
C
      DO 60 IA=1,NEX(IC)
!     READ(KI,1150) JEX(1,IC,IA),ICOP(1),BAND(1,IC,IA),ENEX(1,IC,IA),
!    X         (JEX(1+I,IC,IA),I=2,4,2),  CPOT(IC,IA),
!    X         JEX(2,IC,IA),ICOP(2),BAND(2,IC,IA),ENEX(2,IC,IA),
!    X         (JEX(2+I,IC,IA),I=2,4,2), (EX(I),I=1,2),IOFAM
      call readst(KI,JEX(1,IC,IA),ICOP(1),BAND(1,IC,IA),ENEX(1,IC,IA),
     X         JEX(1+2,IC,IA),JEX(1+4,IC,IA),  CPOT(IC,IA),
     X         JEX(2,IC,IA),ICOP(2),BAND(2,IC,IA),ENEX(2,IC,IA),
     X         JEX(2+2,IC,IA),JEX(2+4,IC,IA), EX(1),EX(2),
     X         IOFAM(1,IC,IA),IOFAM(2,IC,IA))
 1150 FORMAT(F4.1,2I2,F8.4,2F4.1,I4,2X,F4.1,2I2,F8.4,2F4.1,3A2)
      IT = IT + 1
      ITC(IC,IA) = IT
      IF(CPOT(IC,IA).EQ.0) CPOT(IC,IA) = IC
      IF(ICOP(1).LT.0) ICOP(2) = ICOP(1)
      IF(ICOP(1).EQ.0.AND.ICOP(2).LT.0) ICOP(2) = -IC
      DO 50 IN=1,2
         EXTRA(IN,IC,IA) = EX(IN)
         IAP = ICOP(IN)
         IF(IAP.EQ.0 .AND. BAND(IN,IC,IA).EQ.0) IAP = IA-1
         COPY(IN,IC,IA,1) = IAP
       IF(IAP.EQ.0) THEN
            J = NINT( 1.5 + 0.5*SIGN(1,BAND(IN,IC,IA)) )
           IF(BHEAD(IN,J,ABS(BAND(IN,IC,IA))).NE.0) THEN
            JEX(IN+2,IC,IA)=JEX(IN+2,IC,BHEAD(IN,J,ABS(BAND(IN,IC,IA))))
            JEX(IN+4,IC,IA)=JEX(IN+4,IC,BHEAD(IN,J,ABS(BAND(IN,IC,IA))))
           ELSE
             BHEAD(IN,J,ABS(BAND(IN,IC,IA))) = IA
             JEX(IN+2,IC,IA) = MIN(ABS(JEX(IN+2,IC,IA)),JEX(IN,IC,IA))
             JEX(IN+2,IC,IA) = MAX(JEX(IN+2,IC,IA) ,JEX(IN,IC,IA)
     X                                      - AINT (JEX(IN,IC,IA)) )
             JEX(IN+4,IC,IA) = MAX(JEX(IN+4,IC,IA),
     X          ABS(NINT(MASS(IN,IC))-2*MASS(IN+2,IC))) * 0.5
           ENDIF
       ELSE
C ICOP(IN) NE 0
        IF(ICOP(IN).LT.0) THEN
C                          COPY SAME-NUMBERED PAIR OF EXCITED STATES
C                          FROM PARTITION IC2=/ICOP(1)/
C                               WITH PROJECTILE & TARGET EXCHANGED
         IC2 = ABS(ICOP(IN))
         IN2 = 3 - IN
         IAP = IA
         FLIP = .TRUE.
       ELSE
C                          COPY STATE FROM THIS PARTITION & NUCLEUS
         IC2 = IC
         IN2 = IN
       ENDIF
         IF(COPY(IN2,IC2,IAP,1).NE.0) IAP = COPY(IN2,IC2,IAP,1)
      JEX(IN,IC,IA) = JEX(IN2,IC2,IAP)
      BAND(IN,IC,IA) = BAND(IN2,IC2,IAP)
      ENEX(IN,IC,IA) = ENEX(IN2,IC2,IAP)
      JEX(IN+2,IC,IA) = JEX(IN2+2,IC2,IAP)
      JEX(IN+4,IC,IA) = JEX(IN2+4,IC2,IAP)
C
C     THE ARRAY COPY(IN,IC,IA,1) > 0   OR   COPY(IN,IC,IA,2) > 0
C               (this partition)        (other partition)
C           INDICATES COPIED EXCITED STATE.
         COPY(IN,IC,IA,1) = IAP
         IF(IA.EQ.IAP)         COPY(IN,IC,IA,1) = 0
C     THE ARRAY COPY(IN,IC,IA,2) > 0 INDICATES EXCHANGE COPIES ONLY.
         COPY(IN,IC,IA,2) = IC2
         IF(ICOP(IN).GE.0)     COPY(IN,IC,IA,2) = 0
       ENDIF
   50 CONTINUE
      WRITE(KO,1160) IA,   JEX(1,IC,IA),PSIGN(SIGN(1,BAND(1,IC,IA))+2),
     X    BAND(1,IC,IA), ENEX(1,IC,IA),(JEX(1+I,IC,IA),I=2,2,2),
     X  CPOT(IC,IA),  JEX(2,IC,IA),PSIGN(SIGN(1,BAND(2,IC,IA))+2),
     X    BAND(2,IC,IA), ENEX(2,IC,IA),(JEX(2+I,IC,IA),I=2,2,2)
 1160 FORMAT(/1X,I5,':  ',2('J=',F4.1,A1,' (B#',I2,'), E=',
     X       F8.4,', K=',F4.1,'    ',4X  ,:,'      Potl#',I3,7X) )
      IF(COPY(1,IC,IA,1).NE.0) WRITE(KO,1170) COPY(1,IC,IA,1)
      IF(COPY(2,IC,IA,1).NE.0) WRITE(KO,1180) COPY(2,IC,IA,1)
      IF(COPY(1,IC,IA,2).NE.0) WRITE(KO,1190) COPY(1,IC,IA,2)
      IF(COPY(2,IC,IA,2).EQ.IC)  WRITE(KO,1195)
 1170 FORMAT(' ',33X,'(Just State',I3,')')
 1180 FORMAT(' ',87X,'(Just State',I3,')')
 1190 FORMAT(' ',33X,'(Just Exchange copy from partition',I3,')')
 1195 FORMAT(' ',87X,'(Just Exchange copy of projectile)')
      IF(EXTRA(1,IC,IA)) WRITE(KO,1200) 
 1200 FORMAT('+',87X,14X,' Make cross sections for 180-theta')
         EXTRA(2,IC,IA) = EX(2)
      IF(EXTRA(2,IC,IA)) PRINT 1201
 1201 FORMAT(' ',87X,20X,'IGNORE ITERATED CHANGES')
      IF(IOFAM(1,IC,IA)/=0) 
     x   WRITE(KO,1202) abs(IOFAM(1,IC,IA)),IOFAM(1,IC,IA)>0
 1202 FORMAT('+',87X,14X,' Read in EXTRA angle amplitudes file',i4/
     X       ' ',87X,14X,' With ALL Spin Cases =',L2)
      IF(IOFAM(2,IC,IA)/=0) 
     x   WRITE(KO,1203) abs(IOFAM(2,IC,IA)),IOFAM(2,IC,IA)>0
 1203 FORMAT('+',87X,14X,' Write out EXTRA angle amplitudes file',i4/
     X       ' ',87X,14X,' With ALL Spin Cases =',L2)

      NCHPMAX(IC) = NCHPMAX(IC) + 
     X		nint(JEX(1,IC,IA)+1)*nint(JEX(2,IC,IA)+1)
   60 CONTINUE
!	write(KO,*) '     Expect up to ',NCHPMAX(IC),' partial waves'
      DM = MASS(1,IC)+MASS(2,IC) - MASS(1,PEL) - MASS(2,PEL)
     X         + (QVAL(IC) - QVAL(PEL))/931.481
      DZ = MASS(2+1,IC)+MASS(2+2,IC) - MASS(2+1,PEL) - MASS(2+2,PEL)
      IF(ABS(DM).LT.0.001 .AND. ABS(DZ).LT.0.001.or..not.final) GO TO 70
      WRITE(KO,1210) PEL,DM,DZ
 1210 FORMAT('0****** WARNING : COMPARED WITH PARTITION',I3,', THIS',
     X' PARTITION HAS TOTAL MASS & CHARGE EXCESSES',F10.5,' &',F9.4)
C     IF(ABS(DM).GT.0.10 .OR. ABS(DZ).GT.0.10) STOP
      IF(ABS(DM).GT.0.10 .OR. ABS(DZ).GT.0.10) CALL ABEND(4)
c     IF(                     ABS(DZ).GT.0.10) CALL ABEND(4)
C
   70 IF(INH.GT.0.AND.IC.GT.1) WRITE(KO,1220) INH,HP(IC)
 1220 FORMAT('0INH =',I2,' => radial step sizes =',F8.5,' in this ',
     X 'partition and its potentials')
      IF(PWFLAG(IC)) WRITE(KO,1225)
 1225 FORMAT('0Acceleration of long-range Coulomb couplings in this',
     X       ' partition using CRCWFN')
C
   75 IC = IC + 1
      GO TO 30
   80 NCHAN = IC - 1
      ITCM = IT
      CALL CHECK(ITCM,MXPEX,21)
      RETURN
      END

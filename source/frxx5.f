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

*****INFORM*************************************************************
      SUBROUTINE INFORM(FORML,FORMC,NSP,FORMF,NF,PTYPE,TR,QVAL,BPHASE,
     &   QNF,D0,BEE,NCHAN,NEX,ENEX,JEX,BAND,MASS,AFRAC,ITC,COPY,
     &   N,NLN,MINT,NNN,RNN,RMI,HCM,RIN,MR,NNU,ERANGE,DK,PI,NAME,
     &   NBINS,EMID,DELTAE,KMINX,NKBIN,NORBIN)
	use parameters
	use io
	use drier
	use searchdata
	use searchpar
	use fresco1, only: FATAL,cxwf,sumccbins
      IMPLICIT REAL*8(A-H,O-Z)
C				mp should be defined here the same as in fread
	parameter (mp=200)
      REAL*8 FORMR(MAXM),FORML(MAXNLR,MSP,2),FFR4,KMINX(2,MSP)
      REAL*8 EMID(MSP),DELTAE(MSP),autowid
      COMPLEX*16 FORMF(MAXM,MLOC),TC,FORMC(MAXNLC,MSP,2),FRMC(MAXM),FFC4
      INTEGER QNF(19,MSP),BAND(2,MXP,MXX),NEX(MXP),TR,TNT(4,mp),
     &        PTYPE(8,MLOC),PARITY,POTK,ITC(MXP,MXX),COPY(2,MXP,MXX,2),
     & 	      NKBIN(2,MSP),NORBIN(MSP)
      REAL*8 D0(MSP,2),JEX(6,MXP,MXX),MASS(4,MXP),WL(20),WLD(20),
     &       J,JN,JNMAX,JNMIN,JCORE,JCOM,KCORE,KCOM,MA,MB,
     &       JNA,JNB,JCOREA,JCOREB,KCOREA,KCOREB,K,QVAL(MXP)
      REAL*8 ENEX(2,MXP,MXX),AFRAC(MXPEX,MXPEX,2,MSP),BEE(MSP,5),
     &       TRITON(8),COEF(mp),BPHASE(NKMAX,MSP),FORMIN(2*MAXM),HCM
      real*8,allocatable:: PSI(:,:),CC(:,:,:,:)
      complex*16,allocatable:: PSIC(:,:)
      CHARACTER*3 VAR(18),ADJ
      CHARACTER*8 NAME(2,MXP),TFACTS,TBIN*6,NORMS*4,KFA*2
      CHARACTER CH1,CH2,PSIGN(3)
      data PSIGN / '-','?','+' /
      CHARACTER*120 COMMENT
      LOGICAL FAIL3,FRAC,LAST1,EIGEN,TRES,TDEL,op,keep,cmb,rlb
      LOGICAL TKNRM,PWF,TKMAT,REN(MLOC),FUSED(1000),LIAP,VPOT
     	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN
	namelist/overlap/ kn1,kn2,ic1,ic2,in,kind,ch1,nn,l,lmax,sn,
     &         ia,j,ib,kbpot,krpot,be,isc,ipc,nfl,nam,ampl,keep,
     & 	       dm,nk,er,e,rsmin,rsmax,rsalpha,nlag,phase,autowid
	namelist/dalitz/ triton
	namelist/twont/ tnt,coef
      DATA VAR /'BE','VR','WR','VSO','USO','VTR','UTR','T2R','DEF',
     &          'DEF','DF0','DF1','DF2','DF3','DF4','DF5','DF6','DF7' /
      DATA TRITON /  .224, 0.50, 1.00, 1.38, 5.5, 0.0, 1.38, 4.2405 /
      FRAC(X) = ABS(X-NINT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
C     the 'QNF' array gives the channel & quantum numbers of
C                the form factors (both 1 & 2 particle bound states)
C         for each KN=KN1,KN2   QNF(i,KN) gives
C
C     QNF(1 : KN1 = no. of a parent form factor for which
C                   fractional parentages etc are defined
C         2 : ICR =  core partition
C         3 : IA  =  core excitation pair (or zero, if not specified)
C         4 : ICP =  compound nucleus partition
C         5 : IB  =  compound nucleus excitation pair (or zero)
C         6 : IN = 1 for projectile state, = 2 for target state
C
C         7 : KIND of bound state (0-3 = 1N  and  6-9 = 2N)
C                  = 0 :  ln,sn; jn                (any IA,IB)
C                  = 1 :  ln, (sn,Jcore)jn; Jcom   (fixed IA,IB)
C                  = 2 :  ln,sn; jn part of sn/K/parity deformed state
C                  = 3 : (ln,sn)jn, Jcore; Jcom    (fixed IA,IB)
C                  = 4 :  (Dalitz-Thacker)
C                  = 5 :
C                  = 6 :  ln, (l,s)sn; jn            & (.5,5)T
C                  = 7 : (ln,l)jn, (s,Jcore)sn; Jcom & (.5.5)T,Tcor;Tcom
C                  = 8 :
C                  = 9 : (ln,(l,s)sn)jn,Jcore; Jcom  & (.5.5)T,Tcor;Tcom
C
C         8 : NN = no. of nodes (incl origin & not infinity, so NN.ge.1)
C         9 : ln = L values relative to core nucleus
C        10 : 2*sn= 2 * total spin of bound cluster  (excl. ln & Jcore)
C                               (but for KIND 7, excluding l too)
C        11 : 2*jn= 2 * (ln + sn as vectors)
C                               (but for KIND 7, jn = ln + l = 'lambda',
C                                and for KIND 1, jn = sn + Jcore = 'S')
C        12 : n  =  0 for single-particle states
C                =  index of nucleon-nucleon separation for 2N states
C        13 : l  =  l value between two nucleons        for 2N states
C        14 : 2*S=  2 * combined NN spin  = 2*(0. or 1.)for 2N states
C                  so (vec) sn = l + S                  for 2N states
C        15 : 2*T=  2 * combined NN isospin=2*(0 or 1)  for 2N states
C                  so l + S + T is odd                  for 2N states
C        16 : Type = character identifier for otherwise-structerless
C                    bound clusters.  A-M for + parity, O-Z for -
C        17 : Number of coupled channels (1, for single channel)
C        18 : IL = Incoming wave for multichannel continuum bins
C        19 : KN1 = no. of parent formfactor for IL=1
C
C        where  Jcore = spin of core nucleus = JEX(IN,ICR,IA)
C               Jcom  = spin of compound nuc = JEX(IN,ICP,IB)
C                  so (vec) Jcom = Jcore + jn
C               Kcore = K-value of core nuc  = JEX(IN+2,ICR,IA)
C               Kcom  = K-value of compound  = JEX(IN+2,ICP,IB)
C                  and      K = Kcom - Kcore
C               Tcore = isospin of core nuc  = JEX(IN+4,ICR,IA)
C               Tcom  = isospin of compound  = JEX(IN+4,ICP,IB)
C                  so (vec) Tcom = Tcore + T
C          and  Parity = sign(1, BAND(IN,ICP,IB) * BAND(IN,ICR,IA))
C
C     the 'BEE' array gives the energies and norms of
C                the form factors (both 1 & 2 particle bound states)
C         for each KN=KN1,KN2   BEE(KN,i) gives
C
C      BEE(KN,1 : BE  = binding energy of state (negative if unbound)
C           ,2 : NORM= root-mean-square norm of wavefunction
C           ,3 : BEE(KN,3) = rms radius of state
C           ,4 : BEE(KN,4) = ANC of state
C           ,5 : ETAP=2*k*ETA      "     "   "    "
C  but BEE(KN,1 : k**2= BE*CONV  (during call to EIGCC)
C
      EPSCON= 1E-4
      PWF = .false. 	! may need in the future
      SMALL = 1D0/SQRT(FPMAX)
      MAXC = 40
      Z = 0.0
       AFRAC(1:MXPEX,1:MXPEX,1:2,1:MSP) = 0.0
C
C
C    single-particle form factors and their parameters
C    -------------------------------------------------
      DO 20 KN=1,MSP
      DO 20 I=1,17
20    QNF(I,KN) = 0
      NSP = 0
      NBINS = 0
      LAST1 = .FALSE.
      FUSED(:) = .false.
	FORMR(:) = 0.
C
	inquire(iolength=iol) FORMR
	inquire(8,opened=op) 
      if(.not.op) then
      IF(MACH==8.or.MACH==3.or.MACH==4) THEN
           OPEN(8,RECL=4*iol,FILE=TMPN(1:lnbl(TMPN))//'08',
     X     ACCESS='DIRECT',FORM='UNFORMATTED')
        ELSE
           OPEN(8,ACCESS='DIRECT',RECL=4*iol,STATUS='SCRATCH',
     X     FORM='UNFORMATTED')   ! 4 rather than 2, to allow for complex bins
        ENDIF
	endif
!	write(6,*) ' FILE 8 opened, recl =',2*iol
21    FORMAT(//' ',132('*'),//,
     &' The following SINGLE-PARTICLE FORM FACTORS are constructed :'/)
22    FORMAT(/' No.   P1 P2 IN KIND T N  L  S1 IA J/S  IB ',
     &  'BIND XFER  BE    SC PC FIL AFRAC Adjust  to',
     &  '  Z  Mass   K     Norm    rms      D0     D   ANC/Gsp')
      DO 900 KNP=1,1000
!      READ(KI,852) KN1,KN2,IC1,IC2,INI,KIND,CH1,NN,L,LMAX,SN,IAK,J,IB,
!     &         KBPOT,KRPOT,E,ISC,IPC,NFL,NAM,AMPL
	KN1=0;IC1=0;IC2=0;IN=0;BE=0.0;IA=0;NK=0
	read(KI,nml=overlap)
852   FORMAT(2I3,4I2,1X,A1,3I2,F4.1,I2,F4.1,I2,2I3,F8.4,4I3,3F8.4)
	IAK=IA
	E = BE
      IF(TR.GE.4)
     &WRITE(KO,852) KN1,KN2,IC1,IC2,IN,KIND,CH1,NN,L,LMAX,SN,IAK,J,IB,
     &         KBPOT,KRPOT,E,ISC,IPC,NFL,NAM,AMPL,rsmin,rsmax
      KN2 = MAX(KN2,KN1)
      KN2I= KN2
      INI = ABS(IN)
      IN = ABS(IN)
      IF(KN1*IC1*IC2*IN.LT.1)GO TO 1000
      CALL CHECK(KN2,MSP,3)
C     if(knp.eq.1)  write(KO,8858) 8,n,'complex*8'
C8858       FORMAT(/' File ',I3,' needs RL =',I3,' ',A9,' numbers')
      IF(KNP.EQ.1) WRITE(KO,21)
      DZ = MASS(2+IN,IC1) - MASS(2+IN,IC2)
      DMM = MASS(IN,IC1) - MASS(IN,IC2)
      ICR = IC2
      ICP = IC1
      IF(DMM.LT.0) ICR = IC1
      IF(DMM.LT.0) ICP = IC2
      IF(IC1.GT.NCHAN.OR.IC2.GT.NCHAN.OR.(IN.NE.1.AND.IN.NE.2))THEN
            WRITE(KO,23) KN1,KN2,IC1,IC2,INI
23    FORMAT(/' UNUSABLE COMBINATION OF CHANNEL AND/OR QUANTUM NUMBERS :
     & '   //,' ',2I4,':',3I3,1X,'?'/)
            GO TO 890
            ENDIF
C     IF(IN.EQ.INH) HDN = HP(ICR)
C     IF(IN.NE.INH) HDN = HP(ICP)
C             STEP SIZES HERE ARE INDEPENDENT OF PARTITION!
      HDN = HCM
C     z1z2 = mass(2+in,icr) * sign(dz,dm)
      ZC = MASS(2+IN,ICR)
	isearch=0
      NKHERE = 0
      ERHERE = ERANGE
	if(nk/=0) NKHERE = nk
	if(abs(er)>1e-9) ERHERE = er
       AMPL = AMPL * SQRT(REAL(max(NAM,0)))
      RM = DM * MASS(IN,ICR) / (DM + MASS(IN,ICR))
      CONV = FMSCAL * RM
      IF(TR.GE.4) WRITE(KO,*) NF,ICR,ICP,DZ,ZC*ABS(DZ),DM,RM,CONV
      IAI = MIN(MAX(IAK,1),NEX(ICR))
!      IB  = MIN(IB,        NEX(ICP))
      IF(IN==2)THEN
       IB  = MIN(IB,        NEX(ICP)+NEX(ICR))
      ELSE
       IB  = MIN(IB,        NEX(ICP))
      ENDIF
         AB = 1.0
         IF(L.LT.0) AB = -1.0
         L = MAX(L,0)
      IAMIN = 1
      IAMAX = NEX(ICR)
      LMAX = MAX(LMAX,L)
      IF(KIND.LE.2 .OR. KIND.EQ.6) THEN
         IAMIN= IAI
         IAMAX= IAI
         ENDIF
	LMX1 = 0
C
      IF(KIND.LE.5) THEN
      IF(KNP.EQ.1.OR..NOT.LAST1) WRITE(KO,22)
C                          one-nucleon bound state(s)
C                          --------------------------
      LMIN = 0
      JNMIN = 0.0
      JNMAX = J
      JCOM = JEX(IN,ICP,MAX(IB,1))
C
      IF(KIND.EQ.0) THEN
C                          (ln,sn) jn  coupling order
         LMIN = L
         LMAX = L
         JNMIN= J
         JNMAX= J
C
      ELSE IF(KIND.EQ.1) THEN
C                             (ln, (sn,Jcore)s; Jcom) coupling order
         IF(IB.LE.1) IB = 1
         TCOM = JEX(IN,ICP,IB)
C+++++++++++++++++++++++++++++++ PATCH MARCH 1993
         JNMIN = J
C
      ELSE IF(KIND.EQ.2) THEN
C                             (ln,sn) jn part of a sn/k/parity state
C                             in a deformed potential.
         IB = MAX(IB,1)
         KCOM = JEX(IN+2,ICP,IB)
         KCORE= JEX(IN+2,ICR,IAI)
         K = KCOM - KCORE
         JNMIN = ABS(K)
         JNMAX = LMAX + SN
C
      ELSE IF(KIND.EQ.3) THEN
C                             (ln,sn)jn, jcore; jcom coupling order
         IB = MAX(IB,1)
         JNMAX = LMAX + SN
      ELSE IF(KIND.EQ.4) THEN
C                              construct form for dalitz-thacker triton
!         IF(NN.GE.2) READ(KI,*) TRITON
         IF(NN.GE.2) READ(KI,nml=dalitz)
         IF(IPC.GE.1) WRITE(KO,44) TRITON
44       FORMAT(' Dalitz-Thacker parameters =',8F9.4)
         BEE(KN1,1) = TRITON(8)  * CONV
         NU = TRITON(2)/HDN + 1.1
	 allocate(PSI(MINT,1))
         DO 45 I=1,NU
45       PSI(I,1) = 0.0
         DO 46 I=NU+1,MINT
            X = (I-1)*HDN
            XR= X - TRITON(2)
46     PSI(I,1) = (EXP(-TRITON(1)*XR) + (TRITON(4)-1)*EXP(-TRITON(3)*XR)
     &             - TRITON(7)*EXP(-TRITON(5)*XR) )  * SQRT(X)
         M = 1
         KN2 = KN1
         VARY = 0.0
         ADJ = 'DZT'
         EIGEN = .TRUE.
          QNF(1,KN1) = KN1
          QNF(2,KN1) = ICR
          QNF(4,KN1) = ICP
          QNF(6,KN1) = IN
          QNF(7,KN1) = KIND
          QNF(8,KN1) = 1
          QNF(9,KN1) = 0
          QNF(17,KN1) = 1
         IF(MOD(IPC,2).EQ.1) WRITE(KO) 826,HDN*MR,(PSI(I,1),I=1,MINT,MR)
826       FORMAT(' DT wavefunction at intervals of',F6.3,' from 0 is',
     &           /(1X,12F10.5))
         GO TO 50
      ENDIF
C
      IF(E.EQ.0.0) E = 10.0
      EIGEN = E.GT.0.0
C                      if eigen then bound state, else continuum bin
         IF(ISC.EQ.0.AND.EIGEN) THEN
            THETA = CONV
            VARY = E
         ELSE
            THETA = 0.0
            VARY = 1.0
         ENDIF
      IF(.NOT.EIGEN) THEN
C                      continuum bin 
	    if(ISC==0)  ISC=2
            I10 = mod(abs(ISC),10)
            TRES  = I10.EQ.3 .OR. I10.EQ.4 .OR. I10.EQ.7 .OR. I10.EQ.8
            TDEL  = I10.EQ.1 .OR. I10.EQ.2
            TKMAT = I10.GE.5 .AND. I10.LE.8
            TKNRM = ISC.ge.10
            TFACTS = ' without'
            IF(TRES) TFACTS = '    WITH'
            IF(TDEL) TFACTS = 'phase of'
            IF(TKMAT)TFACTS = ' (1-iK) '
            IF(TKMAT.and.TRES)TFACTS = 's(1-iK) '
            TBIN =           '  Real'
!           if(M>1) TBIN =   '~~Real'     ! M only known below
            if(ISC<0) TBIN = 'Complx'
         ENDIF
      KN = KN1-1
      PARITY = 0
       IK = 0
      DO 30 IA=IAMIN,IAMAX
c      DO 30 IA=IAMAX,IAMIN,-1
         JCORE = JEX(IN,ICR,IA)
         KCORE = JEX(IN+2,ICR,IA)
C        tcore = jex(in+4,icr,ia)
         IF( COPY(IN,ICR,IA,1).NE.0 ) GO TO 30
         IF(MOD(KIND,2).EQ.1)
     &      PARITY = SIGN(1,BAND(IN,ICP,IB)*BAND(IN,ICR,IA))
         IF(IPC.GE.5) WRITE(KO,*) 'TRY IA =',IA,' SO PARITY =',PARITY
      DO 29 LN=LMIN,LMAX
c      DO 29 LN=LMAX,LMIN,-1
         IF(IPC.GE.5) WRITE(KO,*) 'TRY LN =',LN
          IF(MOD(KIND,2).EQ.1 .AND. (-1)**LN.NE.PARITY  .OR.
     &       MOD(KIND,2).EQ.0 .AND. (-1)**(L+LN).NE.1)  GO TO 29
!      DO 28 JN=JNMIN,JNMAX,0.5
      NJN=NINT(JNMAX-JNMIN)*2
      DO 28 IJN=0,NJN
c      DO 28 IJN=NJN,0,-1
      JN=JNMIN+IJN*0.5
C
      IF(KIND.NE.1 .AND. FAIL3(LN+Z,JN,SN)) GO TO 28
         IF(IPC.GE.5) WRITE(KO,*) 'TRY JN =',DBLE(JN)
      IF(KIND.EQ.1 .AND.(FAIL3(SN,JCORE,JN) .OR.
     &                   FAIL3(LN+Z,JN,JCOM))) GO TO 28
      IF(KIND.EQ.3 .AND. FAIL3(JN,JCORE,JCOM)) GO TO 28
      KN = KN + 1
         IF(KN.GT.KN2) GO TO 28
      QNF(1,KN)  = KN1
      QNF(2,KN)  = ICR
      QNF(3,KN)  = IA
      QNF(4,KN)  = ICP
      QNF(5,KN)  = IB
      QNF(6,KN)  = IN
      QNF(7,KN)  = KIND
      QNF(9,KN)  = LN
      QNF(10,KN) = NINT(2.*SN)
      QNF(11,KN) = NINT(2.*JN)
         DO 25 I=12,15
25       QNF(I,KN) = 0
      QNF(16,KN) = ICHAR(CH1)
      IF(IA.EQ.IAI .AND. L.EQ.LN .AND. ABS(JN-J).LT..1) IK=KN
	LMX1 = max(LMX1,LN+1)
      BEE(KN,1) = E * CONV
      IF(KIND.EQ.1 .OR. KIND.EQ.3)
!     &   BEE(KN,1) = (E + ENEX(IN,ICR,IA) - ENEX(IN,ICR,IAI)) * CONV
     &   BEE(KN,1) = (E + ENEX(IN,ICR,IA)) * CONV
      BEE(KN,1) = BEE(KN,1) - THETA * VARY
      BEE(KN,5) = ZC*ABS(DZ) * ETACNS * 2 * RM  * SQRT(FMSCAL) * AB
       if(ia>0.and.ib>0.and.icr>0.and.icp>0) 
     #       AFRAC(ITC(ICP,IB),ITC(ICR,IA),IN,KN) = AMPL
      IF(IPC.GE.5) WRITE(KO,27) KN,(QNF(I,KN),I=1,16),BEE(KN,1)/CONV
27    FORMAT(/' BOUND PARTIAL WAVE @',I3,' IS',16I4,', ASYMP BE =',F8.4)
28    CONTINUE
29    CONTINUE
30    CONTINUE
      IF(KN.GT.KN2) THEN
           WRITE(KO,31) KN-KN2
31         FORMAT(/' NO ROOM FOR',I4,' MORE CHANNELS: INCREASE KN2')
           KN = KN2
        ENDIF
      KN2 = KN
      M = KN2 - KN1 + 1
      QNF(17,KN1:KN2) = M
         IF(M.EQ.0) THEN
          WRITE(KO,*) '  NO SUITABLE CHANNELS FOUND ????  INPUT DATA WER
     XE:'
        WRITE(KO,852) KN1,KN2,IC1,IC2,INI,KIND,CH1,NN,L,LMAX,SN,IAK,J,IB
     &        ,KBPOT,KRPOT,E,ISC,IPC,NFL,NAM,AMPL
               GO TO 890
             ENDIF
      IL= IK  - KN1 + 1
        IF(IL.le.0.or.IL.gt.M) IL=1
      QNF(18,KN1:KN2) = IL
      QNF(19,KN1:KN2) = QNF(1,KN1)
      IF(IL>1)THEN
       DO KNA=KN1-1,1,-1   ! search back for parent with IL=1
        LIAP=.TRUE.
        DO IQNF=1,11     ! parent will have identical quantum numbers
         IF(IQNF==3 .OR. IQNF>=9)THEN
          IF(QNF(IQNF,KNA)/=QNF(IQNF,KN1) .OR.
     &       abs(JEX(QNF(6,KNA),QNF(4,KNA),QNF(5,KNA))
     &           -JEX(QNF(6,KN1),QNF(4,KN1),QNF(5,KN1)))>0.001 .OR.
     &       abs(BEE(KNA,1)-BEE(KN1,1)/CONV)>0.001 .or.
     &       QNF(18,KNA)/=1 ) LIAP=.FALSE.
         ENDIF
        ENDDO
        IF(LIAP)THEN
         QNF(19,KN1:KN2) = QNF(1,KNA)
         EXIT
        ENDIF
       ENDDO
      ENDIF
      IF(IPC.GE.5)
     &WRITE(KO,*)' IEXT=',QNF(5,KN1),' ALPHA=',QNF(5,QNF(19,KN1))
C
      IF(NFL.GT.0) THEN
C                       Just read in wave functions from file NFL
        ADJ = 'FIL'
!        VARY = 0.0
!	rlb = .true.; cmb = .false.; EIGEN=.true. ! only real for now
	cmb = .not.EIGEN
	rlb = EIGEN
        GO TO 50
       ENDIF
C                                          couplings
C                                          ---------
	allocate(CC(M,M,NF,3))
      DO 40 INA=1,M
         KNA = KN1 + INA - 1
         LNA = QNF(9,KNA)
         JNA = QNF(11,KNA)*0.5
         JCOREA = JEX(IN,ICR,QNF(3,KNA))
         KCOREA = JEX(IN+2,ICR,QNF(3,KNA))
      DO 40 INB=1,M
         KNB = KN1 + INB - 1
         LNB = QNF(9,KNB)
         JNB = QNF(11,KNB)*0.5
         JCOREB = JEX(IN,ICR,QNF(3,KNB))
         KCOREB = JEX(IN+2,ICR,QNF(3,KNB))
C
      REN(:) = .false.
      DO 35 JF=1,NF
         KP = PTYPE(1,JF)
         POTK = PTYPE(2,JF)
         DO 32 I=1,3
32       CC(INA,INB,JF,I) = 0.0
          IF(PTYPE(4,JF).LT.0) GO TO 35
          IF(PTYPE(3,JF).LT.0.OR.(KP.NE.KBPOT.AND.KP.NE.KRPOT)) GO TO 35
          IF(PTYPE(2,JF).GE.10 .AND.
     X       PTYPE(3,JF).EQ.0  .AND. INA.NE.INB) GO TO 35
      IF(KIND.EQ.0)
     &   T = TENS0 (POTK,LNA,SN,JNA,JCOREA,LNB,SN,JNB,JCOREB,
     &              PTYPE(3,JF),ABS(DZ),ZC)
      IF(KIND.EQ.1)
     &   T = TENSLS(POTK,LNA,SN,JNA,JCOREA,JCOM,LNB,SN,JNB,JCOREB,
     &           PTYPE(3,JF),SN,KCOREA,SN,KCOREB,ABS(DZ),ZC,PTYPE(5,JF))
      IF(KIND.EQ.2)
     &   T = TENDEF(POTK,LNA,SN,JNA,LNB,JNB,PTYPE(3,JF),K,ABS(DZ),ZC)
      IF(KIND.EQ.3)
     &   T = TENSOR(POTK,LNA,SN,JNA,JCOREA,JCOM,LNB,SN,JNB,JCOREB,
     &           PTYPE(3,JF),SN,KCOREA,SN,KCOREB,ABS(DZ),ZC,PTYPE(5,JF))
      IF(ABS(T).LT.SMALL) GO TO 35
      T = - T * CONV
     X      * (-1)**NINT((JCOREA-JCOREB - ABS(JCOREA-JCOREB))/2.)
C           The above phase factors with JCORE etc
C           are there only because of definition of M(Ek) matrix element

         IF(KP.EQ.KBPOT) THEN
            IF(ISC.EQ.0 .OR. POTK.EQ.0 .OR. .NOT.EIGEN) THEN   ! bin potls not adjusted yet
                I = 1
            ELSE IF(abs(ISC).EQ.POTK) THEN
                I = 2
            ELSE IF((abs(ISC).EQ.8.OR.abs(ISC).EQ.9) .AND.
     & 		    PTYPE(4,JF).GT.0) THEN
                I = 2
            ELSE IF(abs(ISC).GE.10.AND.PTYPE(4,JF).GT.0
     &                       .AND.PTYPE(3,JF).EQ.abs(ISC)-10) THEN
                I = 2
            ELSE
                I = 1
            ENDIF
C                               i=1 is the fixed  & i=2 the varied part
         CC(INA,INB,JF,I) = T
	 REN(JF) = REN(JF) .or. I==2
         ENDIF
         IF(KP.EQ.KRPOT)     CC(INA,INB,JF,3) = T
      IF(IPC.GE.5) WRITE(KO,34) KNA,LNA,JNA,JCOREA,KNB,LNB,JNB,JCOREB,
     &      JF,(PTYPE(I,JF),I=1,5),(-CC(INA,INB,JF,I)/CONV,I=1,3)
34    FORMAT(' Between',I4,' (',I3,2F4.1,') &',I4,' (',I3,2F4.1,') ',
     &       ' for potl at',I3,' (',5I3,') --- get',3F10.4)
35       CONTINUE
40     CONTINUE
C-------------------------------------------------------------------
	cmb = .not.EIGEN
	rlb = EIGEN
	if(nlag>0) then
	rlb = .true.; cmb = .false.  ! complex (cc) bins not allowed yet
	  write(KO,*) 
!	  write(KO,*) ' ************ '
	  write(KO,*) ' Find overlap with ',nlag,' Lagrange mesh basis'
	  VPOT = phase/=0.
	  if(VPOT) then 
           write(KO,401) phase
401	  format('   Fix potential giving phase shift of',f8.3,' deg')
           if(abs(autowid)>1e-5) write(KO,4011) autowid
4011      format('    and set bin width =',f8.3,'*Gamma')
	  endif
	  open (900,file='pluto.in',form='formatted')
	  write(KO,*) '  Read potential parameters from file ',900+KBPOT
	  rewind(900+KBPOT)
	  do i=1,100
	  read(900+KBPOT,'(A)',end=402) COMMENT
	  write(900,'(A)') COMMENT
	  enddo
402	  continue
	IA = IAMIN ! guess
	 if(NFL==0)  NFL = -33
	 write(900,403)  'particle', DM,abs(DZ),SN,
     x           NAME(IN,ICR),MASS(IN,ICR),ZC,
     X           1,JEX(IN,ICR,QNF(3,KN1)),ENEX(IN,ICR,IA),abs(NFL)
403	 format('  name1=''',a8,'''  mass1=',f8.4,' z1=',f8.3,
     x           ' spin1=',f5.1,/
     x           '  name2=''',a8,'''  mass2=',f8.4,' z2=',f8.3,/
     x           '  ncore=',i1,' icore=',f5.1,' encore=',f8.3,/
     x           '  nbin =',i4)
	 IF(MOD(KIND,2).NE.1)  then
	   PARITY= (-1)** QNF(9,KN1) !  else keep above parity
	   JCOM = QNF(11,KN1)*0.5    ! just  j=l+s
	   endif
 	 write(900,404) nlag,0,isc
404	 format('  nlag=',i3,' njt=',i2,'  plot(:)=0, bin=',i3)

	if(BEE(KN1,1)>0) then  			! bound states
	 write(900,405) emin,emaxlist,JCOM,PARITY
405	 format('  emin=',f8.2, ' emaxlist =',f8.3, 
     x          ' jtot=',f6.1,' parity=',i2)

	else					! bin states
	 de= max(abs(ER)/NK,1d-10)              ! precision of eminscat printing!
	 ! E = -BEE(KN1,1)
	if(VPOT) then 				! vary potential to phase
	 isearch=2
	 write(900,406) isearch,1,NN,-E,phase,autowid,
     x                  JCOM,JCOM,PARITY,PARITY
406	 format('  search=',i2,' sjt=',i1,' enodes=',i2,
     x     ' eigen=',f10.6,' bloch=T phase=',f6.2,' autowid=',f9.3,/,
     x     '  jtot=',2f6.1,' parity=',2i3)
	else					! simply get phases from given potential
	 isearch=0
	 write(900,407) JCOM,PARITY,1
407	 format('  jtot=',f6.1,' parity=',i3,' plot(4)=',i1)
	endif

	write(900,408) -E-abs(er)*0.5,-E+abs(er)*0.5,de
408	format('   eminscat=',f15.9,' emaxscat=',f15.9,' de=',1p,e10.2)
	endif
	 write(900,*) '/'
	 close(900)
	 call system('pluto < pluto.in > pluto.out')
	 Adj='Plu'
!	 rlb = .false.; cmb = .true.
	 allocate(PSI(MINT,M))
         FK = SQRT(FMSCAL * (-E) * RM)
	 go to 50   ! read in new wf from file fort.33


      else IF(EIGEN) THEN
	 allocate(PSI(MINT,M))
         ADJ = VAR(abs(ISC)+1)
        IF(IL.le.0.or.IL.gt.M) IL=1
      CALL EIGCC(PSI,FORMF,CC,NF,QNF(1,KN1),BEE(KN1,5),AB,MR,
     &           BEE(KN1,1),THETA,VARY,MAXM,IL,NN,MINT,HDN,M,IPC,
     &           2*M+1,MAXC,EPSCON,IFAIL,LMX1,MINT)
       if(IFAIL>0) then
        if(number_calls>5) then
	 penalty = penalty + fine
	 write(6,41) number_calls,KN1,IFAIL,penalty
41	 format('  At call ',i5,', sp state ',i3,': IFAIL =',i3,:,
     & 		' so  penalty =',1p,e10.1)
	else
	 write(6,41) number_calls,KN1,IFAIL
	 if(FATAL)  stop 'EIGCC FAILURE'
	endif
       endif
       ELSE
C               CONTINUUM BINS
	 allocate(PSIC(MINT,M))
         FK = SQRT(FMSCAL * (-E) * RM)
        IF(ERHERE.GT.0) THEN
C                            E RATIO = ERANGE
          IF(ERHERE.LT.1.) ERHERE = 1./ERHERE
          T = ERHERE ** 0.25
          MB = FK * T
          MA = FK / T
        ELSE
C                   E DIFFERENCE = ERANGE
         EMIN = -E - ABS(ERHERE)*0.5
         EMIN = MAX(EMIN,1D-5)
         EMAX = EMIN + ABS(ERHERE)
         MB = SQRT(FMSCAL * EMAX * RM)
         MA = SQRT(FMSCAL * EMIN * RM)
        ENDIF
          NK = MAX(10, NKHERE)
	  if(abs(DK)>1e-9) NK = MAX( INT((MB-MA)/DK), NK)
	  NK = min(NK,NKMAX)
	  NBINS = NBINS+1
	  NKBIN(1,NBINS) = KN1
	  NKBIN(2,NBINS) = NK
	  NORBIN(NBINS) = ISC
	  EMID(NBINS) = -E
          DELTAE(NBINS) = EMAX - EMIN
	  KMINX(1,NBINS) = MA
	  KMINX(2,NBINS) = MB
         NORMS = '  no'
         IF(MOD(ISC,2).EQ.1) NORMS = 'unit'
         KFA = '  '
         IF(TKNRM) KFA = '*k'
         if(M>1) TBIN =   '~~Real'
         T = 1.0 / (FMSCAL * RM)
!        IF(IL.le.0.or.IL.gt.M) IL=1         ! moved earlier
      WRITE(KO,49) TBIN,MA**2 * T, MB**2 * T, NK, TFACTS, KFA,NORMS
49    FORMAT(/' ',A6,' Bin from',F11.6,' MeV to',F11.6,' MeV (cm) in',
     X  I4, ' steps, ',A8,' T-matrix (or D0) factors',
     X A2,' and ',A4,' normalising:')
      if(M.gt.1) WRITE(KO,491) IL
491   FORMAT(' Incoming waves in channel ',I3)
!      IF(IPC.ge.3) then
        write(43,492) KN1
492	format('# KN1=',2i4)
        write(44,492) KN1,IPC
!      endif
      CALL BINCC(PSIC,FORMF,CC,NF,M,BEE(KN1,1),IL,CONV,BPHASE(1,KN1),
     &  ISC,MA,MB,NK,QNF(1,KN1),BEE(KN1,5),MINT-1,HDN,IPC,TRES,TDEL,
     &  TKMAT,TKNRM,LMX1,MAXM,ANC)
         ADJ = 'BIN'
         QNF(8,KN) = -1
      ENDIF
C-------------------------------------------------------------------
50       WNORM = 0.0
         WRMS = 0.0
         WD0 =0.0
         WD = 0.0
         IF(NFL.NE.0) then
	    call openif(ABS(NFL))
	    FUSED(abs(NFL)) = .true.
	    endif
        I = 0
          if(IB>0) I = SIGN(1,BAND(IN,ICP,IB))
          IF(M.GE.2) WRITE(KO,51) M,KN1,JCOM,PSIGN(I+2)
51     FORMAT(/' The following',I3,' components are in a group labelled'
     &        ,I3,' for ',f5.1,a1)
      KM = KN1 - 1
      DO 860 KN=KN1,KN2
         IC = KN - KN1 + 1
         BEE(KN,1) = BEE(KN,1)/CONV
         IF(ISC.EQ.0.AND.EIGEN) BEE(KN,1) = BEE(KN,1) + VARY
         FK = SQRT(CONV * ABS(BEE(KN,1)))
              L = QNF(9,KN)
         BEE(KN,2) = 0.0
         BEE(KN,3) = 0.0
         BEE(KN,4) = 0.0
      IF(NFL.GT.0.or.nlag>0) THEN
      IF(EIGEN) then
          if(.not.allocated(PSI)) allocate(PSI(MINT,M))
         else
          if(.not.allocated(PSIC)) allocate(PSIC(MINT,M))
         endif
        
!OLD         READ(NFL,*) (FORMR(I),I=1,MINT)
         READ(abs(NFL),'(a)')  COMMENT
	if(nlag==0.or.isearch==0) then
         READ(abs(NFL),*) NPOINTS,RSTEP,RFIRST
        else
         READ(abs(NFL),*) NPOINTS,RSTEP,RFIRST,VARY  ! get also potential scaling for eigensolution
        endif
         WRITE(KO,71) COMMENT,NPOINTS,RSTEP,RFIRST
71      FORMAT('  Input:',a120,/'  Reading ',I4,'*2 wf+vertex points',
     X   ' at ',F8.4,' intervals, starting from r =',F8.4)
	if(NPOINTS>2*MAXM) then
	  write(KO,*) ' Only room to read in ',2*MAXM,' potential',
     x     ' points, not ',NPOINTS
	  stop
	  endif
          READ(abs(NFL),*) (FORMIN(I),I=1,NPOINTS)
!	 write(300,'(f8.3,g12.4)') (RFIRST+(I-1)*RSTEP,FORMIN(I),
!     x              I=1,NPOINTS)
!	 write(300,*) '&'
          RSTEPI = 1./RSTEP
          DO I=1,MINT
            R = (I-1)*HDN
            FORMR(I) = FF4(R-RFIRST,RSTEPI,FORMIN,NPOINTS)
	    PSI(I,IC) = FORMR(I)  ! only real for now
	  enddo
        ENDIF
        if(final) then
	  call openif(46)
        IF(mod(IPC,2).eq.1.and.EIGEN)  then
              ETA  = 0.5*BEE(KN,5)/FK
          WRITE(46,58) KN,QNF(9,KN),QNF(11,KN)*0.5,QNF(3,KN),FK,ETA
          written(46) = .true.
	 endif
  58   FORMAT('# For Channel #',I3,': l =',I2,', j =',F5.1,', c#',I3,
     x     ';  k,eta =',2f10.5)
        endif

        IF(IC==1.and.IPC.GE.1) WRITE(KO,22)

	  if(PWF) then
	  inquire(55,opened=op)
	  if(.not.op) then
	      call rewop(55); call rewop(59)
            WRITE(55,*) real(HDN),MINT,' : H,N'
	  endif
	  endif
	if(RSMIN>HDN .or. RSMAX<(MINT-2)*HDN) then ! trim radial wfs
        DO I=1,MINT
         R = (I-1)*HDN
	 if(R<RSMIN .or. R>RSMAX) then
	   if(rlb) PSI(I,IC) = 0.0
	   if(.not.rlb) PSIC(I,IC) = 0.0
	 endif
	 ENDDO
	 if(RSMIN>HDN) WRITE(KO,591) 'below',RSMIN
	 if(RSMAX<(MINT-2)*HDN) WRITE(KO,591) 'above',RSMAX
591	 format(/' *** The next WFN is set to zero ',a5,f8.4,' fm ***')
	endif
	
        if(rsalpha>0.) then ! regularize radial wfs
        DO I=1,MINT
         R = (I-1)*HDN
           if(rlb) PSI(I,IC) = PSI(I,IC) * exp(-rsalpha*R)
           if(.not.rlb) PSIC(I,IC) = PSIC(I,IC) * exp(-rsalpha*R)
         ENDDO
         WRITE(KO,592) rsalpha
592      format(/' *** The next WFN is regularized by exp(-',f7.3,'*R)')
        endif


        IF(mod(IPC,2)>0.and.final)  then
		call openif(58)
	 if(PWF)
     x    WRITE(55,*) KN1,IC,M,QNF(9,KN),QNF(11,KN)*.5,-real(BEE(KN,1)),
     X		' : KN,in-set,total-set,l,j,E'
          WRITE(58,155) KN1,-real(BEE(KN,1)),QNF(9,KN),QNF(11,KN)*0.5,
     X		 QNF(11,KN)+1,MINT
	if(PWF)
     x    WRITE(59,155) KN1,-real(BEE(KN,1)),QNF(9,KN),QNF(11,KN)*0.5,
     X		 QNF(11,KN)+1,MINT
155          format('# KN1,E,l,j,2j+1,N =',i4,f10.5,i4,f5.1,2i4)
          if(PWF) written(55) = .true.
          written(58) = .true.
          if(PWF) written(59) = .true.
	 endif
              ETA  = 0.5*BEE(KN,5)/FK
        DO 60 I=1,MINT-1
         R = MAX(I-1.,0.01)*HDN
	   if(rlb) then
           IF(NFL.LE.0.and.nlag==0) FORMR(I) = PSI(I,IC) / R
	    if(mod(IPC,2)>0.and.final) write(58,*) real(R),FORMR(I)*R
           WWW = FORMR(I)
	   else
           IF(NFL.LE.0.and.nlag==0) FRMC(I) = PSIC(I,IC) / R
	    if(mod(IPC,2)>0.and.final) 
     x          write(58,*) real(R),real(FRMC(I))*R,aimag(FRMC(I))*R
           WWW = sign(abs(FRMC(I)),real(FRMC(I)))  ! abs only, with sign of real part
	   endif
156		format(f8.3,2g15.5)
           W2 = WWW**2 * R*R * HDN
         BEE(KN,2) = BEE(KN,2) + W2
         BEE(KN,3) = BEE(KN,3) + W2 * R*R
         IF(EIGEN.and.((mod(IPC,2).eq.1.and.I.ge.10)
     x             .or.I.eq.MINT-1)) THEN
             	IE = 0
       	    WWW = WWW*R
             	CALL WHIT(ETA,R,FK,E,L,WL,WLD,IE)
		T = exp(dble(IE))
           	ANC = WWW*T/WL(L+1)
	       if(final) then
		call openif(46)
	   	RHOA = R*FK
	   	EX = ETA*LOG(2.0*RHOA)+RHOA
             	write(46,59) R,ANC,WWW,WL(L+1)/T,WLD(L+1)/WL(L+1),IE,
     x 		   exp(-EX)
		written(46) = .true.
59       	format(f9.3,4g15.6,i4,g15.6)
	       endif
         ENDIF
60    CONTINUE
	    if(mod(IPC,2)>0.and.final) write(58,*) '&'
	if(mod(IPC,2)>0.and.NFL.le.0.and.PWF) then
	   if(rlb)  write(55,64) (PSI(I,IC),I=1,MINT)
	   if(cmb) write(55,64) (PSIC(I,IC),I=1,MINT)
	   call flush(55)
	   endif
        BEE(KN,4) = ANC
         IF(mod(IPC,2).eq.1.and.EIGEN.and.final) then
	   RHOA = R*FK
	   T = ETA*LOG(2.0*RHOA)+RHOA
            write(KO,61) R,KN,ANC,ETA,T
61           format(/'    At R=',F9.3,', ANC in ch',I3,' is',g14.5,
     x          ',   eta,eta*ln(2*rho)+rho =',2g14.5)
           write(46,*)  '&'
         endif
      IF(rlb.and.NFL.LE.0) FORMR(1) = 2. * FORMR(2) - FORMR(3)
      IF(cmb.and.NFL.LE.0) FRMC(1) = 2. * FRMC(2) - FRMC(3)
C     if(abs(be(kn,2)).lt.1e-20) go to 860
      KM = KM + 1
      IF(KM.NE.KN) THEN
         DO 62 I=1,17
62       QNF(I,KM) = QNF(I,KN)
         DO 63 I=1,4
63       BEE(KM,I) = BEE(KN,I)
          AFRAC(ITC(ICP,IB),ITC(ICR,QNF(3,KN)),IN,KN) = 0.0
       ENDIF
       IF(NFL.LT.0.and.nlag==0) THEN
          WRITE(ABS(NFL),635) KM,IC,M,NAME(IN,ICR),NAME(IN,ICP),
     x       qnf(9,km),qnf(10,km)*0.5,qnf(11,km)*.5,
     x       JEX(QNF(6,KM),QNF(2,KM),QNF(3,KM)),JCOM,-BEE(KM,1)
635	  format(' Overlap wf',i3,'#',i3,'/',i3,' of <',a8,'|',a8,'>:',
     x       ' lsj,I;J =',i3,2f5.1,',',f5.1,';',f5.1,' @',f9.5,' MeV')
          WRITE(ABS(NFL),*) MINT,real(HDN),0
          if(rlb) WRITE(ABS(NFL),64) (FORMR(I),I=1,MINT)
          if(cmb) WRITE(ABS(NFL),64) (FRMC(I),I=1,MINT)
	  written(abs(NFL)) = .true.
         ENDIF
64     FORMAT(1P,6E12.4)
         if(rlb) WRITE(8,REC=KM) (FORMR(I),I=1,N)
         if(cmb) WRITE(8,REC=KM) (FRMC(I),I=1,N)
         WNORM = WNORM + BEE(KM,2)
         WRMS  = WRMS  + BEE(KM,3)
         BEE(KM,3)= SQRT(BEE(KM,3)/(BEE(KM,2) + SMALL))
         BEE(KM,2) = SQRT(BEE(KM,2))
         IF(EIGEN) BEE(KM,2) = BEE(KM,2) * SIGN(1D0, FORMR(3))
         QNF(8,KM) = 1
       DO 65 I=1,MINT
       if(rlb) FORMR(I) = FORMR(I)/(MAX(I-1.,0.01)*HDN)**QNF(9,KM)
       if(cmb) FRMC(I) = FRMC(I)/(MAX(I-1.,0.01)*HDN)**QNF(9,KM)
       IF(I.GE.3.and.rlb) then
          if(FORMR(I)*FORMR(I-1).LT.0..and.I<MINT)  then
	  	QNF(8,KM) = QNF(8,KM) + 1
!		if(IPC>2) write(KO,*) ' Ch ',KM,' node at ',I
		endif
        endif
65       CONTINUE

        if(rlb) FORMR(1) = 2.0 * FORMR(2) - FORMR(3)
        if(cmb) FRMC(1) = 2.0 * FRMC(2) - FRMC(3)
c
	if(rlb) then   !  put into FRMC in any case
	 do I=1,MINT
	 FRMC(I) = FORMR(I)
	 enddo
	endif
       DO 70 I=1,NLN
          R = (I-1)/RIN / HDN
	TC = FFC4(R,FRMC,MINT)
       if(.not.cxwf) FORML(I,KM,1) = TC   ! real part only
70     if(     cxwf) FORMC(I,KM,1) = TC
       FORMR(:) = 0.0
       FRMC(:) = 0.0
C
      IF(NFL.LE.0 .AND. KIND.NE.4.and.nlag==0) THEN  ! get VERTEX functions
      DO 86 JF=1,NF
      DO 86 IP=1,M
      IF(KRPOT.GT.0) THEN
         T = - CC(IC,IP,JF,3)/CONV
       ELSE
         T = - (CC(IC,IP,JF,1) + VARY*CC(IC,IP,JF,2)) / CONV
      ENDIF
      IF(ABS(T).LT.SMALL) GO TO 86
         DO 80 I=2,MINT
         R = (I-1)*HDN
         if(rlb) FORMR(I) = FORMR(I)+PSI(I,IP)*FORMF(I,JF)*T/R
         if(cmb) FRMC(I)  = FRMC(I)+PSIC(I,IP)*FORMF(I,JF)*T/R
80       CONTINUE
86    CONTINUE
       ELSE IF(NFL.LE.0 .AND. KIND.EQ.4.and.nlag==0) THEN
         DO 87 I=2,MINT-1
         R = (I-1)*HDN
         FORMR(I) =-(PSI(I+1,1) - 2*PSI(I,1) + PSI(I-1,1)) /(HDN*HDN *R)
 87      CONTINUE
       ELSE
C             NFL GT 0 .or. nlag>0
          READ(NFL,*) (FORMIN(I),I=1,NPOINTS)
!	 write(301,'(f8.3,g12.4)') (RFIRST+(I-1)*RSTEP,FORMIN(I),
!     x              I=1,NPOINTS)
!	 write(301,*) '&'
          DO I=1,MINT
            R = (I-1)*HDN
            FORMR(I) = FF4(R-RFIRST,RSTEPI,FORMIN,NPOINTS)
	  enddo
       ENDIF
       IF(NFL.LT.0.and.nlag==0) THEN  
           if(rlb) then
              FORMR(1) = 2.0 * FORMR(2) - FORMR(3)
	      WRITE(ABS(NFL),64) (FORMR(I),I=1,MINT)
	      endif
           if(cmb) then
              FRMC(1) = 2.0 * FRMC(2) - FRMC(3)
	      WRITE(ABS(NFL),64) (FRMC(I),I=1,MINT)
	      endif
	   written(abs(NFL)) = .true.
        ENDIF
         D0(KM,1) = 0.0
         D0(KM,2) = 0.0
      DO 88 I=1,MINT
         R = MAX(I-1.,1E-5)*HDN
         R1 = R ** QNF(9,KM)
         if(rlb) TC = FORMR(I)
         if(cmb) TC = FRMC(I)
	    if(mod(IPC,2)>0.and.PWF.and.rlb) write(59,156) R,FORMR(I)*R
	    if(mod(IPC,2)>0.and.PWF.and.cmb) write(59,156) R,FRMC(I)*R
         D0(KM,1) = D0(KM,1) + R * TC *         R  * R1
         D0(KM,2) = D0(KM,2) + R * TC * SINH(FK*R) * R1 / FK
         IF(rlb.and.R1.NE.0.0) FORMR(I) = FORMR(I) / R1
         IF(cmb.and.R1.NE.0.0) FRMC(I) = FRMC(I) / R1
88	continue
	    if(mod(IPC,2)>0.and.PWF) write(59,*) '&'
c      
      IF(rlb.and.QNF(9,KM).GT.0) FORMR(1) = 2.0 * FORMR(2) - FORMR(3)
      IF(cmb.and.QNF(9,KM).GT.0) FRMC(1) = 2.0 * FRMC(2) - FRMC(3)
c
         R = HDN * SQRT(4*PI)
            DO 89 I=1,QNF(9,KM)
89          R = R / ((2.*I+1.)*(2.*I))
         D0(KM,1) = D0(KM,1) * R
         D0(KM,2) = D0(KM,2) * R
            WD0 = WD0 + D0(KM,1)
            WD  = WD  + D0(KM,2)
       DO 90 I=1,NLN
          R = (I-1)/RIN / HDN
          if(rlb) TC = FFR4(R,FORMR,MINT)
          if(cmb) TC = FFC4(R,FRMC,MINT)
       if(.not.cxwf) FORML(I,KM,2) = TC
       if(     cxwf) FORMC(I,KM,2) = TC
            R = MAX(I-1.,1E-5)/RIN
            R1 = R ** (QNF(9,KM)+1)
      IF(IPC.GE.6) THEN
          if(I==1) WRITE(132,*) 'Vertex function @',KM,QNF(9,KM)
       if(.not.cxwf) write(132,'(f8.3,3g12.4)') R,FORML(I,KM,:)*R1
     x           , FORML(I,KM,2)/FORML(I,KM,1)
       if(     cxwf) write(132,'(f8.3,6g12.4)') R,FORMC(I,KM,:)*R1
     x           , FORMC(I,KM,2)/FORMC(I,KM,1)
      endif
90     continue
      IF(IPC.GE.6) THEN
          WRITE(KO,992) KM,QNF(9,KM)
          if(.not.cxwf) WRITE(KO,997) ((I-1)/RIN,
     x         ((I-1)/RIN)**(QNF(9,KM)+1)*FORML(I,KM,2),I=1,NLN)
          if(     cxwf) WRITE(KO,998) ((I-1)/RIN,
     x         ((I-1)/RIN)**(QNF(9,KM)+1)*FORMC(I,KM,2),I=1,NLN)
         ENDIF
! 992 FORMAT(/' V.Psi/R**L @ channel #',I3,' with L =',I2,
!    &       //'  R   V.Psi/R**L')
  992 FORMAT(/' V.Psi @ channel #',I3,' with L =',I2,
     &       //'  R   V.Psi')
  997 FORMAT(4(0p,F6.2 ,1p,E12.4,2X))
  998 FORMAT(3(0p,F6.2 ,1p,2E12.4,2X))
C
	CH2 = ' '
	if(M>1.and.IC==IL) CH2 = '*'
        IF(NFL.GT.0) VARY=0.
      if(abs(BEE(KM,4))<99) then
      WRITE(KO,858) KM,ICR,ICP,IN,KIND,CH1,CH2,QNF(8,KM),QNF(9,KM),SN,
     &    QNF(3,KM),QNF(11,KM)*.5,QNF(5,KM),KBPOT,KRPOT,BEE(KM,1),
     &    ISC,IPC,NFL,AMPL,  ADJ,VARY,
     &    ABS(DZ),DM,FK,BEE(KM,2),BEE(KM,3),D0(KM,1),D0(KM,2),BEE(KM,4)
        else
      WRITE(KO,8581) KM,ICR,ICP,IN,KIND,CH1,CH2,QNF(8,KM),QNF(9,KM),SN,
     &    QNF(3,KM),QNF(11,KM)*.5,QNF(5,KM),KBPOT,KRPOT,BEE(KM,1),
     &    ISC,IPC,NFL,AMPL,  ADJ,VARY,
     &    ABS(DZ),DM,FK,BEE(KM,2),BEE(KM,3),D0(KM,1),D0(KM,2),BEE(KM,4)
        endif
858   FORMAT(/1x,I4,':',4I3,2X,2A1,i2,I3,F4.1,I3,F4.1,I4,I4,I4,F9.4,2I3,
     &  I3,F7.4,1X,A3,F7.4,';',F3.0,F6.3,F7.4,F7.4,F7.3,F7.1,F7.1,F8.4)
8581  FORMAT(/1x,I4,':',4I3,2X,2A1,i2,I3,F4.1,I3,F4.1,I4,I4,I4,F9.4,2I3,
     &  I3,F7.4,1X,A3,F7.4,';',F3.0,F6.3,F7.4,F7.4,F7.3,F7.1,F7.1,
     &  1p,E10.2)

      NSP = MAX(NSP,KM)
860   CONTINUE
         DO 861 KN=KM+1,KN2
         DO 861 I=1,16
861      QNF(I,KN) = 0
      IF(KN2.GT.KN1) WRITE(KO,859) WNORM,SQRT(WRMS/MAX(WNORM,Z+SMALL)),
     &                         WD0,WD
859   FORMAT(/' Overall rms norm & radius =',2F8.4,
     &     ', & Overall D0 & D =',2F9.2)
 	 if(allocated(PSI)) deallocate(PSI)
 	 if(allocated(PSIC)) deallocate(PSIC)

         IF(KM.GT.KN1) WRITE(KO,*) ('-',I=1,132)
      		LAST1 = .TRUE.
		do id=1,datasets
		if(data_type(id)==5 .or. data_type(id)==9) then
		 do ip=1,datalen(id)
 		  if(KN1==bs_no(ip,id)) then
		   if(data_type(id)==5) T = VARY
		   if(data_type(id)==9) T = BEE(KM,4)
	 	   EE = (T-datavals(ip,id))/dataerr(ip,id)
                   data_chisq(id) = data_chisq(id) + EE**2
                   theoryvals(ip,id) = T
		   endif
	          enddo
		  endif
	         enddo
      ELSE !if(kind.gt.5)
C
C    two-particle form factors
C    -------------------------
       if(NNN.lt.5) then
        write(6,*) ' Only ',NNN,' points for two-particle transfers!'
        stop 'RNN'
       endif
      IT = KBPOT
      KNZR = KRPOT
      EP2 = E
      LMIN = L
      SMAX = 1.0
      SMIN = SN
      JN = J
      IF(KIND.EQ.9) THEN
         JN = 0.0
         IB = MAX(IB,1)
         ENDIF
      NK = NN
      IF(ABS(EP2).LT.1E-5) EP2 = 1.0
      IF(KNP.EQ.1.OR.LAST1) WRITE(KO,8625)
      WRITE(KO,8621) KN1,KN2,IC1,IC2,INI,KIND,CH1,J,IT,NK,KNZR,LMIN,LMAX
     &       ,SMIN,NFL,ISC,IPC,EP2,AMPL
8621  FORMAT(/1X,6I4,3X,A1,F5.1,',',I1,3I5,I7  ,F5.1,3I4,F6.2,'%',f8.4)
8625  FORMAT( /' The following TWO-PARTICLE FORM FACTORS are constr'
     & ,'ucted :'//'  KN1 KN2  P1  P2  IN KIND TYP J12,T  PAIRS KNZR ',
     &  'LMIN LMAX SMIN FILE ISC IPC  EP2%  AMPL')
      TNT(1,1)=0
      FK= 0.
      IF(NK.GT.0) THEN
!            READ(KI,*) ((TNT(I,JJ),I=1,4),COEF(JJ),JJ=1,NK)
            READ(KI,nml=twont)
8626  FORMAT(3(4I3,F8.4))
	EGS = 0.0
      IF(TNT(1,1).GE.1) THEN
      EGS = BEE(TNT(1,1),1) + BEE(TNT(2,1),1)
         IF(KNZR.NE.0) EGS = EGS - BEE(KNZR,1)
      ENDIF
           WRITE(KO,8627) ((TNT(I,JJ),I=1,4),COEF(JJ),JJ=1,NK)
8627       FORMAT(/5(4I4,F8.4,','))
        ENDIF
	if(cxwf) then
	do 8629 KN=1,MSP
	do 8629 ip=1,2
	do 8629 I=1,NLN
8629	FORML(I,KN,ip) = FORMC(I,KN,ip)   ! to get real single-nucleon states
	endif
	if(IPC>0.and.TNT(1,1)>0) then
	KN = TNT(1,1)
	write(6,*) 'SP state :',KN
	write(6,'(10f9.4)') FORML(1:10,KN,1)
	endif

      CALL TWONN(KN1,KN2,HDN,MR,RNN,ISC,IPC,NLN,NNN,QNF,EGS,
     &       KIND,IN, IAK,IAMIN,IAMAX,IB,JEX,BAND,MXP,MXX,ICR,ICP,COPY,
     &       KNZR,LMIN,LMAX,SMIN,SMAX,JN,J,IT,D0,RIN,DM,MASS(IN,ICR),
     &   NK,TNT,COEF,NFL,cxwf,FORML,MAXNLN,MSP,FORMR,RMI,FMSCAL,
     &   NNN*2,NNU,BEE,FK,N,EP2,NAME(IN,ICP),MAX(NK,1))
      NSP = MAX(NSP,KN2)
      DO 867 KN=KN1,KN2
            IF(QNF(12,KN).EQ.1) WRITE(KO,865)
      WRITE(KO,865) KN,QNF(8,KN),QNF(9,KN),QNF(12,KN),QNF(13,KN),
     &      QNF(14,KN)*.5,QNF(10,KN)*.5,QNF(11,KN)*.5,(BEE(KN,I),I=1,3),
     &      (D0(KN,I),I=1,2)
865   FORMAT(27X,I4,':',' NN LL,(n l s)j; J12 =',2I3,',(',2I2,F4.1,')',
     &  F4.1,';',F4.1,'  BE,Norm =',F8.3,F8.4,',',F6.3,F8.1,F7.1)
      QNF(2,KN) = ICR
      QNF(4,KN) = ICP
      QNF(5,KN) = IB
      QNF(6,KN) = IN
      QNF(7,KN) = KIND
      IA = QNF(3,KN)
      IF(IA.GE.1 .AND. IB.GE.1)THEN
       AFRAC(ITC(ICP,IB),ITC(ICR,IA),IN,KN) = AMPL
      ENDIF
      QNF(16,KN) = ICHAR(CH1)
      QNF(17,KN) = KN2-KN1+1
	if(cxwf) then
	do 866 ip=1,2
	do 866 I=1,NLN
866	FORMC(I,KN,ip) = FORML(I,KN,ip)   ! to get back two-nucleon states
	endif
867   CONTINUE
      IF(IPC.GT.0) WRITE(KO,871)
871   FORMAT(' ',132('-'))
      LAST1 = .FALSE.
      ENDIF !if(kind<=5 or kind>5)

890     if(allocated(PSI)) deallocate (PSI)
        if(allocated(CC)) deallocate (CC)
        if(EIGEN.and.ISC<0) THEN ! renormalise all varied potentials permanently!
         DO JF=1,NF
         if(REN(JF)) then
           FORMF(:,JF) = FORMF(:,JF) * VARY
           ADJ = VAR(abs(ISC)+1)
           write(KO,895) ADJ,JF,VARY
895        format('  ### ',a3,' potential at ',i3,
     x            ' permanently scaled by ',f8.4)
           endif
         ENDDO
         ENDIF
!	 if(NFL/=0) write(6,*) 'REWIND ',abs(NFL)
900   CONTINUE
C
1000      inquire(55,opened=op)
	  if(op) close(55) 
	  inquire(46,opened=op)
	  if(op) close(46)
	  do i=1,1000
	  if(FUSED(i)) rewind i
	  enddo
      RETURN
      END
*****EIGCC**************************************************************
      SUBROUTINE EIGCC(PSI,FORMF,CCF,NFM,QNF,ETAP,AB,MR,KAP2,THETA,P,
     &  MAXN,MC,NODES,NP,H,M,PCON,MM2,MAXC,EPS,IFAIL,LMX1,MINT)
	use io
	use factorials
	use drier
	! NP here is MINT
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PSI(MINT,M),F(NP,M,M),CCF(M,M,NFM,3),
     &       KAP2(M),ETAP(M),MAT(MM2,MM2)
      COMPLEX*16 FORMF(MAXN,NFM)
      INTEGER QNF(19,M),PCON,COUNT,BC,CCP,CC
      REAL*8 CENT(M),ZI(M,M),ZM(M,M),ZP(M,M)
     &             ,COUT(M),COUTP(M),COUPL(M,M)
      REAL*8 ETA,K,WL(LMX1),WLD(LMX1)
      LOGICAL SING
C
      COUNT = 1
      CCP = 1
      BC = 0
      PP=P
      N = NP-1
      RN = (NP-1)*H
      HP = H*H
      R12= 1./12.
      HP12=HP/12.
      NR = 2*M+1
      IFAIL = 0
      SMALL = 1D0 / FPMAX
      SMALLR = SQRT(SMALL)
      SMALLQ = SQRT(SMALLR)
C     IF(NP.GT.MAXN.OR.NR.GT.MM2) STOP 101
      IF(NP.GT.MAXN.OR.NR.GT.MM2) CALL ABEND(8)
      IF(PCON.GE.3) WRITE(KO,207) NODES,MC
  207 FORMAT(/'    Parameters      Mismatch   Nodes ->',I2,' in ch',i3/)
102   DO 21 J=1,M
      CENT(J) = QNF(9,J) * (QNF(9,J)+1)
      K = SQRT(ABS(KAP2(J) + THETA*P)) * AB
      ETA  = 0.5*ETAP(J)/K
      L = QNF(9,J)
      IE = 0
      CALL WHIT(ETA,RN+H,K,E,L,WL,WLD,IE)
      COUTP(J) = WL(L+1)
      CALL WHIT(ETA,RN,K,E,L,WL,WLD,IE)
21    COUT(J) = WL(L+1)
      I   =(NP + 1) *3/4
         IMAX = NP/2
         DEL = -SQRT(FPMAX)
22    I=I-1
      UDIAG = -KAP2(MC) - THETA*P
      DO 23 JF=1,NFM
23    UDIAG = UDIAG + DBLE(FORMF(I,JF)) *
     &           ( CCF(MC,MC,JF,1) + P * CCF(MC,MC,JF,2))
      DEN = -CENT(MC)/((I-1)*H)**2 + UDIAG
         IF(DEL.LT.DEN) THEN
            DEL = DEN
            IMAX = I
            ENDIF
      IF(DEN.LT.0 .AND. I.GT.10) GO TO 22
      MAM = I
      IF(I.EQ.10) MAM = IMAX
      MAP = MAM+1
      IF(COUNT.EQ.1 .AND. PCON.GE.7)
     &WRITE(KO,*)MC,QNF(9,MC),K,ETA,RN+H,COUTP(MC),RN,COUT(MC),MAM,MAP
     &             ,CENT(MC)
103   DO 30 J=1,M
      DO 25 IT=1,M
      ZM(IT,J) = 0.0
25    ZI(IT,J) = 0.0
      ZM(J,J) = COUTP(J)
30    ZI(J,J) = COUT(J)
C
C      outer integration,  zi from np to map
      NT = -1
      NF = NP
      NO = NT * (MAP-NP) + 1
40    DO 60 III=1,NO
         I = NF + (III-1)*NT
         II= (I + MIN0(NT,0)) + 1
         RRI= 1.0/(II-1.)**2
      DO 42 IT=1,M
      DO 42 J=1,M
42    F(I,J,IT) = ZI(IT,J) * (1. + CENT(J)*RRI*R12 )
       DO 45 J=1,M
       DO 45 L=1,M
         C = 0.0
         DO 425 JF=1,NFM
            T = CCF(L,J,JF,1)  + P * CCF(L,J,JF,2)
            IF(T.EQ.0.0) GO TO 425
         C = C + T * DBLE(FORMF(II,JF))
425      CONTINUE
         IF(L.EQ.J) C = C - KAP2(J) - THETA * P
       COUPL(L,J) = C
       IF(C.EQ.0) GO TO 44
       C = C * HP12
          DO 43 IT=1,M
43         F(I,L,IT) = F(I,L,IT) - C * ZI(IT,J)
44    CONTINUE
45    CONTINUE
      DO 49 IT=1,M
         DO 49 L=1,M
49       MAT(IT,L) = 0.0
      DO 54 L=1,M
      DO 54 J=1,M
      C = COUPL(L,J)
      IF(C.EQ.0.0) GO TO 54
      C = C * HP
      DO 53 IT=1,M
53    MAT(IT,L) = MAT(IT,L) + C * F(I,J,IT)
54    CONTINUE
      DO 55 J=1,M
      DO 55 IT=1,M
      ZP(IT,J) = 2*ZI(IT,J) - ZM(IT,J) - MAT(IT,J)
     &                      + F(I,J,IT) * CENT(J) * RRI
      ZM(IT,J) = ZI(IT,J)
55    ZI(IT,J) = ZP(IT,J)
C                              now check for incipient overflows:
      DO 59 IT=1,M
      C = 0.
      DO 57 J=1,M
57    C = MAX(C,ABS(F(I,J,IT)))
      IF(C .LT. FPMAX) GO TO 59
      C = SMALLR
      DO 58 J=1,M
         ZI(IT,J) = ZI(IT,J) * C
         ZM(IT,J) = ZM(IT,J) * C
         DO 58 L=1,III
            I = NF + (L-1)*NT
58       F(I,J,IT) = F(I,J,IT) * C
59    CONTINUE
60    CONTINUE
C
      NT = -NT
      IF(NT.LE.0) GO TO 3
      DO 65 J=1,M
      DO 64 IT=1,M
      ZM(IT,J) = 0.0
64    ZI(IT,J) = 0.0
65    ZI(J,J) = 1E-10 * H**(QNF(9,J)+1) / EXP(0.5 * FACT(QNF(9,J)+1))
C      inner integration,  zi from 1 to mam
      NF = 1
      NO = NT * (MAM-1) + 1
      GO TO 40
3     COUNT = COUNT + 1
C   now calc. derivatives at matching pt. inner(zm) & outer(zp)
      DO 80 J=1,M
         MAT(J,NR) = 0.0
         MAT(J+M,NR) = 0.0
      DO 75 IT=1,M
      DEL =147.0*F(MAM,J,IT)-360.0*F(MAM-1,J,IT)+450.0*F(MAM-2,J,IT)
     1 -400.0*F(MAM-3,J,IT)+225.0*F(MAM-4,J,IT)
     2 -72*F(MAM-5,J,IT)+10.0*F(MAM-6,J,IT)
      ZM(IT,J) = 1  * DEL / (60.0 * H)
      DEL =147.0*F(MAP,J,IT)-360.0*F(MAP+1,J,IT)+450.0*F(MAP+2,J,IT)
     1 -400.0*F(MAP+3,J,IT)+225.0*F(MAP+4,J,IT)
     2 -72*F(MAP+5,J,IT)+10.0*F(MAP+6,J,IT)
      ZP(IT,J) =(-1)* DEL / (60.0 * H)
C     WRITE(KO,*) j,it,zm(it,j)/f(mam,j,it)-zp(it,j)/f(map,j,it)
      MAT(J,IT) = F(MAM,J,IT)
      MAT(J,IT+M) = -F(MAP,J,IT)
      MAT(J+M,IT) = ZM(IT,J)
75    MAT(J+M,IT+M) = - ZP(IT,J)
80    CONTINUE
      DO 78 IT=1,MM2
78    MAT(M+MC,IT) = 0.0
      MAT(M+MC,M+MC) = 1.0
      MAT(M+MC,NR)  = 1.0
c     do 81 j=1,nr
c     WRITE(KO,*) j,(mat(j,it),it=1,nr)
c81        continue
       CALL GAUSSR(2*M,MM2,MAT,SING,DET,SMALL,.false.)
       IF(SING) GO TO 200
       DO 88 L=1,M
          DO 85 I=1,NP
85        PSI(I,L) = 0.0
       DO 88 IT=1,M
       DO 87 I=2,MAM
87     PSI(I,L) = PSI(I,L) + F(I-1,L,IT) * MAT(IT,NR)
       DO 88 I=MAP,NP
88     PSI(I,L) = PSI(I,L) + F(I,L,IT) * MAT(IT+M,NR)
       DERIN = 0.0
       DEROUT= 0.0
       DO 90 IT=1,M
       DERIN = DERIN + ZM(IT,MC) * MAT(IT,NR)
       DEROUT= DEROUT+ ZP(IT,MC) * MAT(IT+M,NR)
90     CONTINUE
       NCO = 1
       DO 91 I=3,NP
C 91     IF(PSI(I,MC)*PSI(I-1,MC).LT.0.) NCO = NCO + 1
       IF(PSI(I,MC).GT.0 .AND. PSI(I-1,MC).LT.0.) NCO = NCO + 1
91     IF(PSI(I,MC).LT.0 .AND. PSI(I-1,MC).GT.0.) NCO = NCO + 1
       DEN = 0.0
       DO  96 J=1,M
       DO  96 L=1,M
         DO 94 JF=1,NFM
         T = CCF(J,L,JF,2)
         IF(ABS(T).LT.SMALL) GO TO 94
            DO 92 IMAX=NP,1,-1
C92          IF(ABS(DBLE(FORMF(IMAX,JF))).GT.SMALL) GO TO 925
92          IF(ABS(DBLE(FORMF(IMAX,JF))).GT.SMALLR .AND.
     X         ABS(PSI(IMAX,L)).GT.SMALLQ) GO TO 925
925         DO  93 I=2,IMAX
93          DEN = DEN + PSI(I,J) * T * DBLE(FORMF(I,JF)) * PSI(I,L)
94       CONTINUE
         IF(THETA.EQ.0.0 .OR. J.NE.L) GO TO 96
         DO 95 I=2,NP
95       DEN = DEN - PSI(I,J) * THETA * PSI(I,L)
96    CONTINUE
       DEL =(DEROUT - DERIN) / PSI(MAP,MC)
       IF(PCON.GE.3) WRITE(KO,217) P,DEL,NCO,MAM
C 217    format(e14.6,e15.3,i6)
  217    FORMAT(1X,F13.6,F15.8,2I6)
       POLD = P
       IF(COUNT.GT.MAXC) GO TO 210
       IF(abs(P).lt.1d-4) GO TO 210
       IF(NCO.NE.NODES) GO TO 6
       IF(ABS(DEL).LT.EPS) GO TO 7
         IF(DEN.EQ.0.0) GO TO 220
       Q =-DEL * PSI(MAP,MC)**2 / (DEN * H)
       IF(P.LT.Q) Q=P
       IF(P.LT.-Q) Q = -0.5*P
       P = P+Q
       BC = 1
       Q = 2*Q
99     IF(THETA.NE.0.0) GO TO 102
       if(P>POLD*1.8) GO TO 102  ! recalculate matching point MAM for big potl changes!
       if(P<POLD*0.6) GO TO 102  ! recalculate matching point MAM for big potl changes!
       GO TO 103
6     IF(BC.NE.0) GO TO 61
      CC = 1
      IF(NODES.LT.NCO) CC = -1
      IF(DEN.LT.0.0) CC = -CC
      IF(CC.LT.0) CCP = -IABS(CCP)
      IF(CCP.GT.0) PP = P
      IF(CCP.LT.0) PP = 0.5*PP
      P = PP*CC + P
      GO TO 99
61    Q = 0.5*Q
      P = P - 0.5*Q
      GO TO 99
7     Q = 0.0
      DO 700 J=1,M
      DO 700 I=1,NP
700   Q = Q + PSI(I,J)**2
         PIC = 1.0
      Q = SIGN(PIC,PSI(3,MC)) / SQRT(Q * H)
      DO 710 J=1,M
      DO 710 I=1,NP
710   PSI(I,J) = PSI(I,J) * Q
      IF(PCON.EQ.0) RETURN
      IF(MOD(PCON,2).EQ.0)      RETURN
      DO 996 J=1,M
      WRITE(KO,992) J,QNF(9,J),QNF(11,J)*0.5,QNF(3,J)
  992 FORMAT(/' For Channel #',I3,' with l =',I2,',  j =',F5.1,
     &  ', around core state #',I3,             //,'   R    Y(R)')
  996 WRITE(KO,998) ((I-1)*H,PSI(I,J),I=1,NP,MR)
  998 FORMAT(6(1X,0P,F6.2,2X,1P,E11.4,', '))
      RETURN
200   IFAIL = 1
      WRITE(KO,202) DET
202   FORMAT(' ***** FAIL IN EIGCC : DETERMINANT =',1P,2E12.4)
      STOP
210   IFAIL = 2
      WRITE(KO,212) DEL,MAXC
212   FORMAT(' ***** FAIL IN EIGCC : DISCREPANCY =',1P,E12.4,
     &       ' EVEN AFTER',I3,' ITERATIONS'/)
      STOP
220   IFAIL = 3
      WRITE(KO,222)
222   FORMAT(' ***** FAIL IN EIGCC : VARIABLE POTENTIAL PART ZERO'/)
      STOP
      END
      SUBROUTINE GAUSSR(N,NR,A,SING,DET,EPS,SHOW)
C
C    solve by Gaussian elimination sum(j): A(i,j).P(j) = A(i,N+1)
C             where A(j) is left in A(j,N+1)
C
	use io
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(M = 1)
       REAL*8 A(NR,N+M),DET,RA
       LOGICAL SING,SHOW
       SING  = .FALSE.
      NPM = N + M
C      DO 201 I=1,N
C201   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
      DET = 1
      DO 9 K = 1, N
      DET = DET * A(K,K)
      IF (ABS (A(K,K)) .GT. EPS ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET,A(K,K)
    3    FORMAT(//' THE MATRIX IS SINGULAR AT',I3,'  determinant is ',
     &  2E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
402   FORMAT( 1X,20F6.3/(1X,20F6.3))
         RETURN
    5 KP1 = K + 1
         RA = 1.0/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = 1
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
Cdir$ ivdep
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = 0
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET
15    FORMAT(/' The determinant is ',E16.8)
      RETURN
      END
*****BINCC*************************************************************
      SUBROUTINE BINCC(Y,FORMF,CC,NF,M,K2,IL,CONV,BPHASE,ISC,KMIN,KMAX,
     & NK,QNF,ETAP,NMAX,H,PCON,TRES,TDEL,TKMAT,TKNRM,LMX1,MAXN,ANC)
	use factorials
	use io
	use drier
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KMIN,KMAX,K,KP,INTP,K2(M),CC(M,M,NF,3),BPHASE(NK),HSP(M)
      INTEGER PCON,QNF(19,M),IOP(M)
      LOGICAL SING,TRES,TKNRM,TDEL,TKMAT,TRA
      COMPLEX*16 F(NMAX+1,M,M),CH,CHD,TMAT,TMATI,YFAC,
     &           ZI(M,M),ZM(M,M),ZP(M,M),
     &           MAT(2*M,3*M),C,T,COUPL(M,M),CI
      COMPLEX*16 FORMF(MAXN,NF),Y(NMAX+1,M),W(NMAX+1,M),WY
      REAL*8 ETA,KJ,CF(LMX1),CG(LMX1),CFP(1),CGP(1),KEIG,ETAP(M)
C
C     pcon    trace iterations     final iteration    list wave function
C      0          no                    no               no
C      1          no                    yes              yes
C      2          no                    yes              no
C      3          yes                   yes              yes
C      4          yes                   yes              no
C
C  If TDEL, multiply wave functions by phase of T before integrating
C  If TRES, then multiply wave functions by abs(T) before integrating.
C  If TKNRM, then multiply wave functions by k before integrating.
C
      NMAXP=NMAX+1
	MAXNR = 2*M
	MAXNR1 = 3*M
      SMALL = 1D0/FPMAX
         CALL CHECK(NMAXP,MAXN,7)
	TRA = PCON>4 .and. TKMAT
      HP=H*H
      R12 = 1./12.
      CI = (0.,1.)
      SQFPI = SQRT(4.*PI)
      RADEG = 180.0/PI
      HP12 = HP*R12
       II = MIN(NK/10,1)
      MATCH = NMAX-10
      DK = (KMAX - KMIN)/(NK-1)
C     INTP = DK/SQRT(KMAX-KMIN) * SQRT(2.0/PI)
      INTP = DK * SQRT(2.0/PI)
      NP = 2*M + 1
       DO 12 L=1,M
         DO 10 N=1,NMAXP
10       Y(N,L) = 0.0
12     CONTINUE
        YINT = 0.0
        IF(PCON.GE.5) then
       WRITE(KO,998)
       write(42,*) DK,NK,KMIN
       write(42,*) H,NMAXP,0
       write(42,*) M,(QNF(9,L),l=1,M)
       endif
C
	ANC = 1e+6
	DELTAL = 0.0
	deltalast = 0.
	deltaiadd = 0.
	ELAST = 0.0
      DO 90 IK=1,NK
         K = KMIN + (IK-1)*DK
         KP= K*K
         NOP = 0
         do 20 J=1,M
         TKJ = K2(IL) + KP-K2(J)
         if(TKJ>0.) then
            NOP = NOP+1
            IOP(NOP) = J
            endif
20	  continue
C
      DO 30 L=1,M
          X = (K2(IL) + KP-K2(L))/CONV
 	  if(Abs(X).lt.0.002) then
 		write(KO,15) KP/CONV,L,X
 15		format(' Skipping energy',f8.4,' as ch',i3,' is ',
     X    ' at energy',f8.4,': too close to threshold!')
    		go to 90
 		endif
      DO 28 J=1,M
      ZI(J,L) = 0.0
28    ZM(J,L) = 0.0
30    ZI(L,L) = H**(QNF(9,L)+1) / EXP(0.5 * FACT(QNF(9,L)+1))
C
      F(1,:,:) = 0d0
      DO 60 I=2,NMAXP
         RRI= 1.0/(I-1.)**2
      DO 42 IT=1,M
      DO 42 J=1,M
42    F(I,J,IT) = ZI(IT,J) * (1. + QNF(9,J)*(QNF(9,J)+1)*RRI*R12 )
       DO 45 J=1,M
       DO 45 L=1,M
         C = 0.0
         DO 425 JF=1,NF
            T = CC(L,J,JF,1)
            IF(ABS(T).LT.SMALL) GO TO 425
         C = C + T * FORMF(I,JF)
425      CONTINUE
         IF(L.EQ.J) C = C + K2(IL) + KP-K2(J)
       COUPL(L,J) = C
       IF(C.EQ.0) GO TO 44
       C = C * HP12
          DO 43 IT=1,M
43         F(I,L,IT) = F(I,L,IT) - C * ZI(IT,J)
44    CONTINUE
45    CONTINUE
      DO 49 IT=1,M
         DO 49 L=1,M
49       MAT(IT,L) = 0.0
      DO 54 L=1,M
      DO 54 J=1,M
      C = COUPL(L,J)
      IF(C.EQ.0.0) GO TO 54
      C = C * HP
      DO 53 IT=1,M
53    MAT(IT,L) = MAT(IT,L) + C * F(I,J,IT)
54    CONTINUE
      DO 55 J=1,M
      DO 55 IT=1,M
      ZP(IT,J) = 2*ZI(IT,J) - ZM(IT,J) - MAT(IT,J)
     &                      + F(I,J,IT) * QNF(9,J)*(QNF(9,J)+1) * RRI
      ZM(IT,J) = ZI(IT,J)
55    ZI(IT,J) = ZP(IT,J)
60    CONTINUE
      X = (NMAXP-1) * H
      DO 65 J=1,MAXNR1
      DO 65 L=1,MAXNR
65       MAT(L,J) = 0.0
      DO 75 J=1,M
         DO 70 IT=1,M
            MAT(J,IT) = F(NMAXP,J,IT)
70          MAT(J+M,IT)=F(NMAX ,J,IT)
         TKJ = K2(IL) + KP-K2(J)
         KJ = SQRT(ABS(TKJ))
         L= QNF(9,J)
         XL= L
         ETA = ETAP(J) * 0.5 / KJ
         IF(TKJ.GT.0.0) THEN
            CALL COULFG(KJ*X,ETA,0D0,XL,CF,CG,CFP,CGP,2,0,I,M1)
             CH  = CMPLX(CG(L+1),CF(L+1)) * (0.,.5)
            CALL COULFG(KJ*(X-H),ETA,0D0,XL,CF,CG,CFP,CGP,2,0,I,M1)
             CHD = CMPLX(CG(L+1),CF(L+1)) * (.0,.5)
             HSP(J) = -atan2(CF(L+1),CG(L+1))*RADEG
         ELSE 
         IE = 0
           CALL WHIT(ETA,X,KJ,E,L,CF,CG,IE)
           CH = CF(L+1) * (0.,.5)
           CALL WHIT(ETA,X-H,KJ,E,L,CF,CG,IE)
           CHD = CF(L+1) * (0.,.5)
         ENDIF
!     IF(PCON.GE.5) WRITE(142,*) J,TKJ/CONV,KJ,X,L,ETA,CH,CHD
!    x       			,MAT(J,J),MAT(J+M,J)
          MAT(J,J+M) = CH
          MAT(J+M,J+M)=CHD
!          IF(J.EQ.IL) THEN
!             MAT(J,NP) = - CONJG(CH)
!             MAT(J+M,NP)=- CONJG(CHD)
!             ENDIF
	   MAT(J,2*M+J) = - CONJG(CH)
	   MAT(J+M,2*M+J)=- CONJG(CHD)
75      CONTINUE
C
       CALL GAUSS5(2*M,MAXNR,MAT,M,SING,T,SMALL,.FALSE.)
C         IF(SING) STOP 'SINGULAR BIN'
          IF(SING) THEN
             WRITE(KO,*) 'SINGULAR ENERGY SOLUTION OMITTED FROM BIN,'
             WRITE(KO,77) ((K2(IL) + KP-K2(J))/CONV,J=1,MIN(M,5))
77           FORMAT(' at channel energies =',5F10.4)
             GO TO 90
           ENDIF
C
C**********************************************************************
C*************Inserted code for core ex********************************
C**********************************************************************
C
      W(:,:) = 0.
      IOUT = IL
      IN = IOP(IOUT)
      NP = 2*M + IN

      DO 810 J=1,M		! so W = scattering wf at 1 energy
      DO 810 IT=1,M
      DO 810 N=1,NMAXP
      W(N,J) = W(N,J) + F(N,J,IT) * MAT(IT,NP)  
810    continue    

      TMATI = (MAT(IN+M,NP)-1.)/(0.,2.)
      DELTAI = 0.0
       X = 0.0
           IF(abs(MAT(IN+M,NP)).gt.1e-20) then
             X = (0.,-.5) * LOG(MAT(IN+M,NP))
             DELTAI = RADEG * X
            endif
      WY = DELTAI
	 if(IK>1) then
           if(deltai<deltalast-90) deltaiadd=deltaiadd+180
           if(deltai>deltalast+90) deltaiadd=deltaiadd-180
	 endif
          deltalast=deltai
          deltai = deltai + deltaiadd
	
	BPHASE(IK) = DELTAI/RADEG
      IF(PCON.GE.3)  WRITE(KO,999) IK,IN,K,KP/conv,
     X   MAT(il+M,NP), DELTAI,-sqfpi*tmati / conv

      YFAC = 1.0
        IF(TDEL) YFAC = EXP(-CI*BPHASE(IK))
        IF(TRES) YFAC = CONJG(TMATI)
        IF(TKNRM) YFAC = YFAC*K
      YINT = YINT + ABS(YFAC)**2 * DK
C
      DO 815 J=1,M
      DO 815 N=1,NMAXP
      Y(N,J) = Y(N,J) + W(N,J) * YFAC * INTP
815    continue
!	 N = nint(10./H) + 1
!	write(301,816) KP/conv,W(N,1),YFAC*INTP/DK,DELTAI
!816     format(f10.5,2f10.5,2f8.4,',',f8.4)

       DERIV = 1e-6
       if(IK>1) DERIV = (DELTAI-DELTAL)/(RADEG * (KP/conv - ELAST))
	 DELTAL = DELTAI
	 ELAST =  KP/conv
	 ANC = min(ANC,2d0/DERIV)
      IF(PCON.GE.3)  then
	 WRITE(43,9997) K,WY,KP/conv,-tan(DELTAI/RADEG)/K
         ETA = ETAP(1) * 0.5 / KJ
	 PHRES = DELTAI-HSP(IL)
	 STR = SIN(PHRES/RADEG)
         if(IK<2) WRITE(44,9998) KP/conv,DELTAI,0d0,ETA,
     x                           HSP(IL),PHRES,STR
         if(IK>1) WRITE(44,9998) KP/conv,DELTAI,2d0/DERIV,ETA,
     x                           HSP(IL),PHRES,STR
         WRITE(144,'(i10,1p,e12.4)') IK,DELTAI
	 written(43) = .true.
	 written(44) = .true.
	 endif
C
      IF(PCON.GE.5) then			! for information only
      DO 80 J=1,M
         BE = (K2(IL) + KP-K2(J))/CONV
         KJ = SQRT(ABS(BE)*CONV)
       IN = 0
       IF(J.eq.IL) IN =1
      TMAT = (MAT(J+M,NP)-IN)/(0.,2.)
      DELTA = 0.0
      IF(abs(MAT(J+M,NP)).gt.1e-20)
     .DELTA = RADEG * (0.,-.5) * LOG(MAT(J+M,NP))
      IF(DELTA.LT.0.0) DELTA = DELTA + 180.0
      C =-SQFPI * TMAT / CONV
      IF(IK/II*II.EQ.IK) THEN
       if(J.eq.IL) then
          IF(BE.GT.0) WRITE(KO,999) IK,J,K,BE,MAT(J+M,NP),DELTA,C
     x    , real(Y(10,J)),0. , W(10,J) , YFAC , INTP, YINT
          IF(BE.LE.0) WRITE(KO,9999) IK,J,K,BE,MAT(J+M,NP)
       else
          IF(BE.GT.0) WRITE(KO,999) IK,J,K,BE,MAT(J+M,NP)
          IF(BE.LE.0) WRITE(KO,9999) IK,J,K,BE,MAT(J+M,NP)
       endif
      ENDIF
 80   continue
      endif

90    CONTINUE
      IF(PCON.ge.3) then
        write(43,*) '&'
        write(44,*) '&'
      endif
      YINT = 1.0/SQRT(YINT)
C
      WOLD = -1.
      IF(MOD(ISC,2).EQ.1) THEN
       X = 0.0
C                        NORMALISE WAVE FUNCTION TO UNITY EXACTLY!!
       DO 92 J=1,M
       DO 92 N=2,NMAX
  92   X = X + abs(Y(N,J))**2
       C = 1.0 / SQRT(X*H)
       WOLD = YINT/ABS(C)
      ELSE
C                   Normalise according to YINT
C                    (if not TRES, then YINT = KMAX - KMIN)
       C = YINT
      ENDIF
       DO 95 N=1,NMAX
       DO 93 J=1,M
93     Y(N,J) = Y(N,J) * C
95     CONTINUE
!        write(6,*) ' Normalise bin wf by ',C

C
      K = (KMIN+KMAX) * 0.5
      WNORM = 0.0
      RMS   = 0.0
      CD0   = 0.0
      CD1   = 0.0
      DO 100 J=1,M
	if(QNF(9,J)>0) Y(1,J) = 0.0
      DO 100 N=2,NMAX
      X=(N-1)*H
      WNORM = WNORM + abs(Y(N,J))**2
      RMS   = RMS   + abs(Y(N,J))**2 * X*X
         WY = 0.0
         DO 98 L=1,M
         DO 98 JF=1,NF
98       WY = WY + CC(L,J,JF,1) * FORMF(N,JF) * Y(N,L)
      CD0 = CD0 + X**(QNF(9,J)+1) * WY *H
      CD1 = CD1 + X**(QNF(9,J)+1) * Y(N,J)*H
      IF(N.EQ.MATCH.AND.J.EQ.IL) CWN = Y(N,J) * EXP(X*K) / SQFPI
100   CONTINUE
      WNORM = SQRT(WNORM*H)
      RMS   = SQRT(RMS*H)/WNORM
      CD0 =-CD0 * SQFPI      /CONV
      CD1 =-CD1 * SQFPI * K*K/CONV
      D1 = CD1
      D0 = CD0
      Q =CWN * 4*PI / CONV
      IF(WOLD.LT.0.) WOLD=WNORM
      IF(PCON.GT.0) WRITE(KO,991) WNORM,WOLD,RMS,CD0,CD1
      RMS = WNORM
      Q = D1
      IF(PCON.EQ.0)RETURN
C
  991 FORMAT(  / '  Norm =',F8.4,
     &  ' from =',F8.4,', RMS =',F8.4,' & D0 =',F12.2,' or',1p,e12.2/)
C
      IF(MOD(PCON,2).EQ.0)      RETURN
      DO 996 J=1,M
      WRITE(KO,992) J,QNF(9,J)
  992 FORMAT(/' For channel #',I3,' with L =',I2//,'   N    Y(R)')
!      if(ISC>0)  WRITE(KO,997) (N,Y(N,J),N=1,NMAXP,5)
!      if(ISC<0) WRITE(KO,998) (N,Y(N,J),N=1,NMAXP,5)
      WRITE(KO,998) (N,Y(N,J),N=1,NMAXP,5)
  996 continue
  997 FORMAT(4(I4,2X,E11.4,3X))
  998 FORMAT(3(I4,2X,2E11.4,3X))
 9997 FORMAT(F11.7,2F9.4,1x,F11.6,f9.4)
 9998 FORMAT(F11.7,1x,F9.4,1x,F11.6,4f9.3)
 9999 FORMAT(' For #',i6,i2,', K,BE =',2F11.7,': S =',2F12.2)
  999 FORMAT(' For #',i6,i2,', K,BE =',2F11.7,': S =',2F9.4,:,
     & ' & Del =',F8.3,',  D0 =',2F9.2,2g12.5,10f10.5)
      RETURN
      END
      SUBROUTINE GAUSS5(N,NR,A,M,SING,DET,EPS,SHOW)
C
C    SOLVE BY GAUSSIAN ELIMINATION SUM(J): A(I,J).P(J) = A(I,N+1)
C             WHERE P(J) IS LEFT IN MAT(J,N+1)
C
       IMPLICIT REAL*8(A-H,O-Z)
	parameter(KO=142)
       COMPLEX*16 A(NR,N+M),DET,RA
       LOGICAL SING,SHOW
       SING  = .FALSE.
      NPM = N + M
      DO 201 I=1,N
201   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
      DET = 0.
      DO 9 K = 1, N
      IF (abs(A(K,K)) .NE. 0.0 ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET/LOG(10.)
    3  FORMAT(//' THE MATRIX IS SINGULAR AT',I3,', Log10 determinant is
     &  ',2E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
C402   FORMAT( 1X,20F6.1/(1X,20F6.3))
402   FORMAT( 1X,1P,14E9.1/(1X,24E9.1))
         RETURN
    5 KP1 = K + 1
      DET = DET + LOG(A(K,K))
         RA = 1.0/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = 1
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
CDIR$ IVDEP
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = 0
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET/LOG(10.)
15    FORMAT(/' Log10 determinant is ',2F10.5)
      IF(SHOW) WRITE(KO,402) (A(I,N+1),I=1,N)
      RETURN
      END

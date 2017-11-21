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

***readcp**************************************************************
      subroutine readcp(ki,icto,icfrom,kind,ip1,ip2,ip3,ip4,ip5,
     X                  p1,p2,jmax,rmax,infile)
	real*8 jmax,p1,p2,rmax
	integer infile
	namelist/coupling/icto,icfrom,kind,ip1,ip2,ip3,ip4,ip5,
     X			  p1,p2,jmax,rmax,kfrag,kcore,nforms,infile
      	icto=0;icfrom=0;kind=0;ip1=0;ip2=0;ip3=0;ip4=-1;ip5=-1
	p1=0;p2=0;jmax=0;rmax=0;kfrag=0;kcore=0;nforms=0;infile=4
	read(KI,nml=coupling)
	if(kfrag>0) p1 = kfrag ; if(kcore>0) p2 = kcore
	return
	end
      SUBROUTINE SUMX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGT,SIGJ,SIGR,PART,
     X           EXCIT,AJUMP,FUSL,CSIG,JEX,ITC,K,RMASS,PEL,EXL,LVAL,
     X         NSA,NJA,JVAL,JPROJ,JTARG,JTOTAL,TOTFUS,FUSUM,
     X         CORFUS,SIGEL,SIGTOT,SMATS,ITCM,LXS,IF1,IF2)
	use io
	use parameters
	use searchpar, only: final
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(LXSEC=41)
      INTEGER PEL,EXL,C,ITC(MXP,MXX),EL,PART(NCH),EXCIT(NCH),LVAL(NCH),
     X        SMATS
      COMPLEX*16 SMAT(NCH)
      REAL*8 XS(MAXCH),SIGJ(0:MXPEX),SIGR(MXP,MXX),
     X       K(MXP,MXX),CSIG(LMAX1,MXPEX),SIGTOT(3),SIGEL(3)
      REAL*8 JCOEF,TOTFUS(3),CORFUS(3,NFUS1),FUSL(1+NFUS),SIGT(3)
      REAL*8 RMASS(MXP),JEX(6,MXP,MXX),JTOTAL,MSA,MJA
      REAL*8 JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH)
      COMPLEX*16 FUSUM(MAXCH,NSA*NJA),C6,S
C
      DATA  Z / 0D0 /
C
      SIGINEL = 0.0
      SIGTRAN = 0.0
      IFT = 1
      FUSL(1) = 0.0
      SMT = 2.*(1. - REAL(SMAT(EL)))
      SMEL = abs(1. - SMAT(EL))**2
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
!      write(6,*) 'C,JCOEF,S =',C,real(JCOEF),real(S)
!        exactly the same changes as frxx3.f
!
!        S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
!@@
      XS(C) = JCOEF* AMDSQS * S
      X = XS(C)*AJUMP
      IF(C.EQ.EL) THEN
         SIGT(1) = SIGT(1) + X
         SIGT(2) = SIGT(2) + X*JTOTAL
         SIGT(3) = SIGT(3) + X*JTOTAL**2
           XT = JCOEF * SMT * AJUMP
         SIGTOT(1) = SIGTOT(1) + XT
         SIGTOT(2) = SIGTOT(2) + XT*JTOTAL
         SIGTOT(3) = SIGTOT(3) + XT*JTOTAL**2
           XE = JCOEF * SMEL * AJUMP
         SIGEL(1) = SIGEL(1) + XE
         SIGEL(2) = SIGEL(2) + XE*JTOTAL
         SIGEL(3) = SIGEL(3) + XE*JTOTAL**2
         SIGJ(0) = SIGJ(0) + XS(C)
         FUSL(1) = FUSL(1) + XS(C)
        ELSE
         SIGR(IC,IA) = SIGR(IC,IA) + X
         SIGJ(IT) = SIGJ(IT) + XS(C)
         FUSL(1) = FUSL(1) - XS(C)
         IF(IC.EQ.PEL) SIGINEL = SIGINEL + XS(C)
         IF(IC.NE.PEL) SIGTRAN = SIGTRAN + XS(C)
        ENDIF
      IFT = IF2
C      ................CALCULATE SUMS FOR FUSION POLARISATIONS
        S = SMAT(C)*SQRT(S)*EXP((0.,1.)*CSIG(LVAL(EL)+1,ITC(PEL,EXL)))
        IAM = 0
        DO 685 I=1,NSA
           MSA = I-1 - JEX(1,PEL,EXL)
        DO 685 J=1,NJA
           MJA = J-1 - JEX(2,PEL,EXL)
           IAM = IAM + 1
        C6 = S  * CLEB6(JVAL(EL),-MSA,JPROJ(EL),MSA,LVAL(EL)+Z,Z)
     X          * CLEB6(JTOTAL,-(MSA+MJA),JTARG(EL),MJA,JVAL(EL),-MSA)
         FUSUM(C,IAM) = FUSUM(C,IAM) + C6
  685    continue
  690 CONTINUE
C
      X = FUSL(1)*AJUMP
      TOTFUS(1) = TOTFUS(1) + X
      TOTFUS(2) = TOTFUS(2) + X*JTOTAL
      TOTFUS(3) = TOTFUS(3) + X*JTOTAL**2
      DO IT=1,NFUS
      X = FUSL(1+IT)*AJUMP
      CORFUS(1,IT) = CORFUS(1,IT) + X
      CORFUS(2,IT) = CORFUS(2,IT) + X * JTOTAL
      CORFUS(3,IT) = CORFUS(3,IT) + X * JTOTAL**2
      ENDDO
      IF(LXS.LT.0) WRITE(LXSEC,100) JTOTAL,LVAL(EL),SIGJ(0),SIGINEL,
     X                              SIGTRAN,FUSL(1)
      IF(LXS.LT.0) written(LXSEC) = .true.
 100  FORMAT(F6.1,I6,1P,4E10.3,' (SIG FOR J/L)')
 1450 FORMAT(' Reaction Xsec',F8.1,A1,'/',I2,' @',I2,' =',F8.3,A1,',',
     X   ' Out:', 9(F8.3,A1),:,/(11X,'Xsec',23X,10(F8.3,A1)))
C
      IF(SMATS==2.and.final)
     X  WRITE(7,719) SMAT(EL),LVAL(EL),JVAL(EL),JTOTAL
      IF(SMATS==3.and.final) then
	do C=1,NCH
	if(PART(C)==PEL.and.EXCIT(C)==EXL) 
     X    WRITE(7,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
     X    LVAL(EL),JVAL(EL)
	enddo
	endif
      IF(SMATS.GE.4.and.final) THEN
      DO 720 C=1,NCH
        IF(ABS(SMAT(C)).GT.0E-10)
     X WRITE(7,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
     X    LVAL(EL),JVAL(EL)
     *  ,jproj(c),jtarg(c),jproj(el),jtarg(el)
     *  ,C,el   ! added by neil summers  Nov 2009
     *  ,K(PART(C),EXCIT(C))/RMASS(PART(C))/(K(PEL,EXL)/RMASS(PEL))

!	write(79,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!     X    LVAL(EL),JVAL(EL) ,jproj(c),jtarg(c),jproj(el),jtarg(el)
!	call flush(79)
!	write(80,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!     X    LVAL(EL),JVAL(EL) ,jproj(c),jtarg(c),jproj(el),jtarg(el)
!      IC = PART(C)
!      IA = EXCIT(C)
!      if(JPROJ(C)/= JEX(1,IC,IA).or.JTARG(C)/= JEX(2,IC,IA)) then
!	write(78,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!     X    LVAL(EL),JVAL(EL) ,jproj(c),jtarg(c)
!      endif
720    CONTINUE
        ENDIF
	if(SMATS.ge.2.and.final) written(7) = .true.
!	if(SMATS.ge.2) call flush(7)
 719  FORMAT(2F15.10,I6,2F6.1,' : S(L,J,JT)')
C750  FORMAT(2F15.10,I6,2F6.1,I6,F6.1,' S(L,J,JT,L-IN,J-IN)')
C750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,' S-MAT')
C750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,4F4.1)
 750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,4F5.1
     x       ,2I4,f12.8        ! Neil Summers
     x      )
C		Print out cross section information in file 13
      RETURN
      END
	
****READIN**************************************************************
	SUBROUTINE READIN(WOUT,WOUT0,NF0,NF,STREN,ISNONO,
     X     NAME,MASS,NEX,QVAL,RMASS,HP,COPY,EXTRA,IOFAM,GIVEXS,
     X     JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X     PWFLAG,FORMF,FORML,FORMC,PTYPE,CPSO,EPS,MMXQRN,FILE,ffreal,
     X  HPOT,LAMBDA,MATRIX,MEK,NIX ,NSP, QNF,D0,BE,AFRAC,BPHASE,
     x  ICTO,ICFROM,KIND,QQ,IREM,KPCORE,BETAR,BETAI,JMAX,RMAX,REV,
     x  NFI,TWOW,OPN,NPRIOR,NPOST,REW14,LOCAL,NIB,LTRANS,FORMDEF,
     x  MLCALL,BPROJ,XA,XB,XP,XQ,SPINTR,FPT,NKP,VARYL,LSHAPE,LDEP,
     x  NLL,LOCF,KLT,COUPLE,ICOM,ICOR,NCP,IEXCH,GPT,POTCAP,INFILE,
     x  NSA, NJA, XCOEF, DISCIT, DISC8, LAMAX,NBINS,EMID,KMINX,NKBIN,
     x  NORBIN,STOP19,WID,CENTR)
	use parameters
	use factorials
	use drier
	use kcom
	use trace
	use io
	use fresco1
	use searchdata, only: energy_list,num_energies
	use searchpar, only: final
	IMPLICIT NONE
C
C    FORM FACTORS AND THEIR PARAMETERS
C    ---------------------------------
      COMPLEX*16 FORMF(MAXM,MLOC),FORMC(MAXNLC,MSP,2)
      REAL*8 FORML(MAXNLR,MSP,2),VIMX,VARYL(2,MLOC),EMID(MSP),
     X		KMINX(2,MSP),VRMX,FORMDEF(7,MLOC),DELTAE(MSP)
C
      INTEGER NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     X        ITC(MXP,MXX),CPOT(MXP,MXX),PTYPE(8,MLOC),NORBIN(MSP)
      INTEGER QNF(19,MSP),BHEAD(2,2,0:9),LSHAPE(MLOC),NKBIN(2,MSP)
      INTEGER CP,ICTO(MAXCPL+1),ICFROM(MAXCPL+1),KIND(MAXCPL+1),
     X        LOCF(MAXCPL),FPT(6,MAXQRN,MAXCPL),NKP(2,MAXCPL),
     X        IREM(MAXCPL+1),KPCORE(MAXCPL+1),ICOR(MAXCPL,2),
     X        GPT(2,MAXQRN,MAXCPL),FILE(MAXCPL),NFI(3),
     X        NLL(MAXCPL),KLT(MAXCPL),NIB(MAXCPL),
     X        QQ(MAXCPL+1),ICOM(MAXCPL,2),NBINS,
     X        QCM(MAXCPL+1),LAM(MAXCPL+1),INFILE(MAXCPL+1)
C
C    CONTROL VARIABLES
C    -----------------
      LOGICAL GIVEXS(MXP),EXTRA(2,MXP,MXX),FLIP
      LOGICAL LOCAL(MAXCPL),REV(MAXCPL),USED(2,MSP),COUPLE(MAXCPL),
     X        NPRIOR,NPOST,OPN(2),ffreal(MAXCPL),REW14,WOUT,WOUT0
      LOGICAL LTRANS(MAXQRN,MAXCPL),MLCALL(2,MAXCPL),HASO(KPMAX)
      LOGICAL DISC8,DISCIT,PWFLAG(MXP),CPSO(MAXCPL),STOP19
	integer NF,NF0,LAMAX,KN,NSP,NCP,IEXCH,NSA,NJA,LDEP(MLOC)
        real*8 XCOEF,JAP,R
	integer NIX,NCHAN,IC,ITCM,L,JF,J,IN,I,IITER,IOFAM(2,MXP,MXX)
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      REAL*8 MASS(4,MXP+1),RMASS(MXP),HP(MXP),JEX(6,MXP,MXX)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1),CONV
      REAL*8 BE(MSP,5),D0(MSP,2),BPROJ(MSP,MXP),HPOT(MLOC)
      INTEGER NCHPMAX(MXP)
C
C    COUPLINGS
C    ---------
      REAL*8 BETAR(MAXCPL+1),BETAI(MAXCPL+1),JMAX(MAXCPL+1),
     X	RMAX(MAXCPL+1),   BPHASE(max(1,NKMAX),MSP),
     x  WID(MAXCPL+1),CENTR(MAXCPL+1),
     X  XA(MAXCPL),XB(MAXCPL),XP(MAXCPL),XQ(MAXCPL),STREN(MLOC)
      REAL*8 AFRAC(MXPEX,MXPEX,2,MSP),MEK(MPAIR),SPINTR(2,MPAIR)
      INTEGER MATRIX(6,MPAIR),LAMBDA(MLOC),MMXQRN,POTCAP(MLOC,2,MAXCPL)
	LOGICAL TWOW,ISNONO
C
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
      real*8 EPS
C
      CHARACTER*8 NAME(2,MXP+1)
      CHARACTER*3 WAYS(2),W
      CHARACTER*1 PSIGN(3)
      DATA WAYS/'ONE','TWO'/
C
      CALL CHECK(N,MAXN,7)
      CALL CHECK(MINT,MAXM,7)
      CALL CHECK(NLN,MAXNLN,8)
      CALL CHECK(INT(JTMAX),LMAX,11)

C
	if(melfil.eq.1) then
            open (53,access = 'sequential',
     x         form='unformatted',status = 'unknown')
            open (54,access = 'sequential',
     x         form='formatted',status = 'unknown')
	endif

      WOUT0 = VEFF.NE.0 .OR. MOD(WAVES,2).NE.0 .OR. WDISK.NE.0
      WOUT = WOUT0 .OR. NFUS.NE.0 

	if(final) then
	if(SMATS>=2) then
	   open(45,form='formatted',access='sequential')
	   rewind 45
	   endif
	if(SMATS>=2) then
	   open(38,form='formatted',access='sequential')
	   rewind 38
	   endif
	open(39,form='formatted',access='sequential'); rewind 39
	open(56,form='formatted',access='sequential'); rewind 56
	open(40,form='formatted',access='sequential',recl=210);rewind 40
	open(13,form='formatted',access='sequential'); rewind 13
	endif
	HASO(:) = .false.

C	NOW READ PRE-DIGESTED INPUT IN FILE 3:
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      CALL PARTEX(NAME,MASS,NEX,QVAL,RMASS,HCM,HP,COPY,EXTRA,
     X     GIVEXS,JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X     NCHPMAX,PWFLAG,FCWFN,IOFAM)
C
C    ENTER THE VARIOUS TYPES OF POTENTIALS
C    -------------------------------------
      NF0 = 0
	STREN(:)=0d0
   

      CALL POTENT(FORMF,NF0,PTYPE,MAXM,HCM,TRENEG,MR,M,HASO,
     X            HP,HPOT,STREN,LAMBDA,NCHPMAX,VARYL,LDEP,LSHAPE,
     X            MATRIX,MEK,NIX, JEX,MASS,CPOT,NEX,NCHAN,COPY,STOP19,
     x 		  FORMDEF)

!  debugging only!:
!	STREN(:)=0d0
!	print *,' ALL ASYMPTOTIC STRENGTHS == ZERO!!!!'
C
C
C    ONE AND TWO - PARTICLE FORM FACTORS AND THEIR PARAMETERS
C    --------------------------------------------------------
C
      CALL INFORM(FORML,FORMC,NSP,FORMF,NF0,PTYPE,TRENEG,QVAL,BPHASE,
     X  QNF,D0,BE,NCHAN,NEX,ENEX,JEX,BAND,MASS,AFRAC,ITC,COPY,N,
     X  NLN,MINT,NNN,RNN,RMIN,HP(1),RIN,MR,NNU,ERANGE,DK,PI,NAME,
     X	NBINS,EMID,DELTAE,KMINX,NKBIN,NORBIN)
	do i=43,44
	if(written(i)) close(i)
	enddo
C
C
C    COUPLINGS
C    ---------
      CP = 1
      NF = NF0
      NFI(1) = 1
      NFI(3) = 1
      TWOW = .false.
      ISNONO = .false.
      LREC66 = 0
      OPN(1) = .FALSE.
      OPN(2) = .FALSE.
      NPRIOR =.FALSE.
      NPOST = .FALSE.
      REW14 = .FALSE.
      BPROJ(:,:) = 0.
      WID(:) = 0.; CENTR(:) = 0.
	IITER=ITER
	if(IBLOCK<0) IITER=1
	MCALLS = .false.
	ALLPART = .false.
   90 continue
!      READ(KI,1220) ICTO(CP),ICFROM(CP),KIND(CP),QQ(CP),IREM(CP),
!     X         KPCORE(CP),BETAR(CP),BETAI(CP),JMAX(CP),RMAX(CP)
      CALL readcp(KI,ICTO(CP),ICFROM(CP),KIND(CP),
     X         QQ(CP),IREM(CP),KPCORE(CP),QCM(CP),LAM(CP),
     X         BETAR(CP),BETAI(CP),JMAX(CP),RMAX(CP),INFILE(CP))
!1220 FORMAT(3I4,3I2,2F8.2,2F4.1)
  100 IF(ICTO(CP).EQ.0.OR.ICFROM(CP).EQ.0) GO TO 120
      if(KIND(CP)>8) KIND(CP) = KIND(CP)-8   ! temporary compatibility
      CALL CHECK(CP,MAXCPL,4)
      LOCAL(CP) = KIND(CP).LE.6.and..not.(KIND(CP)==1.and.QQ(CP)==1)
      REV(CP)   = ICTO(CP).GT.0
	TWOW = TWOW .or. (REV(CP).and.KIND(CP).ge.5)
	ISNONO = ISNONO .or. KIND(CP).eq.8
      W = WAYS(1)
      IF(REV(CP)) W = WAYS(2)
      IF(JMAX(CP).lt..01) JMAX(CP) = JTMAX + 0.51
      IF(RMAX(CP).lt..01) RMAX(CP) = abs(RMATCH) - 3*HCM
      NLL(CP) = MIN(NLN, NINT(RMAX(CP)/RINTP)+1 )
      LOCF(CP) = 0
      KLT(CP) = 1
      WRITE(KO,1230) W,CP,ICTO(CP),ICFROM(CP),KIND(CP),
     X           QQ(CP),IREM(CP),KPCORE(CP),QCM(CP),LAM(CP),
     X           BETAR(CP),BETAI(CP),JMAX(CP),
!    x           RMAX(CP)
     x           (NLL(CP)-1)*RINTP
 1230 FORMAT(/' ',132('*')//' ',A3,'-way',
     X        ' COUPLING #',I2,' for partitions ',I2,' <- ',I2,
     X' of KIND',I3,', '               ,5I3,
     X' & P1,P2 =',F10.4,F9.4,' : for J <=',F8.1,' & R <',
     X F6.1,' fm.'/)
       CALL FLUSH(KO)
      ICTO(CP)  = ABS(ICTO(CP))
      IC = ICTO(CP)
      COUPLE(CP) = IC.LE.NCHAN .AND. ICFROM(CP).LE.NCHAN
	ALLPART = ALLPART .or. 
     x      (LOCAL(CP).and.ICTO(CP).ne.ICFROM(CP))
!     x      (LOCAL(CP).and.REV(CP).and.ICTO(CP).ne.ICFROM(CP))
      DO 110 IN=1,2
        ICOM(CP,IN) = IC
        ICOR(CP,IN) = ICFROM(CP)
        IF(MASS(IN,ICOR(CP,IN)).LT.MASS(IN,ICOM(CP,IN))) GO TO 110
        ICOM(CP,IN) = ICFROM(CP)
        ICOR(CP,IN) = IC
C       NOW, ICOR(CP,IN) = CORE PARTITION & ICOM(CP,IN) = COMPOSITE.
  110   CONTINUE
C
C    DEAL WITH PARTICULAR KIND OF COUPLING
C    .....................................
      CALL INTER(ICTO(CP),ICFROM(CP),IC,KIND(CP),
     X           QQ(CP),IREM(CP),KPCORE(CP),QCM(CP),LAM(CP),
     X           BETAR(CP),BETAI(CP), KLT(CP),NIB(CP),HPOT,
     X        REV(CP),NLL(CP),LOCF(CP),COUPLE(CP),ICOR,ICOM,CP,IITER,
     X        XA(CP),XB(CP),XP(CP),XQ(CP),RSP,ffreal(CP),FILE(CP),
     X        MASS,HP,JEX,ENEX,QVAL,NEX,BAND,COPY,ITC,CPOT,
     X        PTYPE,NAME,FORMF,FORML,FORMC,BE,D0,BPROJ,POTCAP(1,1,CP),
     X        QNF,AFRAC,FPT(1,1,CP),NKP(1,CP),GPT(1,1,CP),
     X        NFI,USED,NPRIOR,NPOST,OPN,M,HCM,EPS,MMXQRN,
     X        NF,NF0,NSP,NLCN,RIN,WID(CP),CENTR(CP),INFILE(CP),
     X        MATRIX,MEK,NIX,SPINTR, STREN,LAMBDA,HASO,CPSO(CP),
     X        LTRANS(1,CP),REW14,MTMIN,MLCALL(1,CP),MLCALL(2,CP),MCALLS,
     X        NBINS,PSIGN,EMID,DELTAE)
C
      CP = CP + 1
      IF(KIND(CP-1).NE.7.OR.ABS(IREM(CP-1)).NE.2) GO TO 90
C          CONSTRUCT A NON-ORTHOGONALITY SUPPLEMENT FOR ABOVE COUPLING:
      ICTO(CP) = ICTO(CP-1)
          IF(.NOT.REV(CP-1)) ICTO(CP) = - ICTO(CP)
      ICFROM(CP) = ICFROM(CP-1)
      KIND(CP) = 8
      QQ(CP) = QQ(CP-1)
      IREM(CP) = 0
      KPCORE(CP) = 0
      BETAR(CP) = 0.0
      RMAX(CP) = RMAX(CP-1)
      JMAX(CP) = JMAX(CP-1)
      GO TO 100
C
  120 NCP= CP -1
      DO 121 JF=NF0+1,NF
  121 PTYPE(6,JF) = 0
      DO 125 KN=1,NSP
       L = QNF(9,KN)
       IF(L.EQ.0) GO TO 125
       DO 1245 J=1,2
       if(.not.cxwf) then
       DO 124 I=1,NLN
        R = (I-1)*RINTP
  124   FORML(I,KN,J) = FORML(I,KN,J) * R**L
	else
       DO 1241 I=1,NLN
        R = (I-1)*RINTP
 1241   FORMC(I,KN,J) = FORMC(I,KN,J) * R**L
	endif
 1245  CONTINUE
  125  CONTINUE
      IF(LREC66.GT.0) open(66,access='sequential',status='scratch',
     X                    form='unformatted')
      IF(LREC66.GT.0) REWIND 66

C
C    END OF COUPLING TYPES.    NOW FIND INCOMING BEAM
C    ................................................
 
!1240 FORMAT('0File ',I3,' needs NR =',I6)
         IEXCH = 0
       IF(COPY(2,PEL,EXL,2).EQ.PEL) THEN
!         IEXCH = 1
!         IF(mod(nint(2*JEX(1,PEL,EXL)),2)==1) IEXCH = -1   ! Fermion-like for 1/2-int spins
	 IEXCH = (-1)**nint(2*JEX(1,PEL,EXL))   	! symmetry in elastic channel
       ENDIF
      if(final) WRITE(KO,1250) PEL,EXL,LAB,LIN,LEX
 1250 FORMAT(/' Incoming partition',I3,' in excitation state #',I2,
     X ',   Laboratory Energy given for partition',I3,' Nucleus',I2,
     X ' in Excitation pair',I2/)
      IF(IEXCH.NE.0) WRITE(KO,*) 'Scattering of identical particles with
     X SYMMETRISATION FACTOR = ',IEXCH
      DO 140 J=1,NJ+1
         JAP = JBORD(J) + JEX(1,PEL,EXL) + JEX(2,PEL,EXL)
  140    IF(JAP.GT.AINT(JAP+.1)+.1) JBORD(J) = JBORD(J) + 0.5
C
      NSA = 2*JEX(1,PEL,EXL) + 1.1
      NJA = 2*JEX(2,PEL,EXL) + 1.1
      XCOEF = 10.0/(NSA * NJA)
      DISCIT= NPOST .OR. NPRIOR .or. INITWF/=0
      DISC8 = DISCIT
	LAMAX = maxval(LAMBDA(1:MLOC))
	LAMAX = max(LAMAX,1)
C
C    INCOMING ENERGIES
C    -----------------
!      READ(KI,1255) (ELAB(I),NLAB(I),I=1,3),ELAB(4)
!      READ(KI,nml=incident)
!1255 FORMAT(3(F8.4,I8),F8.4)
      if(num_energies==0) then
      IF(NLAB(1)+NLAB(2)+NLAB(3).GT.-1)
     X WRITE(KO,1256) (ELAB(I),ELAB(I+1),NLAB(I),I=1,3)
 1256 FORMAT('0Lab. ENERGY ranges :',/,
     X       3(/'   from ',G14.6,' to',G14.6,' in',i6,' intervals'))
      else
       if(final) WRITE(KO,1257) (energy_list(I),I=1,num_energies)
 1257 FORMAT('0Lab. ENERGY ranges :',/(1x,10f8.3))
      endif
        ELAB(5) = ELAB(4)
        NLAB(4) = 1

 	I = RMATR/HCM+1 - 10
	I = min(I,MINT-5)
	VRMX = 0d0
	VIMX = 0d0
	do JF=1,NF
        VRMX = max(VRMX,abs(real(FORMF(I,JF))))
        VIMX = max(VIMX,abs(aimag(FORMF(I,JF))))
        enddo
	if(final) write(ko,1280) (I-1)*HCM,VRMX,VIMX
1280	format(/' Largest real,imaginary parts of any form factor at R='
     X   ,f8.2,' are ',1p,2E10.2,' MeV')

        RETURN
	END SUBROUTINE READIN

	subroutine fkind(written,KO)
	character*40 flkind(299)
	integer writf(299),nwrit
	logical written(299)
	flkind(:) = ' '
	flkind( 1) = 'two-particle multipoles'
	flkind( 2) = 'KIND=9 nonlocal formfactor'
	flkind( 3) = 'local copy of User input'
	flkind( 4) = 'external form factors & potentials'
	flkind( 5) = 'standard input'
	flkind( 6) = 'standard output'
	flkind( 7) = 'elastic  S-matrix elements'
	flkind( 8) = 'single-particle & scattering wfs'
	flkind( 9) = 'complex transfer multipoles'
	flkind(10) = 'S-matrix elements'
	flkind(11) = 'real transfer multipole'
	flkind(12) = 'transfer kernels'
	flkind(13) = 'total cross sections/state'
	flkind(14) = 'interaction potentials'
	flkind(15) = 'local equivalent potentials'
	flkind(16) = 'tables of cross sections'
	flkind(17) = 'output scattering waves'
	flkind(18) = 'wfns of ''best'' iterate'
!	flkind(19) = 'temporary output amplitudes'
	flkind(20:33) = 'user files'
	flkind(34) = 'input optical potentials'
 	flkind(35) = 'Astrophysics S-factors  / Ecm'
	flkind(36) = 'output scattering AMPL amplitudes'
	flkind(37) = 'output scattering FAM amplitudes'
	flkind(38) = 'cross sections for each J/pi'
	flkind(39) = 'cross sections for each Ecm'
	flkind(40) = 'all cross sectns. for each Elab'
	flkind(41) = 'source terms at each iteration'
	flkind(42) = 'bin wavefunctions for each E'
	flkind(43) = 'bin phase shifts as k functions'
	flkind(44) = 'bin phase shifts as E functions'
	flkind(45) = 'scat phase shift as E functions'
	flkind(46) = 'bs wave functions & ANC ratios'
	flkind(47) = 'reduced matrix elements'
	flkind(48) = 'concurrency log file'
	flkind(49) = 'S-matrix elements concurrent->10'
	flkind(55) = 'Single-particle wfns'
	flkind(56) = 'Fusion for each Jtotal'
	flkind(57) = 'CDCC amplitudes'
	flkind(58) = 'Single-particle wfns'
	flkind(59) = 'Single-particle vertex functions'
	flkind(70) = 'Fusion for each J/pi'
	flkind(75) = 'S-factors for lab energies'
        flkind(90) = 'average polarisation potential'
        flkind(91) = 'polarisation potl for rereading'
	flkind(99) = 'input J/L rescaling factors'
	flkind(201:210) = 'Separate cross sections'
	nwrit = 0
	do i=1,299
	 if(written(i)) then
		flkind(i) = trim(flkind(i))//'.'
	   nwrit = nwrit+1
	   writf(nwrit) = i
	 endif
	enddo
	write(KO,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990	format(/'  The following files have been created:',
     X	/(1x,2(i3,':',a35)))
	return
	end

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

***frxx1.f*********************************************************
	MODULE arrays

C    FORM FACTORS AND THEIR PARAMETERS
C    ---------------------------------
      COMPLEX*16,allocatable:: FORMF(:,:),ELCOEF(:),FORMC(:,:,:)
      REAL*8,allocatable:: FORML(:,:,:),FORMDEF(:,:)
      REAL*8,allocatable:: AFRAC(:,:,:,:),OO(:,:),EVEC(:,:)
C
C    PARTIAL WAVES AND THEIR PARAMETERS
C    ----------------------------------
      COMPLEX*16,allocatable:: CLIST(:,:,:),SRC(:,:)
      COMPLEX*16,allocatable:: PSI(:,:),WNM(:,:,:),FNC(:,:)
      integer,allocatable:: NFLIST(:,:,:),NCLIST(:,:)
      integer,allocatable:: levelp(:),levelt(:)
      CHARACTER*1,allocatable:: INCOME(:)
C
C    COULOMB FUNCTIONS
C    -----------------
      COMPLEX*16,allocatable:: CH(:,:)
      REAL*8,allocatable:: CFMAT(:,:,:),CGMAT(:,:,:),CRCRAT(:)
C
C    SOLVING THE COUPLED EQUATIONS
C    -----------------------------
      COMPLEX*16,allocatable:: SMAT(:),SOMA(:,:),SPAD(:,:,:),SRATIO(:),
     X		FIM(:,:),FIMD(:,:)
      complex*16,allocatable:: ferw(:)
      INTEGER,allocatable:: LVAL(:),PART(:,:),EXCIT(:,:),CUTVAL(:),
     X        INITL(:),CHNO(:,:)
      LOGICAL,allocatable:: finishedcc(:),READWF(:)
      LOGICAL,allocatable:: BLOCKD(:),SAME(:),EMPTY(:),
     X			    SKIP(:,:),SIMPLE(:)
      REAL*8,allocatable:: BPHASE(:,:),EXCH(:,:),alphas(:,:,:)

      REAL*8,allocatable:: XS(:),SIGFUS(:,:)
      COMPLEX*16,allocatable:: FUSUM(:,:)
      REAL*8,allocatable:: JVAL(:),JPROJ(:),JTARG(:),ECM(:,:),LL1(:)
	End MODULE arrays

****fr*********************************************************
      subroutine fr
	use parameters
	use factorials
	use drier
	use kcom
	use trace
	use io
	use fresco1
	use gails, ipmax=>numax, cff=>cf
	use searchdata
	use searchpar
	use arrays
!      implicit real*8(a-h,o-z)
      implicit none
C
C    PARTIAL WAVES AND THEIR PARAMETERS
C    ----------------------------------
      COMPLEX*16 PHI(MAXN,2),VOPT(MAXNLN),VPOT(MAXN,3)
      integer NCHPART(MXP),CHBASE(MXP),nphases,linel(100)
      real*8 phase(20),phasin(100),phaslst(20),phaseadd(20)
      COMPLEX*16 CHL(LMAX1,MXPEX,2),SMATI
      LOGICAL DERIV,INITWFE
      REAL*8 CF,CG,R,RH, JTOTALI,JVALI,JVALE,HPI,ENLABI
      integer NAI,PARITYI,LVALI,LVALE,ITI,NRF


C
C    SOLVING THE COUPLED EQUATIONS
C    -----------------------------
      COMPLEX*16 C6,C7,S,CI
      INTEGER NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     X        ITC(MXP,MXX),CPOT(MXP,MXX),PTYPE(8,MLOC),M1
      INTEGER QNF(19,MSP),BHEAD(2,2,0:9), NSTEPD
      INTEGER CP,ICTO(MAXCPL+1),ICFROM(MAXCPL+1),KIND(MAXCPL+1),
     X        LOCF(MAXCPL),FPT(6,MAXQRN,MAXCPL),NKP(2,MAXCPL),
     X        IREM(MAXCPL+1),KPCORE(MAXCPL+1),ICOR(MAXCPL,2),
     X        GPT(2,MAXQRN,MAXCPL),FILE(MAXCPL),NFI(3),
     X        NLL(MAXCPL),KLT(MAXCPL),NIB(MAXCPL),
     X        QQ(MAXCPL+1),ICOM(MAXCPL,2),PARITY,EL,C,C1,C2,
     X        CUTOFF,APFLG(5),NLN1,ip,id
	real*8 AERR,GMAX,DEGENY,APEPS,WERRMAX,WERR,DIFF,FRACT
	real*8 WID(MAXCPL+1),CENTR(MAXCPL+1)
	logical prwf
C
C    CONTROL VARIABLES
C    -----------------
      INTEGER DONE,DONES,IOFAM(2,MXP,MXX),KOO
      LOGICAL GIVEXS(MXP),EXTRA(2,MXP,MXX),FLIP,LCERW1,LCERW2
      LOGICAL LOCAL(MAXCPL),REV(MAXCPL),WREQ,COUPLE(MAXCPL),FALLOC,
     X        NPRIOR,NPOST,OPN(2),ffreal(MAXCPL),REW14,OPEN12,WOUT,WOUT0
      LOGICAL LTRANS(MAXQRN,MAXCPL),MLCALL(2,MAXCPL),CPSO(MAXCPL),
     X        SCHONP,REPEAT,OPEN18,STOP19
      LOGICAL FAIL,DISC8,DISCIT,PWFLAG(MXP),FJSWTCH,SMALLJ,HERE,SHFL,OP
	integer NF,NF0,LAMAX,NC,NCLREQ
	integer NCP,IEXCH,NSA,NJA,NSB,NJB
        real*8 TCC,TCC0,JUMPER
        real*8 TIME0,R0,XCOEF,TIME0J,RASYMAX,RTURN1,EE,AL,BEST,RERR,R1
	real*8 T,TH1,TH2,ALOSSM,TSTD,RENORM,RTURN,RS,GAP,AMSA,XLMAX,X
	integer JPSET,LUSED,JCCSET,NCH,IEX,MINTL,MEL,NTHRESH,MCCSET
	integer NIX,IAME,IEN,ILEN,NCHAN,MAXF,DROPPED,LEN,NLEN
	integer IC,IA,NA,MAM,IF1,IF2,ITCM,IMA,JBLOCK,L,IPARIT,KS
	integer JF,JIN,IAM,J,NLJ,IN,IPUT,IGET,LAP,IF,I,IT
	integer KN,NSP,L1,IE,NPROJ,NTARG,KINTL,JUMP2
	integer NICH,MAXXCH,IB,JFT,JFTT,NFL,NR,IR0,IR,IC1,KKP
	integer IC2,IN1,ICH,KP,KPB,IX,KK,II,NBEST,NSOLI,ITNL,IPARI
	integer NJDONE,LL,JJ,I0,nbas,CHANSI,JF0,NFD,NEXK
	integer LAST,J2LAST,J2,KNLAST
	integer I1,I2,IX1,IX2,LAM1,LAM2
	integer NJTOTAL,IJTOTAL
! AMM --------------------------------------
      real*8 rturnc
      real*8 x1,ecmi,ecmf,triangle,mp,mt,ep,et,erel,mbig,msmall
!-------------------------------------------

C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      PARAMETER(NTHRESH=5)
      REAL*8 MASS(4,MXP+1),RMASS(MXP),HP(MXP),JEX(6,MXP,MXX)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1),T4,RESM(NTHRESH),THRJ(NTHRESH)
      REAL*8 BE(MSP,5),D0(MSP,2),BPROJ(MSP,MXP),HPOT(MLOC),BEPROJ
      REAL*8 ENEX0(2,MXP,MXX)
C
C    BPM and VEFFPOT
C    ---------------
      REAL*8 BARE,TLE,SIGL,HOM,EKL,EKL1,DER,RSO,POS,BAR,RTURN2
      REAL*8 FUSIO,RSOO
C
C    COUPLINGS
C    ---------
      REAL*8 BETAR(MAXCPL+1),BETAI(MAXCPL+1),JMAX(MAXCPL+1),SCALEL,
     X	     RMAX(MAXCPL+1),MEK(MPAIR),SPINTR(2,MPAIR),VARYL(2,MLOC),
     X       XA(MAXCPL),XB(MAXCPL),XP(MAXCPL),XQ(MAXCPL),STREN(MLOC),
     X       EMID(MSP),KMINX(2,MSP),SPINV,MASSV,CHRGV,TENSOR,TENSOR2
      INTEGER MATRIX(6,MPAIR),LAMBDA(MLOC),LSHAPE(MLOC),LDEP(MLOC),
     X   NKBIN(2,MSP),IBIN(MXX),NBINS,NNJMAX,PARV,NORBIN(MSP),
     X   CHSIGN(MXX),BSIGN(0:MSP),NFUSCH,mag,POTCAP(MLOC,2,MAXCPL),
     x   INFILE(MAXCPL)
	LOGICAL TWOW,ISNONO
C
C    INCOMING ENERGIES
C    -----------------
      REAL*8 K(MXP,MXX),ETA(MXP,MXX),ETOTAL,ECMC(MXP,MXX),ENLAB,EOFF,
     X     ENCOM,RMK,RMKD,CFG(MAXMUL,4),DE,ELAST,CSIG(LMAX1,MXPEX)
      REAL*8 GAM(MXP,MXX)  ! AMoro
      REAL*8 JTOTAL,JAP,JSWITCH,JAL,JN,LJMAX,SSWITCH,JNLAST
C
C    ARRAYS FOR NON-LOCAL COUPLING FORM FACTORS
C    ------------------------------------------
      REAL*8 DNL(MAXNLO)
      REAL*8 NLOC(NLM,MAXCPL)
      COMPLEX*16 FFC
      EXTERNAL FFC
C
C    CALCULATION OF CROSS SECTIONS
C    -----------------------------
      REAL*8 JCOEF,TOTFUS(3),CORFUS(3,NFUS1),SIGT(3),FUSL(1+NFUS1),
     x       SIGTOT(3),SIGEL(3),datanorm,DSPINS(MXX)
      REAL*8 SIGJ(0:MXPEX),SIGR(MXP,MXX),OUTJ,SFAC(MXPEX),rtrace,
     X FUSJPI,CHSIZES(MXPEX),XSC,AMDSQS,FUSJ,TOTJ,CFUSJ(NFUS1),sfa
      integer LCROSS,LFAM,LXSEC,NANGL,SMALLS(MXPEX),NSMALL,
     x        CHPRES(MXPEX),LEG,MAXPLM,DMULTIES(MXX),NMULTIES,LGAM,
     x        DLEVEL(MXX)
C
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
      COMMON /TIMER/ START
      real*8 START,TIME1,SECOND,memuse
      real*4 TM,TM0,TSYNC
      real*8 Z,HALF,EPS,DDEXP
C
      CHARACTER*80 CHME*3,LINE
      CHARACTER*8 NAME(2,MXP+1),NAMEV,DNAME
      CHARACTER*4 W
      CHARACTER*1 SCALE(0:MXPEX),PSIGN(3)
C
      integer NFNL,MMXF,MMXCH,MMXQRN,MMXIT,MAL1,MMXICH,MAXCHT
      DATA NFNL,MMXF,MMXCH,MMXQRN,MMXIT,MAL1,MMXICH/7*0/
      DATA PSIGN / '-','?','+' /
      DATA RESM / 0.01, 0.1, 0.5, 0.90, 0.99 /
      parameter(mag=1)
   
      DDEXP(T) = EXP(MIN(100D0,MAX(-100D0,T)))
      TM(I) = real(SECOND() - TIME0)
C
      START = 0.0
      X=0d0
      TIME0 = SECOND()
      I=0; T = dble(TM(I))
!      write(KOI,*) 'Starting  @ ',real(T)
	koo = 93
        TSYNC=0.0
C
      Z = 0.0
      HALF = 0.5D0
      CI = (0D0,1D0)
      EPS = 1E-10
!	scalx = 1d-100; rscalx = 1d0/scalx

      NFNL=0 ; MMXF=0 ; MMXCH=0 ; MMXQRN=0 ; MMXIT=0 ; MAL1=0
      NCLREQ=0; LUSED = 0; MMXICH=0 ; MAXXCH=0
      NIX=1
	CHSIGN(:) = 0
	numafrac=0
	rewind 3
	call rewop(4)
	if(KO/=6) rewind KO  ! keep only final Fresco run!
	if(rterms) call rewop(60)
	if(rterms.and.pralpha) call rewop(61)
	NANGL = (abs(THMAX) - THMIN)/abs(THINC) + 1 + 0.5	
      data_chisq(:) = 0.0
      STOP19 = .false.
C
      if(final) then
	   inquire(7,opened=op); if(op) close(7)
      OPEN(7,form='formatted',file='fort.7',recl=110)
!      OPEN(80,form='formatted',file='fort.80',recl=90)
	   inquire(16,opened=op); if(op) close(16)
      OPEN(16,form='formatted',file='fort.16')
      call rewop(16)
      call rewop(35)
      call rewop(7)
      endif
C
************************************************************************
      IAME = 0
       I = ichar('0')
      CHME = char(I+mod(IAME/100,10))//char(I+mod(IAME/10,10))//
     X       char(I+mod(IAME,10))
      TMPN = trim(TMP)//'fort.'//CHME//'.'
!	rewind 48
      if(TRENEG>0) then
	open(34,form='formatted',access='sequential')
	rewind 34
	endif
      OPEN12 = .false.
      OPEN18 = .false.
      XA=0.;XB=0.;XP=0.;XQ=0.;MATRIX=0;MEK=0.
	ENEX(:,:,:) = 0.; BE(:,:) = 0.

      if(.not.allocated(AFRAC)) then

      memuse = (MAXM/1024.)*(MLOC/1024.)*(16/1024.)
      if(memuse>0.1) write(KO,'(a,i5,i8,a,f7.3,a)')'Allocating FORMF(',
     &MAXM,MLOC,') in ',memuse,' GB' 
      call flush(KO)

      allocate(FORMF(MAXM,MLOC),ELCOEF(MLOC),FORMDEF(7,MLOC))
      if(.not.cxwf) then
	allocate(FORML(MAXNLN,MSP,2),FORMC(1,MSP,2))
	maxnlr = MAXNLN; maxnlc = 1
	endif
      if(     cxwf) then
	  allocate(FORMC(MAXNLN,MSP,2))
	  maxnlc = MAXNLN
	  if(nn2wf) then
 		allocate(FORML(MAXNLN,MSP,2))
		maxnlr = MAXNLN
	      else
	    	allocate(FORML(1,MSP,2))
		maxnlr = 1
	  endif
	endif
      nforml1 = MAXNLN; nforml2 = MSP

      allocate(AFRAC(MXPEX,MXPEX,2,MSP))

      if(NKMAX>=0) then
	allocate(BPHASE(max(NKMAX,1),MSP))
	BPHASE(:,:) = 0.
	endif
      endif
 
C
C   READ IN MAIN INPUT VARIABLES
C   ----------------------------
C
	CALL READIN(WOUT,WOUT0,NF0,NF,STREN,ISNONO,
     X     NAME,MASS,NEX,QVAL,RMASS,HP,COPY,EXTRA,IOFAM,GIVEXS,
     X     JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X  PWFLAG,FORMF,FORML,FORMC,PTYPE,CPSO,EPS,MMXQRN,FILE,ffreal,
     X  HPOT,LAMBDA,MATRIX,MEK,NIX ,NSP, QNF,D0,BE,AFRAC,BPHASE,
     x  ICTO,ICFROM,KIND,QQ,IREM,KPCORE,BETAR,BETAI,JMAX,RMAX,REV,
     x  NFI,TWOW,OPN,NPRIOR,NPOST,REW14,LOCAL,NIB,LTRANS,FORMDEF,
     x  MLCALL,BPROJ,XA,XB,XP,XQ,SPINTR,FPT,NKP,VARYL,LSHAPE,LDEP,
     x  NLL,LOCF,KLT,COUPLE,ICOM,ICOR,NCP,IEXCH,GPT,POTCAP,INFILE,
     x  NSA, NJA, XCOEF, DISCIT, DISC8, LAMAX,NBINS,EMID,KMINX,NKBIN,
     x  NORBIN,STOP19,WID,CENTR)
!     if(.not.cxwf) write(KO,*) 'WF:',FORML(1:5,2,1)
!     if(.not.cxwf) write(KO,*) 'PT:',FORML(1:5,2,2)
!     if(     cxwf) write(KO,*) 'WF:',FORMC(1:5,2,1)
!     if(     cxwf) write(KO,*) 'PT:',FORMC(1:5,2,2)
	if(NOSOL) WOUT=.false.
	FALLOC = WOUT.or.ITER>0
 	 write(48,*) 'WOUT,ITER,FALLOC =',WOUT,ITER,FALLOC
	LOCFIL = .not.rterms 
	LOCFIL = LOCFIL .and. .not.FCWFN   ! eventually fix PWISER!
	LOCFIL = LOCFIL .and. NF>0		! set some efficiency minimum
	LOCFIL = LOCFIL .and. BPM==0.and.NFUS==0
     x		        .and. melfil==0 .and. .not.STOP19
	LOCFIL = .false.			! ALWAYS SEEMS TO SLOW, NOT SPEED UP!
	if(LOCFIL) write(KO,*)  ' Putting FORMF in file 19'
	 NLOC(:,:) = 0.0
       DO CP=1,NCP
       IF((KIND(CP).EQ.7.OR.KIND(CP).EQ.8).AND..NOT.OPEN12) THEN
	   inquire(12,opened=op)
	   if(op) then
	     WRITE(KO,*) ' FILE ',12,' ALREADY OPEN!!?'
	     endif
           OPEN(12,ACCESS='SEQUENTIAL',STATUS='SCRATCH',
     X           FORM='UNFORMATTED')
           OPEN12 = .TRUE.; written(12) = .true.
           ENDIF
	ENDDO
        T = TM(I)
        if(final)    write(KOI,*) 'Finished all Couplings @ ',real(T),NF
!        write(KOI,*)'NF0,NF,OPEN12,NCHAN,EPS,MMXQRN,NSP,NNN,
!     X               NFI,TWOW,OPN,REW14,NLL,IEXCH,NCP,NSA,NJA,
!     X		     XCOEF,DISCIT,DISC8,LAMAX='
!     X              ,NF0,NF,OPEN12,NCHAN,EPS,MMXQRN,NSP,NNN,
!     X               NFI,TWOW,OPN,REW14,NLL,IEXCH,NCP,NSA,NJA,
!     X		     XCOEF,DISCIT,DISC8,LAMAX

	if(final) then
	if(symm) write(ko,'(/'' Symmetric Hamiltonian'')')
	if(.not.symm) write(ko,'(/'' Non-symmetric Hamiltonian'')')
        if(PLANE>0) write(ko,'(/'' Plane wave treatment ='',i4)') PLANE
	endif
	ryev = 1.
      	if(RMASS(PEL)>1e-5)ryev = 1./(FMSCAL*RMASS(PEL))
	ipmax = LAMAX+1
!	write(KO,*) ' ipmax =',ipmax

	if(LOCFIL) then
	  DO I=1,MAXM
	  open(19,form='unformatted',status='scratch')
	  write(19) (FORMF(I,IF),IF=1,NF)
	  enddo
	  deallocate(FORMF)
	  allocate (FORMF(MLOC,1))
	endif
!	if(.not.rterms.and..not.mcalls) then
!	   write(48,*)  'FORML(1,1,1) =',FORML(1,1,1)
!	   write(48,*) 'FORML:',allocated(FORML); call flush(48)
!	   deallocate(FORML)
!	   allocate(FORML(1,1,2))
!           nforml1 = 1; nforml2 = 1
!	   write(48,*)  'FORML(1,1,1) =',FORML(1,1,1); call flush(48)
!	   endif
	
         ELAST = 0.
	CHANSI = CHANS
       if(final.and.(abs(NLAB(1))>1.or.num_energies>0)) then
	  call rewop(71)
	  write(71,'(''@legend ON'')')
	  endif
	do IC=200,210  ! see IFOUT and IFO in CRISS for the reason for these numbers
	  if(written(IC)) rewind IC
	  enddo
	LFAM=37
        if(written(LFAM)) rewind LFAM
	phaseadd(:) = 0.0
	phaslst(:) = 0.0
	LEN = 0
      DO 950 IEN = 1,4
      IF(num_energies==0.and.ELAB(IEN).lt.1e-8) GO TO 960
         IF(NLAB(IEN).EQ.0) NLAB(IEN) = 1
         DE = (ELAB(IEN+1) - ELAB(IEN)) / abs(NLAB(IEN))
         NLEN=abs(NLAB(IEN))
         if(num_energies>0) NLEN=num_energies
      DO 950 ILEN=0,NLEN
       if(num_energies==0) then
           IF(ILEN.GT.0 .AND. ELAB(IEN+1).lt.1e-20) GO TO 950
	     ENLAB = ELAB(IEN) + ILEN*DE   ! linear
           if(NLAB(IEN)<0) then            ! log
	     CF = (log(ELAB(IEN+1))-log(ELAB(IEN)))
     x               /(ELAB(IEN+1)-ELAB(IEN))
             CG =  log(ELAB(IEN)) - CF * ELAB(IEN)
	     ENLAB = exp(CF*ENLAB + CG)
	     endif
           IF(ABS(ENLAB-ELAST).LT.1E-8) GO TO 950
	 else
	   if(IEN>1.or.ILEN==0) go to 950
	   ENLAB = energy_list(ILEN)
	 endif
         ELAST = ENLAB
      if(final) WRITE(KO,1260) NAME(LIN,PEL),NAME(LIN,LAB),ENLAB
 1260 FORMAT('1',131('*'),/132('*')//8X,'INCOMING ',A8,' ;',
     X        4X,'LABORATORY ',A8,  ' ENERGY =',
     X G13.5,' MeV.' //       1X,131('*'),/1X,131('*') /)
!     if(.not.final) WRITE(KO,1261) NAME(LIN,PEL),NAME(LIN,LAB),ENLAB
!1261 FORMAT(/'  *** In ',A8,' ; Lab ',A8,' at',G13.5,' MeV.' )
      if(.not.final) WRITE(KO,1262) ENLAB
 1262 FORMAT(/'*** ',G13.5,' MeV.' )
	LEN = LEN+1
	CHANS = CHANSI
!	if(.not.final) CHANS = 0
	call flush(KO)
!	if(KO/=6) write(6,*)  'Scattering at lab energy ',real(ENLAB)
	
	if(allocated(finishedcc)) finishedcc(:) = .false.
C
      EOFF= - QVAL(LAB) + ENEX(1,LAB,LEX) + ENEX(2,LAB,LEX)
      ecmrat = MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
      ENCOM = ENLAB * ecmrat
        if(rela=='b') then
         T = MASS(2,LAB)+MASS(1,LAB)
         X = ENLAB*2.0*MASS(3-LIN,LAB) / T**2 / AMU
         ENCOM = T * (sqrt(1.0+x)-1.0) * AMU
         if(abs(x)<0.1) ENCOM = T*(X/2.-X**2/8. + X**3/16.) * AMU
         write(48,*) 'ECMNR,T,X,ENCOM =',ENLAB * ecmrat,T,X,ENCOM
        endif
      ETOTAL = ENCOM + EOFF
      etarat = ETACNS * MASS(2+1,pel)*MASS(2+2,pel) * SQRT(RMASS(pel))
C
C    CHANNEL ENERGIES
C    ----------------
! ----------------------------- AMoro -------------------
      if (.not.allocated(rener)) allocate(rener(mxp,mxx))
      rener(:,:)=0
! -------------------------------------------------------

      MAL1 = MIN(LMAX,INT(JTMAX+20.)) + 1
            XLMAX = MAL1 - 1
	hktarg = min(hktarg,0.20d0)
	FAIL = .false.
	BEST = HCM
      DO 145 IC=1,NCHAN
      NA = NEX(IC)
! AMoro
      sinv=amu**2*(mass(1,lab) +mass(2,lab))**2
     &   + 2*amu*mass(3-lin,lab)*enlab
      ecmi=sqrt(sinv)-(mass(1,lab) +mass(2,lab))*amu- eoff ! total initial kinetic energy
!      write(*,*) 'sqrt(sinv),Ecmi(rel),eoff=',sqrt(sinv),Ecmi,eoff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO 145 IA=1,NA
      IT = ITC(IC,IA)
         ECMC(IC,IA) = ETOTAL+QVAL(IC) - ENEX(1,IC,IA)-ENEX(2,IC,IA)
	GAM(IC,IA) = 1.0
!@@
         IF(RMASS(IC).lt.1e-5) then
            K(IC,IA) = ECMC(IC,IA)/HBC
         ELSE
            EE = FMSCAL*RMASS(IC) * ABS(ECMC(IC,IA))
            K(IC,IA) = SQRT(EE)
           if(rela/=' ') then
           if(rela=='a')
     x      X = (1. + ECMC(IC,IA)/(2.*AMU*RMASS(IC)))
     x        / (1. + 2*ECMC(IC,IA)/((MASS(1,IC)+MASS(2,IC))*AMU))
           if(rela=='b')
     x      X = (1. + ECMC(IC,IA)/(2.*AMU*MASS(LIN,IC)))
     x        / (1. + ECMC(IC,IA)/((MASS(1,IC)+MASS(2,IC))*AMU))
            K(IC,IA) = sqrt(EE*X)
            ECMC(IC,IA) = ECMC(IC,IA)*X
! AMoro
	  if (rela=='c'.or.rela=='3d') then
       mp       =mass(1,ic)*amu+ENEX(1,IC,IA)
       mt       =mass(2,ic)*amu+ENEX(2,IC,IA)
!       ecmnr     =ecmi + QVAL(IC) - ENEX(1,IC,IA)-ENEX(2,IC,IA)
       k(ic,ia) =sqrt(abs(triangle(sinv,mp**2,mt**2)))/2/sqrt(sinv)/hbc
       X        =k(ic,ia)**2/ee
       ep=(sinv + (mp**2-mt**2))/2./sqrt(sinv)      ! total energy of projectile in CM
       et=(sinv - (mp**2-mt**2))/2./sqrt(sinv)      ! total energy of target     in CM
       rener(ic,ia) = ep*et/(ep+et)/amu             ! reduced energy (amu)
! Relatistic energy
       erel = k(ic,ia)**2/(FMSCAL*rener(IC,ia))  ! 3D and LEA convention
! Effective energy for Schrodinger-like equation
       ecmc(ic,ia ) = k(ic,ia)**2/(FMSCAL*rmass(ic))

       gam(ic,ia)= rener(ic,ia)/rmass(ic) ! Erel/mu_rel
       eta(ic,ia)=gam(ic,ia)*fmscal*rmass(ic)*MASS(2+1,IC)*MASS(2+2,IC)
     x   *coulcn/k(ic,ia)/2
      endif ! 3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



           endif
         ENDIF
!        K(IC,IA) = SQRT(FMSCAL*RMASS(IC) * ABS(ECMC(IC,IA)) )
!@@
         ETA(IC,IA) = ETACNS * MASS(2+1,IC) * MASS(2+2,IC)
     X                         * SQRT(RMASS(IC)/ ABS(ECMC(IC,IA)))
        if(IC.eq.PEL.and.IA.eq.EXL) then
         if(mod(PLANE,2)==1) ETA(IC,IA) = 0.  ! elastic channel
        else
         if(    PLANE/2 ==1) ETA(IC,IA) = 0.  ! nonelastic channels
        endif
         SIGR(IC,IA) = 0.0
         NSB = 2*JEX(1,IC,IA) + 1.1
         NJB = 2*JEX(2,IC,IA) + 1.1
            MAM = NSA * NJA * NSB * NJB
            MMXF = MAX(MMXF,MAM)
!		write(161,*) ic,ia,mam,mmxf
!		call flush(161)
 	t = hcm*K(IC,IA)
	NAMEV = 'OK'
	if(t > hktarg) then
	  NAMEV = 'NOT OK'
	  write(ko,1445) ENLAB,IT,t,NAMEV,hktarg
 1445 format(/' Accuracy analysis at',f10.3,' MeV in ch.',i3,':',
     x        '  Elastic h*k =',f6.3,' so ',a8,' compared with',f6.3/)
          endif
	FAIL = FAIL .or. t > hktarg*1.2
	  BEST = min(BEST,hktarg/K(IC,IA))
  145   CONTINUE
	if(FAIL) then
	  write(KO,146) HCM,BEST
  146 	format(' Step size HCM',f8.5,' is too large: REDUCE TO <',
     x              f8.5,', so STOP!')
	  stop
	 endif
	MAXF = MMXF
      if(.not.allocated(SIGFUS)) allocate(SIGFUS(MAXF,3))
C
!      IF(FCWFN.and.final.and.JTMAX>100.) THEN
      IF(          final.and.JTMAX>50.) THEN
C              give warning of limits for Coulomb excitations:
       T = 2. * ETA(PEL,EXL) * 180./PI
        TH1 = T / (JTMAX+0.5)
            JTOTAL = -1.
	    X = -1
           IF(RASYM.lt.-0.01) THEN
		X = abs(RASYM)
                JTOTAL = T/(-RASYM)
                RASYM = max(JTOTAL/K(PEL,EXL),RMATCH+10*HCM)
		RS = RASYM
           ELSE IF(RASYM.gt.0.01) THEN
		RS = RASYM
	   ELSE
		RS = RMATCH
           ENDIF
       TH2 = T / (K(PEL,EXL) * RS)
!@       if(NGAIL>0) TH2 = 0.
C       IF(TH1.gt.THMIN .or. TH2.gt.THMIN)
            if(abs(ETA(PEL,EXL))>1e-10) WRITE(KO,162) TH1,TH2
  162  FORMAT(/'  NOTE: Coulomb excitation cut off ',
     X 'below',F8.3,' deg by JTMAX',:,', and below'
     X ,F8.3,' deg by Rmax'/)
!@        IF(X.gt.0..and.NGAIL==0) WRITE(KO,163) TH2,RS,JTOTAL
        IF(X.gt.0.) WRITE(KO,163) TH2,RS,JTOTAL
  163  FORMAT(' To get scattering to',F5.1,' deg, integrating to',F7.0,
     X   ' fm., and JTMAX should be at least',F6.0/)
!@        IF(X.gt.0..and.NGAIL>=0) WRITE(KO,164) X,JTOTAL
!@  164  FORMAT(' To get scattering to',F5.1,' deg',
!@     X   ', JTMAX should be at least',F8.0/)
 
	RASYMAX = RASYM
       ENDIF
C
C    COULOMB FUNCTIONS
C    -----------------
	DERIV = rterms
	 
	 CHL(:,:,:) = 0.
      DO 161 IC=1,NCHAN
      NA = NEX(IC)
      DO 160 IA=1,NA
      IT = ITC(IC,IA)
!!       if(RMASS(IC).lt.1e-5) go to 160
         RMK = (M-1)*HP(IC) * K(IC,IA)
	 if(DERIV) RMK = (MRM-1)*HP(IC) * K(IC,IA)
         RMKD = (MD-1)*HP(IC)  * K(IC,IA)
         IF(ECMC(IC,IA).GT.0.0) THEN
         CALL PHASES(ETA(IC,IA),MAL1,CSIG(1,IT))
!@       IF (NGAIL>=1) go to 160
          T = K(IC,IA)
       IF (FCWFN.AND.PWFLAG(IC)) THEN
          CALL COULFG(RASYM*T,ETA(IC,IA),Z,XLMAX,CFG,CFG(1,2),
     *         CFG(1,3),CFG(1,4),1,0,I,M1)
           DO 149 L1=1,MAL1
             CF = CFG(L1,1)
             CG = CFG(L1,2)
             CHL(L1,IT,1) = CG+CI*CF
!				Derivatives wrt R:
             CF = CFG(L1,3) * T
             CG = CFG(L1,4) * T
             CHL(L1,IT,2) = CG+CI*CF
 149       CONTINUE
       ELSE
        CALL COULFG(RMK,ETA(IC,IA),Z,XLMAX,CFG,CFG(1,2),
     *              CFG(1,3),CFG(1,4),1,0,I,M1)
           DO 150 L1=1,MAL1
              CF = CFG(L1,1)
              CG = CFG(L1,2)
C 150         CHL(L1,IT,1) = CMPLX(CG,CF)
  150         CHL(L1,IT,1) = CG+CI*CF
	if(DERIV) then
!				Derivatives wrt R:
           DO  L1=1,MAL1
              CF = CFG(L1,3) * T
              CG = CFG(L1,4) * T
	      CHL(L1,IT,2) = CG+CI*CF
	   ENDDO
	else
        CALL COULFG(RMKD,ETA(IC,IA),Z,XLMAX,CFG,CFG(1,2),
     *              CFG(1,3),CFG(1,4),2,0,I,M1)
           DO 151 L1=1,MAL1
              CF = CFG(L1,1)
              CG = CFG(L1,2)
C 151         CHL(L1,IT,2)= CMPLX(CG,CF)
  151         CHL(L1,IT,2)= CG+CI*CF
	endif
       ENDIF
      ELSE
C              CLOSED CHANNELS:
        IE = 0
         R = (M-1)*HP(IC)
      CALL WHIT(ETA(IC,IA),R,K(IC,IA),ECMC(IC,IA),MAL1,CFG,CFG(1,2),IE)
           DO L1=1,MAL1
              CSIG(L1,IT) = 0.0
              CHL(L1,IT,1) = dcmplx(0d0,CFG(L1,1))
	   ENDDO
	if(DERIV) then
           DO  L1=1,MAL1
              CHL(L1,IT,2) = dcmplx(0d0,CFG(L1,2))
	   ENDDO
	else
         R = (MD-1)*HP(IC)
      CALL WHIT(ETA(IC,IA),R,K(IC,IA),ECMC(IC,IA),MAL1,CFG,CFG(1,2),IE)
           DO L1=1,MAL1
  	    CHL(L1,IT,2)= dcmplx(0d0,CFG(L1,1))
	    enddo
	endif
      ENDIF
C
  160 CONTINUE
  161 CONTINUE
C
C
C    CDCC OUTPUT (1)
C    -----------------
!	BEPROJ=0.
	KN=1; 
	if(MSP>0) BEPROJ = BE(KN,1)
	KPB = 1 !	partition assumed for deuteron
	ENEX0(:,:,:) = ENEX(:,:,:)
      if(CDCC/=0) then
	I0 = 0
	do IC=1,NCHAN
	  if(MASS(1,IC)<MASS(1,KPB)) then
		I0 = IC
	        go to 152
	   endif
	enddo
	write(KO,*) 'CDCC: Could not find core partition! Stop'
	stop
152	continue
	do 153 KN=1,MSP
          if(abs(AFRAC(ITC(KPB,1),ITC(I0,1),1,KN))<1e-9) go to 153
          BEPROJ = BE(KN,1)
	  SPINV = QNF(10,KN)*0.5
	  PARV = 1
	  go to 154
153	continue
	write(KO,*) 'CDCC: Could not find projectile g.s! Stop', KPB,I0
	stop
154	continue
!				Try to make bin phases continuous:
        	LAST=-1; J2LAST=-1;KNLAST=-1
	do 156 ib=NBINS,1,-1
	 KN = NKBIN(1,IB)
	 L = QNF(9,KN)
	 J2 = QNF(11,KN)
	 if(ib==NBINS.or.L/=LAST.or.J2/=J2LAST) then
	   BSIGN(IB) = 1
	  else
	   I = nint((BPHASE(NKBIN(2,IB),KN)-BPHASE(1,KNLAST))/PI)
	     BSIGN(IB) = (-1)**I 
	     if(I/=0) 
     X	     write(KO,155) ib,-I,BSIGN(IB) !,nkbin(2,IB),max(NKMAX,1)
155          format('  Increase phase of bin ',i4,' by pi*',i3,
     X			', so  overall sign change = ',i2,2i6)
	     BPHASE(1:NKBIN(2,IB),KN) = BPHASE(1:NKBIN(2,IB),KN) - PI*I
	  endif
	 LAST = L
	 J2LAST = J2
	 KNLAST = KN
156	continue    

	MASSV = MASS(1,KPB)-MASS(1,I0)
	CHRGV = MASS(1+2,KPB)-MASS(1+2,I0)
	NAMEV = 'Valence'
	if(abs(MASSV-1.)<1e-2.and.abs(CHRGV-0.)<1e-2) NAMEV='Neutron '
	if(abs(MASSV-1.)<1e-2.and.abs(CHRGV-1.)<1e-2) NAMEV='Proton  '
	if(abs(MASSV-2.)<1e-2.and.abs(CHRGV-1.)<1e-2) NAMEV='Deuteron'
	if(abs(MASSV-3.)<1e-2.and.abs(CHRGV-1.)<1e-2) NAMEV='Triton  '
	if(abs(MASSV-4.)<1e-2.and.abs(CHRGV-2.)<1e-2) NAMEV='Alpha   '
	IBIN(1) = 0
	BSIGN(0) = 1
	NNJMAX = 0
	do 159 IA=2,NEX(KPB)
	 NNJMAX = max(NNJMAX, nint(2.*JEX(1,KPB,IA)+1.))
	 IB = 0
	 do 157 IB=1,NBINS
	 KN = NKBIN(1,IB)
          if(abs(AFRAC(ITC(KPB,IA),ITC(I0,1),1,KN))>1e-9) go to 158
157	 continue
	 write(KO,*) 'CDCC: Could not bin for projectile state ',IA
!	 stop
	 IB = 0
	 KN = 1
158	 IBIN(IA) = IB
	 CHSIGN(IA) = BSIGN(IB)
!	 write(KO,*) 'Excited state ',ia,' has bin ',ib,
!     X		' with phase ',BSIGN(IB)
         ENEX0(1,KPB,IA) = -BE(KN,1)+BEPROJ
!	 write(48,*) ' State ',KPB,IA,' has energy ',ENEX0(1,KPB,IA)
159	continue
         

        write(57,'(a80)') HEADNG
        write(57,'(F12.4,4F8.4,i2,f8.4)') ENLAB,BEPROJ,1d0/FMSCAL,
     x                                    COULCN,0.,0,0.
	write(57,'(4f8.4)') (MASS(I,KPB),I=1,2),MASS(1,I0),MASSV	! Masses
	write(57,'(4f8.1)') (MASS(I+2,KPB),I=1,2),MASS(1+2,I0),CHRGV  ! Charges
	write(57,'(4A8)') (NAME(I,KPB),I=1,2),NAME(1,I0),NAMEV  	! Names
	write(57,'(4f8.1)') (JEX(I,KPB,1),I=1,2),JEX(1,I0,1),SPINV  	! Spins
	write(57,'(4i8)') (sign(1,BAND(I,KPB,1)),I=1,2),
     X		sign(1,BAND(1,I0,1)),PARV  	! Parities
	write(57,'(4i4)') NBINS,NKMAX,NEX(KPB)-1,NNJMAX
	write(57,'(i4,2f8.4)') NANGL,THMIN,THINC
	do IA=1,NBINS
	  KN = NKBIN(1,IA)
          write(57,'(i2,f4.1,3f8.4,3i4)') QNF(9,KN),QNF(11,KN)*0.5,
     X		 EMID(IA),(KMINX(I,IA),I=1,2),NKBIN(2,IA),KN,NORBIN(IA)
	  write(57,'(10f8.4)') (BPHASE(I,KN),I=1,NKBIN(2,IA))
	  enddo
	call flush(57)
	endif

C
      DO 165 I=1,NSA*NJA
  165 SIGFUS(I,1) = 0.0
      DO 166 I=1,3
      TOTFUS(I) = 0.0
      if(NFUS>0) CORFUS(I,:) = 0.0
      SIGEL(I) = 0.
      SIGTOT(I) = 0.
  166 SIGT(I) = 0.0
      DO 167 I=1,NTHRESH
  167 THRJ(I) = 0.0
	inquire(10,opened=op)
!	if(op) close(10,status='delete')
	if(.not.op) open(10,form='unformatted',status='scratch')
      JSWITCH = 0.0
      SSWITCH = 0.0
      ALOSSM = 0.0
      JTEST = 0
      TSTD = 1.0
      RENORM = 1.0
      IF1 = 1
      IF(ITCM.EQ.1) IF1 = 2
      IF2 = 1+NFUS
       if(SMATS>=2) write(38,*) ITCM,IF1,IF2
	REWIND 10

        allocate (FNC(NLN,NSA),WNM(NLN,NSA,2))
	if(VEFF/=0) then
      DO 170 I=1,NLN
      FNC(I,1) = 0.0
  170 WNM(I,1,2) = 0.0
       open(15,access='sequential',status='scratch',form='unformatted')
      REWIND 15
      DO 180 IMA=1,NSA
  180 WRITE(15) (FNC(I,1),I=1,NLN)
      DO 190 IMA=1,NSA
  190 WRITE(15) (WNM(I,1,2),I=1,NLN)
      written(15) = .true.
      endif
	if(INITWF/=0) rewind abs(INITWF)
	INITWFE = .false.

C
C    FIND EACH J-TOTAL/PARITY COMBINATION
C    ------------------------------------
      JPSET = 0
      DONE = 0
      LJMAX = 0.0
      SMATL = SMATS
      if(.not.final) SMATL=0
      SHFL = SMATS.ge.6
      JCCSET = 0
      MAXB = 0
      nphases = 0
      SMALLS(:) = 0
      SMALLJ = .false.
      TIME0J = TM(I)
      CALL FLUSH(KO)
	if(TRENEG>0) write(89,211) 0,0,0,0
	if(SMATL>0) write(56,195) ENLAB,ETOTAL
	if(SMATL>0) write(156,195) ENLAB,ETOTAL
  195 format('# Fusion cross sections for ELAB=',f8.3,', ECM=',f8.3)
	written(56) = .true.
      NJDONE = 1

      DO 740 JBLOCK=1,NJ
C     IF(CHANS+LISTCC+SMATL.GE.3) WRITE(KO,*) JBLOCK,
C    #       JBORD(JBLOCK),JBORD(JBLOCK+1),JUMP(JBLOCK,1)
         JUMP(JBLOCK,2) = 0
         JUMP(JBLOCK,3) = 0
      IF(JUMP(JBLOCK,1) .LT. 1) GO TO 740
         JAP=JBORD(JBLOCK)
         IF(JBLOCK.GT.1.AND.JUMP(MAX(JBLOCK-1,1),1).GT.1)
     X                                      JAP=JAP+JUMP(JBLOCK,1)
         JAL=JBORD(JBLOCK+1) + 0.1
         IF(JBLOCK.LT.NJ.AND.JUMP(JBLOCK,1).EQ.1) JAL = JAL - 1
C
!      DO 730 JTOTAL=JAP,JAL,JUMP(JBLOCK,1)+Z
      NJTOTAL=NINT((JAL-JAP)/(JUMP(JBLOCK,1)+Z))
      DO 730 IJTOTAL=0,NJTOTAL
      JTOTAL = JAP + IJTOTAL*JUMP(JBLOCK,1)+Z
         SCHONP = .FALSE.
      JCOEF = (2*JTOTAL+1)* PI  * XCOEF / ABS(K(PEL,EXL))**2
	JUMP2 = JUMP(JBLOCK,1)
	if(IJTOTAL==NJTOTAL.and.JBLOCK<NJ) JUMP2=JUMP(JBLOCK+1,1)
	JUMPER = (JUMP(JBLOCK,1)+JUMP2)*0.5
        T = ETA(PEL,EXL)
        T4= K(PEL,EXL)
        L = JTOTAL
        RTURN =(T+SQRT(T**2 + L*(L+1d0)))/T4
        DONES=0
	FUSJ = 0.0
	TOTJ = 0.0
	SIGJ(:) = 0.0; FUSL(:)=0.0
	if(NFUS>0) CFUSJ(:) = 0.0
	JNLAST = JTOTAL
	SSWITCH = SWITCH
	if(JTOTAL>SINJMAX.and.abs(SINJMAX)>.1) SSWITCH = 1e6
      DO 720 IPARIT=2,3
         PARITY = (-1)**IPARIT
        if(PSET.ne.0.and.PARITY.ne.PSET) go to 720
       KS = 10
      call flush(ko)
	FUSJPI = 0.0
!	write(156,'("# J,pi =",f4.1,i3)') JTOTAL,PARITY
C
C   MAKE SET OF COUPLED EQUATIONS
C   -----------------------------
C
C
C      IF(LISTCC.GE.5) WRITE(KO,*) JTOTAL,PSIGN(PARITY+2)
C
C    PARTIAL WAVE SETS FOR EACH J-TOTAL/PARITY,  AND THEIR PARAMETERS
C    ----------------------------------------------------------------
       KINTL = -1
!       KINTL = 1000000
       IF(JSET.GT.0) KINTL = JSET-JPSET
	NSMALL = 2
!	if(SMALLJ.and.ITER==0.and.ITCM>1)  then	! Change from  CC  to DWBA
!	   IEX = 0
!	   ITER = 1
!	   write(KO,*) '  CHANGING TO DWBA!!!!!'
!	   ENDIF
	  

!				Find number of channels:
	CALL NUMCC(NCH,IEX,MINTL,MEL,JTOTAL,PARITY,JTMIN,KINTL,
     X        NEX,NCHAN,GIVEXS,PEL,EXL,LMAX,JEX,COPY,NFUSCH,
     X        RMASS,MAL1,ITC,IBLOCK,BAND,SMALLS,NSMALL,CHPRES)
!          WRITE(48,*) JTOTAL,PSIGN(PARITY+2),'NUMCC:',NCH,MINTL
        IF(NCH.EQ.0) go to 350 
        IF(MINTL==0.and.melfil==0) go to 350 
!        NICH = max(IEX,MEL)
        NICH = MEL
        if(ITER>0)  NICH = max(NICH,IEX)
        if(ITER>1.or.TWOW.or.IBLOCK<0)  NICH = NCH
        if(IBLOCK<0)  NICH = MEL
	 NICH = max(NICH,NFUSCH)
	if(WOUT0) NICH = NCH
        MMXCH=MAX(MMXCH,NCH)

       JCCSET = JCCSET + 1
      IF(.NOT.SCHONP) JUMP(JBLOCK,2) = JUMP(JBLOCK,2) + 1
       SCHONP = .TRUE.
      JUMP(JBLOCK,3) = JUMP(JBLOCK,3) + 1

!#	write(48,*)'NCH>MAXCH.or.NICH>MAXICH.or.IEX+2>MAXB'
!#	write(48,*) NCH,MAXCH,NICH,MAXICH,IEX,MAXB
!!?	FALLOC = WOUT.or.ITER>1 !!??? .or.IEX<NCH
	FALLOC = WOUT.or.ITER>0.or.IEX<NCH
!	 write(48,*) 'B: WOUT,ITER,FALLOC =',WOUT,ITER,FALLOC
!	if(.not.FALLOC) NFDEC=1
	if(.not.FALLOC) NFDEC=N*2*NICH
	if(NCH>MAXCH.or.NICH>MAXICH.or.IEX+2>MAXB) then
!					Deallocate arrays:
	if(MAXCH>0) then
	  deallocate (CLIST,NCLIST,NFLIST,CH,ECM,XS,INCOME,
     X    CFMAT,CGMAT,CRCRAT,CUTVAL,JVAL,JPROJ,JTARG,CHNO,CFF,
     X    LVAL,PART,EXCIT,INITL,BLOCKD,LL1,SMAT,EXCH,SAME,PSI,READWF)
	  if(.not.rterms) deallocate(SRC,SRATIO,SOMA,SPAD,
     X		FIM,FIMD,SIMPLE,EMPTY,SKIP,ferw)
	endif
!				Allocate arrays:
	MAXCH = max(NCH,NSA)
	MAXICH = NICH
	MAXXCH = max(MAXXCH,MAXCH)
	MAXCH = MAXXCH
	MFNL=NCH*NCH
	MAXB = IEX+2
!       NFDEC = 2 * N *(NICH*2+IEX*NICH)  + 4*MAXB*NCH
!       NFDEC = 2 * N *MAXB*NICH  + 4*MAXB*NCH
        NFDEC =  N *MAXB*NICH  !complex
!#	write(48,*) 'Makeset: NFDEC = ',NFDEC,N,MAXB,NICH,' set'
	if(NOSOL) then
	   NFDEC=1; MAXICH=max(IEX,MEL); MAXB=0; MAXIT=0
	   endif
	if(final) write(KO,*) 'Allocate arrays for ',MAXCH,
     X	 ' channels, of which ',MAXICH,' need wfs.'
	if(NFDEC>10 000 000) write(KO,*) 
     x		'Allocate NFDEC =',NFDEC,' in ',NFDEC*16/1000000,' MB'
      memuse = MAXCH*MAXCH*MCLIST*16/dble(1024**3)
      if(memuse>0.5)write(KO,'(a,3i7,a,f7.3,a)')
     &'Allocate CLIST(',MAXCH,MAXCH,MCLIST,') in ',memuse,' GB' 
      
      allocate(ECM(MAXCH,3),XS(MAXCH),INCOME(MAXCH))
      allocate(CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH))
      allocate(CUTVAL(MAXCH),BLOCKD(MAXCH),LL1(MAXCH),READWF(MAXCH))
      allocate(LVAL(MAXCH),PART(MAXCH,3),EXCIT(MAXCH,3),INITL(MAXCH))
      allocate(JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH),EXCH(MAXCH,MAXCH))
      allocate(CLIST(MAXCH,MAXCH,MCLIST),CFF(MAXCH,MAXCH,IPMAX))
      allocate(NCLIST(MAXCH,MAXCH),NFLIST(MAXCH,MAXCH,MCLIST))
      allocate(SMAT(MAXCH),CHNO(MFNL,6),CH(2,MAXCH),SAME(MAXCH))
      if(.not.rterms) then
      	allocate(SRC(MAXN,MAXCH),PSI(MAXN,MAXICH),SRATIO(MAXCH))
	allocate(SOMA(MAXCH,2),SPAD(MAXCH,MAXIT+2,MAXIT+2))
      	allocate(ferw(NFDEC))
      	allocate(FIM(MAXB,MAXCH),FIMD(MAXB,MAXCH))
      	allocate(SIMPLE(MAXCH),EMPTY(MAXCH),SKIP(MAXCH,2))
	else
	allocate(PSI(1,1))
	endif
!	else
!	write(KO,*) 'Have arrays for ',MAXCH,' channels, of which ',
!     X	 MAXICH,' need wfs.'
	endif

!#	write(48,*) 'Makeset: NFDEC = ',NFDEC,N,MAXB,NICH,' proceed'

      CALL CCSET(JTOTAL,PARITY,ETOTAL,JTMIN,KINTL,
     X        NEX,NCHAN,GIVEXS,QVAL,ENEX,PEL,EXL,LMAX,JEX,ECM,ECMC,
     X        LVAL,JVAL,PART,EXCIT,COPY,JPROJ,JTARG,CUTL,CUTR,CUTOFF,
     X        HP,RMASS,INCOME,BLOCKD,LUSED,MAL1,LJMAX,
     X        MINTL,INITL,ITC,IBLOCK,IEX,NCH,NCHPART,CHBASE,
     X        K,ETA,CHL,NAME,BAND,N,ISOCEN,RTURN,SMALLS,NSMALL,DROPPED)
!		write(6,*)  ' CCSET: IEX,MAXB =',IEX,MAXB
!          WRITE(48,*) JTOTAL,PSIGN(PARITY+2),'CCSET:',NCH,MINTL,IEX
C
C    READ IN ANY J/P-DEPENDENT POTENTIALS
C    ----------------------------------
      IB = INT(JTOTAL) + 1
      DO 240 JF=1,NF0
        JFT = PTYPE(2,JF)
          NFL = PTYPE(5,JF)
        IF(JFT.LE.7 .AND. NFL.GT.0.and.NFL<30) THEN
	  OPEN(NFL,ACCESS='DIRECT',FORM='FORMATTED',RECL=72,
     X			STATUS='OLD',file='potent25')
          IF(NFL.GE.24) IB = INT(JTOTAL) + 1
          IF(NFL.LT.24) IB = IPARIT-1
         NR = (N-1)/3 + 1
         IR0 = (IB-1)*NR
       WRITE(KO,*) 'Read in form at ',JF,' for block ',
     X               IB, ' with records from card ',IR0+1
         DO 230 IR=1,NR
           I0 = (IR-1)*3
230       READ(NFL,REC=IR+IR0,FMT=232)
     X                  (FORMF(I+I0,JF),I=1,MIN(3,N-I0))
232    FORMAT(6E12.4)
        IF(TRENEG.GT.0) WRITE(KO,235) ((I-1)*HCM,FORMF(I,JF),I=1,M,MR)
235         FORMAT(5(1X,F5.1,':',2F9.4,1X))
	close(NFL)
       ENDIF
240   CONTINUE
C
      CLIST(:,:,:) = 0.0
      NCLIST(:,:) = 0
      NL = 0
C
C    OPTICAL POTENTIALS FOR EACH CHANNEL & TENSORS BETWEEN CHANNELS
C    --------------------------------------------------------------
!	EPS=-1		! DEBUG ONLY
      ICH = 0
           IF(LISTCC.GE.1) WRITE(48,*) JTOTAL,PSIGN(PARITY+2)
      DO 325 IC=1,MXP
	NEXK = NCHPART(IC)
	I0 = CHBASE(IC)
	 repeat = .false.
         DO 302 C=I0+1,I0+NEXK
         IA = EXCIT(C,1)
         KP = CPOT(IC,IA)
	   do JF=1,NF0
	   if (PTYPE(1,JF)==KP) then
	   if (PTYPE(5,JF)==40) then
	     KP = nint(FORMDEF(IPARIT-1,JF))   ! parity-dependent choice of KP
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses parity-dep potential: ',KP
	   endif
	   if (PTYPE(5,JF)==41) then
	     L = LVAL(C)
	     KP = nint(FORMDEF(min(L,6)+1,JF))   ! L-dependent choice of KP
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses L-dep potential: ',KP
	   endif
	   if (PTYPE(5,JF)==42) then
	     L = int(JTOTAL)
	     KP = nint(FORMDEF(min(L,6)+1,JF))   ! J-dependent choice of KP
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses J-dep potential: ',KP
	   endif
	   endif
	   enddo
         DO 302 C2=I0+1,I0+NEXK
            IB = EXCIT(C2,1)
            KPB = CPOT(IC,IB)

	   do JF=1,NF0
	   if (PTYPE(1,JF)==KPB) then
	   if (PTYPE(5,JF)==40) then
	     KPB = nint(FORMDEF(IPARIT-1,JF))   ! parity-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses parity-dep potential: ',KPB
	   endif
	   if (PTYPE(5,JF)==41) then
	     L = LVAL(C2)
	     KPB = nint(FORMDEF(min(L,6)+1,JF))   ! L-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses L-dep potential: ',KPB
	   endif
	   if (PTYPE(5,JF)==42) then
	     L = int(JTOTAL)
	     KPB = nint(FORMDEF(min(L,6)+1,JF))   ! J-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses J-dep potential: ',KPB
	   endif
	   endif
	   enddo

         DO 300 JF=1,NF0
!     IF(LISTCC.GE.2) WRITE(KO,*) ' TRY ',JF,' @ ',(PTYPE(II,JF),II=1,6)
            IF((PTYPE(1,JF).NE.KP .and. PTYPE(1,JF).ne.KPB) .OR.
     X         PTYPE(3,JF).LT.0 .OR. PTYPE(4,JF).LT.0)    GO TO 300
            JFT = PTYPE(2,JF)
	    JFTT = JFT
	    NFD = PTYPE(4,JF)
	   SCALEL = 1d0
	   if(LDEP(JF)>0) then
	     T = JTOTAL
	     if(LDEP(JF)==2) T = (LVAL(C)+LVAL(C2))/2
	     RH = (T- VARYL(1,JF))/VARYL(2,JF)
	     EE = DDEXP(-RH)
	     if(LSHAPE(JF)==0) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)
	     else if(LSHAPE(JF)==1) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)**2
	     else if(LSHAPE(JF)==2) then
	 	SCALEL = DDEXP(-RH*RH)
	     endif
	    if(LISTCC>50) write(KO,*) 'JL-dependence for form ',JF,
     X		': LDEP,VARY,T,RH,EE,SCALEL=',
     X             LDEP(JF),VARYL(1:2,JF),T,RH,EE,SCALEL
	   endif
      IF(EXCIT(C,3)/=EXCIT(C2,3) .AND.(JFT==10.OR.JFT==12.OR.JFT==14)
     X .OR.
     X   EXCIT(C,2)/=EXCIT(C2,2) .AND.(JFT==11.OR.JFT==13))
     X  GO TO 300
C**********************************************************************
      IF(JFT.GE.10 .AND. JFT.LE.15 .AND.
     X   (EXCIT(C,2).NE.EXCIT(C2,2) .AND. EXCIT(C,3).NE.EXCIT(C2,3)))
     X  GO TO 300
      IF(JFT.LT.10 .AND.
     X   (EXCIT(C,2).NE.EXCIT(C2,2) .OR. EXCIT(C,3).NE.EXCIT(C2,3)))
     X  GO TO 300

       T = 1.0
         if(JFT==0) then
           if(IC.eq.PEL.and.IA.eq.EXL) then
            if(mod(PLANE,2)==1) T = 0.  ! elastic channel
           else
            if(    PLANE/2 ==1) T = 0.  ! nonelastic channels
           endif
         endif
       KK = PTYPE(3,JF)
       IF(KK.eq.7) THEN
C                       Take K=7 as an inelastic monopole
          IF(IA.eq.IB) GO TO 300 ! exclude diagonal monopole
          KK = 0
       ELSE if(KK==0) then
C			Normal KK=0 should exclude off-diagonal monopoles
C                       EXCEPT (March 2021) for spin.spin couplings!!! FIXED
	  if(C/=C2 .and. JFT/=8.and.JFT/=4) GO TO 300
           if(JFT==12) T = sqrt(2*JPROJ(C)+1.)   ! projectile couplings
           if(JFT==13) T = sqrt(2*JTARG(C)+1.)   ! target couplings
       ENDIF
C
       KKP = (-1)**KK
      if (JFTT.eq.10.and.BAND(1,ic,ia).ne.BAND(1,ic,ib)*KKP .or.
     x    JFTT.eq.11.and.BAND(2,ic,ia).ne.BAND(2,ic,ib)*KKP ) go to 300

      IF(JFT>=12.and.JFT<=13.and.PTYPE(3,JF)>0) THEN
C                           LOOK UP TABLE OF ALLOWED COUPLINGS for KK>0 or off-diagonal monopole (7)
         T = 0.0
         DO 294 IX=1,NIX
           IF(MATRIX(1,IX).NE.EXCIT(C,JFT-10) .OR.
     X        MATRIX(2,IX).NE.EXCIT(C2,JFT-10) .OR.
     X        MATRIX(3,IX).NE.PTYPE(3,JF) .OR.
     X        MATRIX(4,IX).NE.JF) GO TO 294
         T = T + MEK(IX)
294      CONTINUE
C
      ELSE IF(JFT>=14.and.JFT<=16) THEN		! Deal with later
         GO TO 302

      ELSE IF(JFT.EQ.17) THEN
	 if(PTYPE(3,JF)/=IC) go to 300
C		This last condition allows only the *first* JF slot.
C                           LOOK UP TABLE OF ALLOWED COUPLINGS :
         T = 0.0
         DO 298 IX=1,NIX
           IF(MATRIX(1,IX).NE.EXCIT(C,1) .OR.
     X        MATRIX(2,IX).NE.EXCIT(C2,1).OR.
     X        MATRIX(5,IX).NE.NFD) GO TO 298
         T = T + MEK(IX)
         KK = MATRIX(3,IX)
	 if(KK==7) KK=0
	 JFTT = 0
         IF(EXCIT(C,3).NE.EXCIT(C2,3)) JFTT=13
         IF(EXCIT(C,2).NE.EXCIT(C2,2)) JFTT=12
	 if(JFTT==0) JFTT = MATRIX(6,IX)+11
298      ENDDO
	if(.not.repeat) then
	 allocate (OO(NEXK,NEXK),EVEC(NEXK,NEXK))
	 OO(:,:) = 0d0
	endif
         JF0 = JF
	 repeat = .true.
       ENDIF
C
          S =  TENSOR(JFTT,LVAL(C),JPROJ(C),JVAL(C),JTARG(C),
     X                JTOTAL,LVAL(C2),JPROJ(C2),JVAL(C2),JTARG(C2),
     X                KK,JEX(3,IC,IA),JEX(4,IC,IA),
     X                 JEX(3,IC,IB),JEX(4,IC,IB),
     X                 MASS(3,IC),MASS(4,IC),PTYPE(5,JF))
	  S = S
C    X      * CI**NINT(JPROJ(C2)+JTARG(C2) - JPROJ(C) - JTARG(C))
     X      * CI**NINT( - ABS(JPROJ(C2)-JPROJ(C))
     X                   +     JPROJ(C2)-JPROJ(C)
     X                   - ABS(JTARG(C2)-JTARG(C))
     X                   +     JTARG(C2)-JTARG(C))
C           The above phase factors with JPROJ & JTARG etc.,
C           are there only because of definition of M(Ek) matrix element
!	if(JFT<=11.or.JFT>=18) 
	if(JFT/=17) S = S * CI**(LVAL(C2)-LVAL(C)) * SCALEL
	if(repeat.and.JFT==17) then
C					Use couplings as deformation lengths
	 OO(C-I0,C2-I0) = OO(C-I0,C2-I0) + T*REAL(S)
!	 OO(C2-I0,C-I0) = OO(C2-I0,C-I0) + T*REAL(S)
           IF(LISTCC.GT.1.and.abs(T)>EPS) WRITE(KO,1330) C,C2,JF,
     X	(PTYPE(II,JF),II=1,6),KK,LDEP(JF),JFTT,T,S,SCALEL
!         IF(LISTCC>2.and.abs(T)>EPS) write(KO,*) OO(C-I0,C2-I0),C,C2,I0
!          IF(LISTCC>2.and.abs(T)>EPS) write(KO,*) OO
	endif
	   NC = NCLIST(C,C2)+1
	if(abs(S*T)>EPS.and.JFT/=17) then
	   if(NC>MCLIST) then
             write(KO,*) 'For channels ',C,C2
             write(KO,*) 'Need NC=',NC,' > MCLIST=',MCLIST,' !!!'
             write(KO,*)'So far use forms ',NFLIST(C,C2,:),' need ',JF
             call check(NC,MCLIST,30)
             stop
             endif
	   CLIST(C,C2,NC) = S*T
	   NFLIST(C,C2,NC) = JF
	   NCLIST(C,C2) =  NC
           IF(C.NE.C2.AND.S*T.NE.0.0) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1) WRITE(KO,1330) C,C2,JF,
     X    (PTYPE(II,JF),II=1,6),KK,LDEP(JF),NC,T,S,SCALEL
 1330      FORMAT(' For',I5,' from',I5,' by form',I5,' (',6I3,')',2i3,
     X	    i10,' get',F8.4,2F15.6,F8.4)
	 endif
  300    CONTINUE
  302  CONTINUE
! 					Deal with 14,15,16 now
         DO 308 C=I0+1,I0+NEXK
         IA = EXCIT(C,1)
         KP = CPOT(IC,IA)
         DO 308 C2=I0+1,I0+NEXK
            IB = EXCIT(C2,1)
            KPB = CPOT(IC,IB)
         DO 307 JF=1,NF0
!       IF(LISTCC.GE.2) WRITE(KO,*) 'To ',C,' from  ',C2,', ',
!     X	' TRY ',JF,' @ ',(PTYPE(II,JF),II=1,6)
            IF((PTYPE(1,JF).NE.KP .and. PTYPE(1,JF).ne.KPB) .OR.
     X         PTYPE(3,JF).LT.0 .OR. PTYPE(4,JF).LT.0)    GO TO 307
            JFT = PTYPE(2,JF)
       IF(JFT<14.or.JFT>17) GO TO 307		! Deal with 14,15,16 only
	    NFD = PTYPE(4,JF)
	   SCALEL = 1d0
	   if(LDEP(JF)>0) then
	     T = JTOTAL
	     if(LDEP(JF)==2) T = (LVAL(C)+LVAL(C2))/2
	     RH = (T- VARYL(1,JF))/VARYL(2,JF)
	     EE = DDEXP(-RH)
	     if(LSHAPE(JF)==0) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)
	     else if(LSHAPE(JF)==1) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)**2
	     else if(LSHAPE(JF)==2) then
	 	SCALEL = DDEXP(-RH*RH)
	     endif
	    if(LISTCC>50) write(KO,*) 'JL-dependence for form ',JF,
     X		': LDEP,VARY,T,RH,EE,SCALEL=',
     X             LDEP(JF),VARYL(1:2,JF),T,RH,EE,SCALEL
	   endif
!      IF(LISTCC.GE.2) WRITE(KO,*) ' TRY2 ',JF,' for ',
!     X   EXCIT(C,3),EXCIT(C2,3),EXCIT(C,2),EXCIT(C2,2)
      IF((EXCIT(C,3)/=EXCIT(C2,3) .AND.JFT==14)  .OR.
     X   (EXCIT(C,2)/=EXCIT(C2,2) .AND.JFT==15)) GO TO 307
!      IF(LISTCC.GE.2) WRITE(KO,*) ' TRY3 ',JF,' * ',SCALEL
C**********************************************************************
!?      IF(JFT==16 .AND.
!?     X   (EXCIT(C,2).eq.EXCIT(C2,2) .or. EXCIT(C,3).eq.EXCIT(C2,3)))
!?     X  GO TO 307

       KK = PTYPE(3,JF)
       if(KK==7) KK=0
C	 I1 and I2 are 1(proj) or 2(targ) for the two transitions
       IF(JFT>=14.and.JFT<=15) THEN		
	 I1 = JFT-13			! Both = 1 for 14, 2 for 15
	 I2 = JFT-13
       ELSE IF(JFT==16) THEN
	 I1=1				! Projectile
	 I1=2				! Target
       ENDIF
!       IF(LISTCC.GE.2) WRITE(KO,*) ' TRY4 ',JF,KK,I1,I2,NFD
C                           LOOK UP TABLE OF ALLOWED COUPLINGS :
         DO 305 IX1=1,NIX
         DO 305 C1=I0+1,I0+NEXK
         DO 305 IX2=1,NIX
           LAM1 = +MATRIX(3,IX1)
           LAM2 = -MATRIX(3,IX2)
!       IF(LISTCC.GE.2) WRITE(KO,*) ' TRY5 ',IX1,IX2,' @ ',LAM1,LAM2,C1
           IF(JFT>=14.and.JFT<=15) LAM1=abs(LAM1)

           IF(MATRIX(1,IX1).NE.EXCIT(C,I1+1) .OR.
     X        MATRIX(2,IX1).NE.EXCIT(C1,I1+1) .OR.
     X        LAM1<0           .OR.
     X        MATRIX(5,IX1).NE.NFD) GO TO 305
         TH1 = MEK(IX1)

           LAM2 = -MATRIX(3,IX2)
           IF(JFT>=14.and.JFT<=15) LAM2=abs(LAM2)
           IF(MATRIX(1,IX2).NE.EXCIT(C1,I2+1) .OR.
     X        MATRIX(2,IX2).NE.EXCIT(C2,I2+1) .OR.
     X        LAM2<0           .OR.
     X        MATRIX(5,IX2).NE.NFD) GO TO 305
         TH2 = MEK(IX2)

	 IF(KK<abs(LAM1-LAM2) .or. KK>LAM1+LAM2) go to 305
!       IF(LISTCC.GE.2) WRITE(KO,*) ' Found ',IX1,IX2,' @ ',LAM1,LAM2,C1
	 T = TH1*TH2 * 0.5 			!0.5 from 1/2!
C
          S =  TENSOR2(JFT,LVAL(C),JPROJ(C),JVAL(C),JTARG(C),
     X                JTOTAL,LVAL(C2),JPROJ(C2),JVAL(C2),JTARG(C2),
     X                KK,LAM1,LAM2,JEX(3,IC,IA),JEX(4,IC,IA),
     X                 JEX(3,IC,IB),JEX(4,IC,IB),
     X                 MASS(3,IC),MASS(4,IC),PTYPE(5,JF))
	  S = S
     X      * CI**NINT( - ABS(JPROJ(C2)-JPROJ(C))
     X                   +     JPROJ(C2)-JPROJ(C)
     X                   - ABS(JTARG(C2)-JTARG(C))
     X                   +     JTARG(C2)-JTARG(C))
C           The above phase factors with JPROJ & JTARG etc.,
C           are there only because of definition of M(Ek) matrix element
	S = S * CI**(LVAL(C2)-LVAL(C)) * SCALEL

	   NC = NCLIST(C,C2)+1
	if(abs(S*T)>EPS) then
	   if(NC>MCLIST) then
		write(KO,*) 'For channels ',C,C2
		write(KO,*) 'Need NC=',NC,' > MCLIST=',MCLIST,' !!!'
		write(KO,*) 'So far forms ',NFLIST(C,C2,:),'; need ',JF
		call check(NC,MCLIST,30)
		stop
		endif
	   CLIST(C,C2,NC) = S*T
	   NFLIST(C,C2,NC) = JF
	   NCLIST(C,C2) =  NC
           IF(C.NE.C2.AND.S*T.NE.0.0) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1) WRITE(KO,1331) C,C1,C2,JF,
     X    (PTYPE(II,JF),II=1,6),KK,LDEP(JF),NC,TH1,TH2,S,SCALEL
 1331      FORMAT(' For',I5,' via',I5,' from',I5,' by form',I5,
     X		' (',6I3,')',3i3,' get',2F8.4,2F15.6,F8.4)
	 endif
  305    CONTINUE
  307    CONTINUE
  308  CONTINUE

	if(repeat) then  
           IF(LISTCC.GT.2) then
               write(KO,*) ' OO matrix:'
              do i=1,NEXK
               write(KO,3101) i,(OO(i,j) , j=1,NEXK)
	      enddo
 3101          format(1x,i3,12f10.3,:,/(4x,12f10.3))
	     endif
C
C    ALL-ORDER COUPLING POTENTIALS BETWEN CHANNELS
C    ---------------------------------------------
	call HDIAG(OO,NEXK,NEXK,0,EVEC,I)

           IF(LISTCC.GT.1) WRITE(KO,*) ' Diagonalised matrix',NEXK,I
           IF(LISTCC.GT.2) then
	  	write(KO,311)(OO(IE,IE),IE=1,NEXK)
311		format('  Eigenvalues are:'/(1x,10f8.4))
	 	endif
	IF = JF0
        do 322 ia=1,NEXK
	do 322 ib=ia,NEXK
	  c = ia+I0
	  c2= ib+I0
!	if(LISTCC>1) write(KO,312) c,c2,IF,NFD
!312	format(30X,'  All-order coupling ',i3,' <->',i3,' is at ',i3,
!    X		', deforming potential at ',i3)
               CALL CHECK(IF,MLOC,24)
C					The monopole is replaced completely!
!               PTYPE(4,NFD) = -JF0
	 FORMF(1:N,IF) = 0d0
	  do 320 I=1,N
	   R1 = (I-1)*HP(IC)
	   do 320 IE=1,NEXK
	   X = OO(IE,IE)*R4PI
	   R = R1 - X
	   S = FFC(R/HP(IC),FORMF(1,NFD),N)
            FORMF(I,IF)=FORMF(I,IF) + S * EVEC(IE,IA)*EVEC(IE,IB)
320	   continue
	L = 1
	if(C/=C2) L=2
C				Do forward (IE=1) and reverse (IE=2) couplings
	Do 321 IE=1,L
	if(IE==2) then
	 I = C
	 C = C2
	 C2 = I
	 endif
C						Put in phase factors:
	  S =     CI**(LVAL(C2)-LVAL(C))
	  
	   NC = NCLIST(C,C2)+1
	   if(NC>MCLIST) then
		write(KO,*) 'For channels ',C,C2
		write(KO,*) 'Need NC=',NC,'>MCLIST=',MCLIST,' !!!'
		write(KO,*) 'So far use forms ',NFLIST(C,C2,:),' need',IF
		call check(NC,MCLIST,30)
		stop
		endif
	   CLIST(C,C2,NC) = S
	   NFLIST(C,C2,NC) = IF
	   NCLIST(C,C2) =  NC
           IF(C.NE.C2) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1) 
     X	WRITE(KO,1330) C,C2,IF,(PTYPE(II,IF),II=1,6),NC,0,0,1.,S
321	continue
      IF(TRENEG*LISTCC.GE.1) THEN
            WRITE(KO,110) IF
110         FORMAT(/' Potential Form at',I4,' is')
            WRITE(KO,120) ((I-1)*HP(IC),FORMF(I,IF),I=1,N,MR)
120         FORMAT(5(1X,F5.1,':',2F9.4,1X))
           ENDIF
	 IF = IF+1
322	 continue
	deallocate(OO,EVEC)
	endif
	 
  325  CONTINUE

C
C    COUPLING COEFFICIENTS BETWEEN PAIRS OF PARTIAL WAVES
C    ---------------------------------------------------
      if(OPEN12) REWIND 12
      IF(REW14) REWIND 14
      DO 290 CP=1,NCP
            IF(.NOT.COUPLE(CP).OR.JTOTAL.GT.JMAX(CP)) GO TO 290
            IF(.NOT.LOCAL(CP).AND.ITER==0.and..not.rterms) GO TO 290
         IC1 = ICTO(CP)
         IC2 = ICFROM(CP)
            IN1 = MOD(ABS(KIND(CP))-1,2)+1
         IF(LISTCC.GT.0) WRITE(KO,1320) CP,ICTO(CP),ICFROM(CP),KIND(CP)
 1320    FORMAT(/' Coupling #',I4,' to',I4,' from',I4,' of KIND=',i2)

      CALL CPAIR(CP,IC1,IC2,KIND(CP),QQ(CP),KPCORE(CP),REV(CP),LOCAL(CP)
     X          ,ICOM,ICOR,LOCF(CP),NKP(1,CP),GPT(1,1,CP),XA(CP),XB(CP),
     X          XP(CP),XQ(CP),IREM(CP),FILE(CP),CHNO,ffreal(CP),KLT(CP),
     X    NCH,LVAL,JVAL,JPROJ,JTARG,JTOTAL,PART,EXCIT,JEX,COPY,BAND,
     X    QNF,AFRAC,CUTOFF,ALOSSM,ACC8,FPT(1,1,CP),NLL(CP),
     X    NKP(1,CP),BE(1,1),CLIST,NFLIST,NCLIST,DNL,EPS,NFORML1,NFORML2,
     X    HP,ITC,NEX(ICOR(CP,IN1)),ICFROM(CP),ICTO(CP),MATRIX,MEK,
     X    CPSO(CP),SPINTR,K,POTCAP(1,1,CP),LTRANS(1,CP),NLOC(1,CP),
     X    MLCALL(1,CP),MLCALL(2,CP),FORML,VPOT,VPOT(1,2),NIB(CP),MASS,
     X    FORMC)

  290 CONTINUE
      IF(NL.GT.0) JTEST = JTEST + 1
      X = ACC8*ALOSSM
      IF(JTEST.EQ.4 .AND.X.GT.1E-4.and..not.DRY) then
         WRITE(KO,291) LOG10(MAX(ALOSSM,1D0)),X*100.
  291 FORMAT('0***** WARNING : SOME LARGE L-TRANSFER GIVES ACCURACY LOSS
     X IN KERNELS OF',F5.1,' DIGITS,',/,
     X '0          SO MUST EXPECT ITS ERROR TO BE',   F10.5,' %',
     X '    MTMIN could be decreased'/)
	 if(X > 0.1) then
		write(6,'(//'' THIS ERROR IS TOO LARGE. STOPPING'')')
		stop
		endif
	 endif
      NFNL = MAX(NFNL,NL)
      NL = MIN(NL,MFNL)


C
C    CHECK COULOMB MATCHING!

      do ip=1,nvars  ! check if NOPOT for this J/pi:
      if(srch_kind(ip)==3.and.abs(JTOTAL-srch_jtot(ip))<.1
     x .and.abs(PARITY-srch_par(ip))<.1 .and. srch_nopot(ip)) goto 336
      enddo
       I = 0
       DO 335 C=1,NCH
          IC = PART(C,1)
          IA = EXCIT(C,1)
          R0 = (M-1)*HP(IC)
          T = MASS(3,IC) * MASS(4,IC) * COULCN / R0
           if(IC.eq.PEL.and.IA.eq.EXL) then
            if(mod(PLANE,2)==1) T = 0.  ! elastic channel
            else
            if(    PLANE/2 ==1) T = 0.  ! nonelastic channels
           endif

          CF = 0.0
	if(LOCFIL) then
	  rewind 19
	  do IB=1,M-1
	  read(19)
	  enddo
	  read(19) FORMF(1:NF,1)
	  endif
       DO 330 NC=1,NCLIST(C,C)
	JF = NFLIST(C,C,NC)
       if(LOCFIL) then
       R = CLIST(C,C,NC) * FORMF(JF,1) 
       else
       R = CLIST(C,C,NC) * FORMF(M,JF) 
       endif
       if(LAMBDA(JF)>1.and.FCWFN.AND.PWFLAG(IC)) 
     X   R = R - CLIST(C,C,NC) * STREN(JF)/R0**(LAMBDA(JF)+1)
        IF(ABS(R).GT.1E-5 .AND. ABS(HP(IC)-HPOT(JF)).GT.1E-5)
     X      WRITE(KO,*) '***** CHANNEL ',C,' USING POTENTIAL AT',JF,
     X          ', WHICH HAS STEP SIZE ',HPOT(JF),', NOT',HP(IC),'!!'
  330  CF = CF + R
       T = T - CF
       IF(final.and.ABS(T).GT.0.020) WRITE(KO,1340) C,IC,EXCIT(C,1),T,R0
       IF(ABS(T).GT.0.10 ) I = I + 1
 1340  FORMAT(/' ****** ERROR : CHANNEL',I4,' PARTITION',I3,' LEVEL',I3,
     X' HAS MATCHING DEFICIENCY OF',F8.4,' MEV. AT R =',F8.2)
  335  CONTINUE
       IF(I>0) then
C        STOP 'COULOMB MONOPOLES'
         write(6,*) ' BAD Coulomb deficiences! ',I !,', found ',real(CF)
         go to 350 
!        CALL ABEND(4)
       endif
336   if(TRENEG>=1) then
      write(89,210) JTOTAL,PARITY,NCH
210	format('#',f7.1,i3,i5,' : J,pi,NCH')
      DO 215 C1=1,NCH
	write(89,212) C1,LVAL(C1),JVAL(C1),JPROJ(C1),JTARG(C1),
     x                PART(C1,1),EXCIT(C1,1)
      DO 215 C2=1,NCH
      IN = NCLIST(C1,C2)
      write(89,211) C1,C2,IN
	do I=1,IN
	write(89,213) NFLIST(C1,C2,I),CLIST(C1,C2,I)
	enddo
211	format('#',4i5)
212	format('##',i4,' LJ:',i5,3f6.1,2i4)
213	format('#',i5,1p,2e14.6)
215   continue
	write(89,210) 0.,0,0
      endif
C
 1358 FORMAT(' ',F9.3,11F11.4/(11X,11F11.4) )
 
	CALL ASYMPTOPIA(NCH,ipmax,CI,STREN,LVAL,MASS,PART,
     X    ryev,LAMBDA,NF,CLIST,NFLIST,NCLIST,CFF)

       If(NOSOL) go to 350 
C
! C   Only use only CRCWFN if R-turn > |RMATCH| + GAP
!         GAP = 4.5
!         if(cutr.lt.0.) GAP=abs(cutr)
C   Only use only CRCWFN if R-turn > CUTOFF*HP(1)
         GAP = 4.5
         if(CUTR.lt.-1.) GAP = max(RTURN - CUTOFF*HP(1),1d0)

      FJSWTCH = FCWFN.AND.JTOTAL.GT.AJSWTCH
     X         .and.(RTURN.gt.RMATR+GAP .or.CUTOFF.gt.MRM*5/6)
      IF(FJSWTCH .and. JSWITCH.le.0.1) JSWITCH=JTOTAL
      IF(CDETR.gt.0) WRITE(KO,*) 'At J=',JTOTAL,
     X               ', Rturn =',RTURN,':',FJSWTCH
C
!	write(KO,*) 'LAMBDA(:) =',LAMBDA(1:MLOC),' so LAMAX =',LAMAX
!	write(106,*) 'STREN(:) =',STREN(1:MLOC)

        IF(FCWFN.and.RTURN.ge.max(RASYM,RMATCH)*1.5) then
            WRITE(KO,*) 'At J=',JTOTAL,', Rturn =',RTURN
            WRITE(KO,*) 'Turning point beyond R limit,',
     X                  ' STOPPING SOON.'
c            DONES = DONES+500
            DONES = DONES+50
c            GO TO 750
         ENDIF
        IF(.not.FCWFN.and.CUTOFF.gt.MD-20) THEN
            WRITE(KO,*) 'At J=',JTOTAL,', CUTOFF =',CUTOFF*HP(1)
            WRITE(KO,*) ' Lower Cutoff approaching RMATCH,',
     X                  ' STOPPING SOON.'
            DONES = DONES+500
c            GO TO 750
         ENDIF

      NSTEPD=0
      IF (FCWFN) CALL PWISER(LVAL,ECM,K,ETA,CLIST,NCLIST,NFLIST,
     *    PART,EXCIT,CHL,ACCRCY,RASYM,CFMAT,CGMAT,CRCRAT,NCH,M,MD,
     *    PWFLAG,LAMAX,STREN,LAMBDA,DERIV,FORMF,NF,MRM,GAP,
     *    ITC,SSWITCH,FJSWTCH,CDETR,SMATL.ge.3,NGAIL>0,NSTEPD,
     *    ALLPART)
      IF (FCWFN.and..not.PWFLAG(PEL)) then
	write(6,*) ' MUST HAVE PWFLAG SET IN ELASTIC CHANNEL!'
	stop
	endif
      IF (FCWFN.and.CDETR>0) then
	if(DERIV) then
	   write(KO,*) ' Find CRCWFN derivatives at',(mrm-1)*real(hcm)
	  else
	   write(KO,*) ' Find CRCWFN at',(m-1)*real(hcm),' and ',
     *				(md-1)*real(hcm)
	  endif
	endif
C
      WRITE(KS) JTOTAL,NCH,MINTL,IPARIT,JUMPER
      WRITE(KS) (LVAL(C),JVAL(C),PART(C,1),EXCIT(C,1),C=1,NCH)
      written(KS) = .true.
c
      rmorto =rmort
      DO 345 C=1,NCH
        CUTVAL(C) = CUTOFF
        IF(CUTL.lt.0) then
C       			Make lower cutoff on L(C) not JTOTAL
               CUTVAL(C) = MAX(abs(CUTL)*LVAL(C), CUTR/HP(1),0d0) + 1
        IC=PART(C,1); IA=EXCIT(C,1); T = ETA(IC,IA)
!        RTURN =(T+SQRT(T**2 + LVAL(C)*(LVAL(C)+1d0)))/K(IC,IA)
               if(CUTR.lt.0.) CUTVAL(C) = max(CUTVAL(C),
     *                  int(10.**(CUTR/(JTOTAL+1))*RTURN/HP(1)) )
!	write(6,*) 'RTURN, rirat=',RTURN,10.**(CUTR/(JTOTAL+1)),
!     X 	                           RTURN*10.**(CUTR/(JTOTAL+1)),
!     #         CUTVAL(C)*HP(1)
c        CUTVAL(C) = 5
! AMoro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(ECM(C,1).GT.0.) THEN
          RTURNC =(T+SQRT(T**2 + LVAL(C)*(LVAL(C)+1d0)))/K(IC,IA)
          if (rturnc.gt.rmorto) rmorto=rturnc
        ELSE
          Rmorto = rmatch
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if(CDETR.gt.0) 
     +  write(KO,*) 'L-cutoff at ',CUTVAL(C),ICUTC,' for C,Lin,J=',
     +       		C,LVAL(C),real(JTOTAL)
      endif
345       continue
       if (hort.gt.0.) then
         write(*,'(" *** Solutions stabilized at intervals of",f5.1,
     +   " fm, and up to",f7.1," fm")') hort,rmorto
       endif
       CUTOFF = N
      DO 346 C=1,NCH
       CUTVAL(C) = max(2,CUTVAL(C))
346       CUTOFF = min(CUTOFF,CUTVAL(C))
       if(CDETR.gt.0) write(KO,*) 'min-cutoff = ',CUTOFF
c
350	continue

      HERE = .true.
      if(NCH==0.or.NOSOL) go to 720
      if(MINTL==0.and.melfil==0) go to 720

	NCLREQ = maxval(NCLIST(1:NCH,1:NCH))
C
	if(.not.rterms) then

	if(rela /= '  ') then
	  DO 85 C1=1,NCH
	  DO 85 C2=1,NCH
	  DO 85 NC=1,NCLIST(C1,C2)
85	  CLIST(C1,C2,NC) = GAM(PART(C1,1),EXCIT(C1,1))*CLIST(C1,C2,NC)
	endif


C    SOLVING THE SET OF COUPLED CHANNELS FOR EACH INCOMING CHANNEL
C    -------------------------------------------------------------
      REPEAT = .FALSE.
	SKIP(:,1) = .false.
      DO 700 JIN=1,MINTL
         EL = INITL(JIN)
         JPSET = JPSET + 1
         IF(JSET.GT.0 .AND. JPSET.GT.JSET) GO TO 700
         DO 360 C=1,NCH
          SOMA(C,1) = 0.0
          SIMPLE(C) = .TRUE.
          SKIP(C,2) = .FALSE.
C                      SKIP(C,1) = SOURCE SAME AS LAST TIME
C                      SKIP(C,2) = FED(C) = SOURCE TERM NON ZERO
	 SRC(1:N,C) = 0.0
         IF(DISCIT) write(8,rec=NCH+C) SRC(1:N,C)
         IF(DISCIT) written(8) = .true.
        IF(C<=NICH) then
	 PSI(1:N,C) = 0.0
         IF(DISC8)write(8,rec=C) PSI(1:N,C)
	endif
  360    CONTINUE
	!write(6,*) 'READWF alloc ',allocated(READWF); call flush(6)
	READWF(:) = .false.
	ELCOEF(:) = 0.
        DO NC=1,NCLIST(EL,EL)
	  JF = NFLIST(EL,EL,NC)
	 ELCOEF(JF) = ELCOEF(JF) + CLIST(EL,EL,NC)
	 enddo
C
      FAIL = .FALSE.
      AL = 0.0
      BEST = 1E5
      NBEST = 0
       NSOLI = NSOL
       if(FJSWTCH) NSOLI = 1
       if(FJSWTCH) NBEST = 1
      DO 520 ITNL=1,NSOLI
       call flush(6)
      IF(.NOT.REPEAT) THEN
C
         DO 438 C=1,NCH
          IT = ITC(PART(C,1),EXCIT(C,1))
         L = LVAL(C)
           IF(ISOCEN.eq.1) L = JTOTAL + 0.1
           IF(ISOCEN.eq.2) L = LVAL(EL)
           LL1(C) = L*(L+1d0)
         DO 437 I=1,2
  437    CH(I,C) = CHL(L+1,IT,I) * CI * HALF
         IF(CDETR.GE.2) WRITE(KO,1300) C,(CH(I,C),I=1,2)
 1300      FORMAT('0',I6, ' MATCH TO',1P,2E13.4, ' AND',2E13.4)
  438    CONTINUE
      ENDIF
C
      WREQ = ITNL.LT.NSOLI .OR. WOUT
      LCERW1 = SMATL.GE.5.AND.PADE.EQ.0.OR.SMATS.GE.6.OR.SMATL.GE.3.AND.
     X ITNL.EQ.NSOLI
      LCERW2 = SMATS.GE.2.AND.ITNL.EQ.1
	RERR=ACC8
C
      ! if(IEX.lt.NCH.and.FALLOC.and.IBLOCK.lt.99) then
      ! if(IBLOCK.lt.99) then
      ! if(IEX.lt.NCH.or.FALLOC) then
      ! if(IEX.lt.NCH.or.WOUT.or.ITER>0) then
      if(FALLOC) then
	
      CALL ERWIN(PSI,ECM(1,1),ECM(1,2),IEX,FORMF,NF,FORMF,SRC,CH,
     X EL,SMAT,LVAL,JVAL,JTARG,LL1,NCH,NICH,N,ECM(1,3),M,MD,TSTD,RENORM,
     X REPEAT,WREQ,BLOCKD,AL,RERR,ICUTC,CDETR ,LCERW1, LCERW2,
     X CUTOFF,SKIP(1,1),SKIP(1,2),PTYPE,ferw,FIM,FIMD,LOCFIL,
     X CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)

      else

c		Erwin version optimised for pure CC S-matrix solutions:
      CALL ERWINCC(ECM(1,1),ECM(1,2),FORMF,NF,FORMF,CH,REPEAT,
     X EL,SMAT,LVAL,JVAL,JTARG,LL1,NCH,N,ECM(1,3),M,MD,TSTD,RENORM,
     X AL,RERR,ICUTC,CDETR ,LCERW1, LCERW2,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     X CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)

      endif
C            READ IN INITIAL WFS, if requested
      IF(INITWF.NE.0.and..not.INITWFE) THEN
	  NRF = abs(INITWF)
	  PSI(:,:) = 0d0     ! all wfs=0 EXCEPT those read in!
!          call openif(NRF)
          NA = (N-1)/2*2 + 1 ! read these
      if(ITNL==1) then ! first iteration
      IF(INITWF>0) READ(NRF,657,END=1692,ERR=1691) NAI,HPI,ENLABI,
     x                                    JTOTALI,PARITYI
      IF(INITWF<0) READ(NRF,END=1692) NAI,HPI,ENLABI,JTOTALI,PARITYI
	if(NA/=NAI) then
 	  write(0,*) 'STOP: read wfs NA ',NA,' neq ',NAI; stop
	 endif
	if(abs(JTOTAL-JTOTALI)>0.1.or.PARITY/=PARITYI) then
 	  write(0,*) 'STOP: read J/pi ',JTOTALI,PARITYI
	  stop ' J/pi in wrong sequence'
	 endif
!657   FORMAT(I4,2F8.4,F8.1,I3,2f12.6,2f8.3)
      DO 1690 C2=1,NCH+1     ! look for at most NCH+1 channel wfs !
        IF(INITWF.GT.0)
     X     READ(NRF,660) ITI,LVALI,JVALI,JTOTALI,LVALE,JVALE,SMATI
        IF(INITWF.LT.0)
     X     READ(NRF) ITI,LVALI,JVALI,JTOTALI,LVALE,JVALE,SMATI
	  if(ITI<0) go to 1695
            IF(LVAL(EL).NE.LVALE .or. abs(JVAL(EL)-JVALE)>0.1) then
 	    write(0,*) 'STOP: wrong elastic I,J: ',LVALE,JVALE
	    stop  ' Elastic channels in wrong sequence'
	    endif
          DO 1670 C=1,NICH  !  look for this channel in current set
            IF(ITC(PART(C,1),EXCIT(C,1)).NE.ITI) GOTO 1670
            IF(LVAL(C).NE.LVALI.or.abs(JVAL(C)-JVALI)>0.1) GOTO 1670
	    write(ko,1696) 1,C,SMATI
	     SMAT(C) = SMATI
	     READWF(C) = .true.
           IF(INITWF>0) then
             READ(NRF,'(6e12.4)') (PSI(I,C),I=1,NA)
!	     write(6,*) ' Last 3 wf values ='
!             write(6,'(6e12.4)') (PSI(I,C),I=NA-2,NA)
	    ELSE
             DO 1665 I=2,NA
 1665        READ(NRF) PSI(I,C)
            ENDIF
 1670      CONTINUE
! 681      FORMAT(1P,6E13.6)
 1690  CONTINUE
       go to 1695
 1691  backspace NRF
	write(6,*) 
	write(6,*) ' Read error on file', NRF,' on line:'
	read(NRF,'(A)') LINE
	write(6,'(/1x,a)') line
	 stop
       
 1692	INITWFE = .true. ! reached EOF
 1695  CONTINUE
       else !  ITNL > 1, on later iterations, just read from file 8
	!write(6,*) 'READWF alloc2',allocated(READWF); call flush(6)
	 do C=1,NICH
         if(READWF(c)) then 
 	   read(8,rec=C) PSI(1:N,C)
  	   SMAT(C) = SOMA(C,1) 
	    write(ko,1696) ITNL,C,SMAT(C)
 1696    format(' Iteration #',i2,': fixed wf for  channel ',i4,
     x          ', so S =',2f12.8)
	   endif
	 enddo  ! C
       endif    ! ITNL==1
       ENDIF    ! INITWF/=0

      IF(NPOST.AND.ITNL.LT.NSOLI.AND.ITNL.GT.1)
     X    CALL RENO(PSI,PHI,N,ECM(1,3),NCH,CUTOFF,NLN,NLO,MR,MLM,0,
     X         NF,  0,NL,EMPTY,MLT,NLC,SMATL,CLIST,NCLIST,NFLIST,CHNO,
     X              NPOST,NPRIOR,SRC,ECM(1,2),SIMPLE,
     X              FORMF,FORMF,LOCFIL,ECM(1,1),LVAL)
C
      IF(FLIP) CALL XCH(PSI,N,NCH,NICH,SMAT,JTOTAL,LVAL,
     X      JVAL,JPROJ,JTARG,PART,EXCIT,COPY,MXP,MXX,MAXN,SMATL,SAME,
     X      .false.,EXCH,MAXCH,EL)
C
C      PADE ACCELERATION
C      -----------------
        DO 470 C=1,NCH
470     SRATIO(C) = SMAT(C)
      IF(PADE.GE.1)    CALL PAD(PADE,SMAT,NCH,ITNL,SPAD,MAXCH,MAXIT+2)
        DO 480 C=1,NCH
        SRATIO(C) = SMAT(C) / (SRATIO(C) + (1E-20,0.0))
 480    IF(.NOT.PSIREN) SRATIO(C) = (1.0,0.0)
C
C               THAT IS THE PADE APPROXIMANT BY THE EPSILON ALGORITHM
C               -----------------------------------------------------
      FRACT = 0.0
      DO 500 C=1,NCH
         T = ABS(SMAT(C))
         EMPTY(C) = T.LE.1E-12
         SAME(C) = .TRUE.
         DIFF = ABS(SMAT(C) - SOMA(C,1)) * 100.
     X            * (2.*JTOTAL+5.)/(2.*JTMAX+5.)
         SAME(C) = DIFF.LE.1E-6
         IF(ECM(C,1).LT.0.) SMAT(C)=0.  ! closed channels give no xs, but can still vary
         IF(WREQ.AND.DISC8.and.C<=NICH)  write(8,rec=C) PSI(1:N,C)
         IF(FRACT.GT.DIFF.OR.EXTRA(2,PART(C,1),EXCIT(C,1))) GO TO 490
            FRACT = DIFF
            C2 = C
  490    SOMA(C,1) = SMAT(C)
  500  CONTINUE
       IF(FRACT.LT.BEST.AND.ITNL-1.GE.ITMIN.and.
     X		(MOD(ITNL,2)==1.or.FRACT<abs(IPS))) THEN
          BEST = FRACT
          NBEST = ITNL
C   Store best SMAT in SOMA(*,2) and wfns in file 18
          IF(WOUT) THEN
           IF(.NOT.OPEN18) THEN
            IF(MACH==6.or.MACH==7.or.MACH==8.or.MACH==3.or.MACH==4)THEN
               OPEN(18,ACCESS='SEQUENTIAL',
     X            FORM='UNFORMATTED',FILE=trim(TMPN)//'18')
            ELSE
               OPEN(18,ACCESS='SEQUENTIAL',
     X            FORM='UNFORMATTED',STATUS='SCRATCH')
            ENDIF
           ENDIF
           OPEN18 = .TRUE.
           REWIND 18
          ENDIF
          DO 501 C=1,NCH
          SOMA(C,2) = SMAT(C)
          IF(WOUT.and.C<=NICH) WRITE(18) (PSI(I,C),I=1,N)
          IF(WOUT            ) WRITE(18) (SRC(I,C),I=1,N)
  501     continue
	  if(WOUT) written(18) = .true.
       ENDIF
      DO 505 I=1,N
  505 IF(ITNL.eq.1) VPOT(I,3) = PSI(I,EL) / SRATIO(EL)
      IF(PSIREN.AND.WOUT.AND.ABS(SRATIO(EL)-1.).GT.0.1.AND.ITNL.GE.3)
     X  WRITE(KO,1395) SRATIO(EL),FRACT
 1395 FORMAT(' Psi(EL) renorm. by ',1P,2E9.2,', Conv =',0P,F9.5,' %')
       IF(NSOLI.GT.4.AND.ITNL-1.GE.1.AND.
     X     (SMATL.GE.3.OR.SMATL.GE.2.AND.IT0.GT.0)) THEN

          FUSIO = 0.0
          DO C=1,NCH
             IC = PART(C,1)
             IA = EXCIT(C,1)
             IT = ITC(IC,IA)
          IF(ECMC(IC,IA) .GT. 0.0) THEN
          AMDSQS = ABS(SMAT(C))**2
          IF(C.EQ.EL) AMDSQS = 1 - AMDSQS
!@@
             IF(RMASS(IC).lt.1e-5) then
                XSC = JCOEF* AMDSQS * 
     X                ( RMASS(PEL)*amu/ (HBC*K(PEL,EXL)) )
             ELSE IF(RMASS(IC).lt.1e-5) then
                XSC = JCOEF* AMDSQS / 
     X                ( RMASS(IC)*amu/ (HBC*K(IC,IA)) )
             ELSE
                XSC = JCOEF* AMDSQS *
     X                K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
             ENDIF
!            see justification on frxx3.f
!
!            XSC = JCOEF* AMDSQS *
!    X             K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
!@@
          IF(C.EQ.EL) THEN
             FUSIO = FUSIO + XSC
            ELSE
             FUSIO = FUSIO - XSC
            ENDIF
          endif
          ENDDO
              IF(PADE.GT.0) WRITE(KO,1400) ITNL-1,FRACT,C2,SMAT(C2),
     X                           FUSIO,SPAD(C2,MIN(ITNL,MAXIT+2),2)
              IF(PADE.EQ.0) WRITE(KO,1400) ITNL-1,FRACT,C2,SMAT(C2),
     X				 FUSIO
        ENDIF
 1400 FORMAT(' So max change at iter. #',I3,' =',F10.4,' %, @ CH',I3,
     X ' with',2F10.5,', Fus=',F10.5,:,' prev S=',1p,2e12.4)
      IF(SMATS.GE.5.AND.PADE.GE.1) WRITE(KO,1430) (SMAT(C),C=1,NCH)
      IF(ITNL-1.GE.ITMIN.AND.FRACT.LT.ABS(IPS) .OR. DRY) GO TO 530
      WERR = 3.* AL*RERR
      IF(FRACT.LT.WERR*0.1.AND.ITNL-1.GT.ITMIN.AND.IPS.GT.0.0) THEN
        WRITE(KO,1401) ITNL-1,FRACT,WERR
1401      FORMAT('0STOPPING AFTER',I4,' ITERATIONS, as change',F9.3,
     X         ' % is less than 0.1 of ERWIN accuracy loss ',F9.3,' %'/)
         GO TO 530
         ENDIF
      IF(FRACT .GT. 100.*BEST .AND. NBEST-1.GT.ITMIN.AND.IPS.GT.0.) THEN
C               Give up, as appears to be diverging
         WRITE(KO,1402) ITNL-1,FRACT,BEST,NBEST-1
1402     FORMAT('0STOPPING AFTER',I4,' ITERATIONS, as change of ',F9.3,
     x ' % is worse than 100 times best change of',F9.3,'% at IT =',I3/)
         GO TO 530
         ENDIF
      IF(ITNL.EQ.NSOLI) GO TO 520
C
      CALL SOURCE(SRC,PSI,N,ECM(1,3),NCH,IEX,FORMF,NF,FORMF,CUTOFF,
     X  ICUTC,SIMPLE,.TRUE.,NLN,NLO,MR,NL,EMPTY,SAME,SHFL,
     X  WAVES,LVAL, MLT,SKIP,SKIP(1,2),NLC,CHNO,
     X  MLM,NICH,MMXICH,PTYPE,LOCFIL,CLIST,NCLIST,NFLIST,
     X  GAM,mxp,mxx,EXCIT,PART)! MGR  & AMoro 
C
      IF(NPRIOR.AND.ITNL.GT.1)
     X CALL RENO(SRC,PHI,N,ECM(1,3),NCH,CUTOFF,NLN,NLO,MR,MLM,1,
     X      NF,  NCH,NL,EMPTY,MLT,NLC,SMATL,CLIST,NCLIST,NFLIST,CHNO,
     X    NPOST.AND.ITNL+1.LT.NSOLI,NPRIOR,SRC,ECM(1,2),SAME,
     X        FORMF,FORMF,LOCFIL,ECM(1,1),LVAL)
C
      DO 510 C=1,NCH
  510  IF(DISCIT) write(8,rec=NCH+C) SRC(1:N,C)
  520 REPEAT = .TRUE.
      ITNL = NSOLI
C
c      IF(IPS.NE.0..AND.ITER.GT.4) THEN
      IF(IPS.NE.0..AND.NSOLI.GT.4) THEN
        WRITE(KO,1410) ITER
 1410   FORMAT(' FAILED TO CONVERGE AFTER',I3,' ITERATIONS !!!!!!!',/)
        IF(FATAL) WRITE(KO,*) ' Set ITER negative to allow continuation'
        FAIL = .TRUE.
        penalty = penalty + fine
       if(number_calls>5) then
	 write(6,41) number_calls,JTOTAL,PSIGN(parity+2),penalty
41	format('  At call ',i5,', cc set ',f6.1,a1,
     & 		' failed to converge, so  penalty =',1p,e10.1)
	else
C       IF(FATAL) STOP
        IF(FATAL) CALL ABEND(2)
       endif
      ENDIF
C
C    NOW HAVE CONVERGED RESULT in itnl iterations
C    --------------------------------------------
C 530 IF(ITNL.NE.NBEST.AND.NBEST.NE.0.OR.FAIL) THEN
  530 IF(IPS.NE.0.0.AND.ITNL.NE.NBEST.AND.ITER.GT.4.OR.FAIL) THEN
C               Find S-matrix elements and wfns from ITNL = NBEST
C               Get SMAT from SOMA(*,2), wfns etc. from file 18.
          FRACT = BEST
          IF(WOUT) REWIND 18
          DO 531 C=1,NCH
          SMAT(C) = SOMA(C,2)
          IF(WOUT.and.C<=NICH) READ(18) (PSI(I,C),I=1,N)
  531     IF(WOUT) READ(18) (SRC(I,C),I=1,N)
        ENDIF
C
      IF(PADE.GE.1.OR.ITER.GT.4) THEN
      ITNL = ITNL - 1
      IF(SMATL.GE.2) WRITE(KO,1420) JTOTAL,PSIGN(PARITY+2),
     X                          NBEST-1,ITNL,FRACT,JTMAX,SMAT(EL)
      IF(SMATL.GE.3) WRITE(KO,1430) (SMAT(C),C=1,NCH)
 1420 FORMAT(' Final S-matrices (',F7.1,A1,') at',I3,' after',I3,
     X       ' accurate to',
     X F9.4,' % of J=',F7.1,' unitarity:',F10.5,' +i*',F8.5,',')
 1422 FORMAT(F10.1,2F12.8,'i: elastic S-matrix  @@',f10.2,i3,i4,l2,i3)
 1430 FORMAT(5(1X,F11.5,' +i*',F9.5,',') )
!1431 FORMAT(5(1X,1p,e11.3,' +i*',e9.1,',') )
      ENDIF
      MMXIT = MAX(MMXIT, ITNL)
C
	if(WOUT) then
        DO 550 C=1,NICH
        DO 550 I=1,N
        PSI(I,C) = PSI(I,C) * SRATIO(C)
550     SRC(I,C) = SRC(I,C) * SRATIO(c)
	endif
C
      IF(SMATL.GE.2) WRITE(KO,1422) JTOTAL,SMAT(EL),TM(I)-TIME0J,
     X			DROPPED,NSTEPD,FJSWTCH,IAME
	phasin(JIN) = log(SMAT(EL))*(0.,-0.5)*180./pi
	linel(JIN) = EL

      NLN1 = NLN
	include 'usescatwf.f'
C
!      IF(CDETR.GT.1) CDETR = CDETR - 1
C                next JIN    :
  700 CONTINUE

        else

C   SOLVE SET OF COUPLED EQUATIONS BY R-MATRICES
C   --------------------------------------------
C
        IF(JSET.GT.0 .AND. JPSET.GT.JSET) GO TO 720

      IF(FLIP) then
        CALL XCH(PSI,N,NCH,NICH,SMAT,JTOTAL,LVAL,
     X      JVAL,JPROJ,JTARG,PART,EXCIT,COPY,MXP,MXX,MAXN,SMATL,SAME,
     X      .true.,EXCH,MAXCH,EL)
       else
         EXCH(:,:) = 0.0
       endif
        nbas = nrbases
        if(FJSWTCH) nbas = 0
        CALL RMATRIX(JTOTAL,NCH,MINTL,INITL,ECM,ECM(1,2),ECM(1,3),ITC,
     X    CHL,PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,MRM,FORMF,NF,ITCM,EOFF,
     X    CUTVAL,ISOCEN,ebeta(IPARIT-1),weak,nbas,nrbmin,nbas*NCH,
     X    KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH,LL1, NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigs,
     X    MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,ETA,
     X    FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,IAME,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X    VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,
     x    linel,IEXCH,FUSJPI)
        JPSET = JPSET + MINTL
        endif
C
	     if(ETOTAL.lt.0) go to 720    
	     if(abs(NLAB(1))>1.or.num_energies>0) then

	     do JIN=1,MINTL
	     I = nphases+JIN
	     if(I<=20) then
		phase(I) = phasin(JIN)
		if(LEN==1.and.final) then
		if(abs(JTOTAL)<100.) then
		write(71,88) I-1,
     x  	 JTOTAL,PSIGN(PARITY+2),mod(LVAL(linel(JIN)),100),JIN
88	    format('@legend string',I3,' "',f4.1,a1,'/',i2,i2,'"')
		else
		write(71,881) I-1,
     x  	 JTOTAL,PSIGN(PARITY+2),mod(LVAL(linel(JIN)),100),JIN
881	    format('@legend string',I3,' "',f8.1,a1,'/',i2,i2,'"')
		endif
		else
		 if(phase(I)<phaslst(I)-90) phaseadd(I)=phaseadd(I)+180
		 if(phase(I)>phaslst(I)+90) phaseadd(I)=phaseadd(I)-180
		endif
		phaslst(I)=phase(I)
		phase(I) = phase(I) + phaseadd(I)
	     endif ! I<=20
	     enddo

		do id=1,datasets
		if(data_type(id)==4
     x	          .and.abs(JTOTAL-data_jtot(id))<.1 
     x            .and.PARITY==data_par(id)) then
		 do ip=1,datalen(id)
 		  if(abs(enlab-data_energies(ip,id))<1e-5) then
		   JIN = data_ch(id)
		   JIN = min(JIN,MINTL); JIN=max(JIN,1)
		   T = phase(nphases+JIN)
	 	   EE = (T-datavals(ip,id))/dataerr(ip,id)
                   data_chisq(id) = data_chisq(id) + EE**2
                   theoryvals(ip,id) = T
!                  theoryplot(ip,id) = T
		   endif
	          enddo
		  endif
	         enddo

	     nphases = nphases+MINTL
	     endif

      WRITE(156,1447) JTOTAL,FUSJ, PARITY
      DONE = DONE + DONES
      T = TM(I)
*     write(48,*) 'DONE incremented by ',DONES,' to ',DONE,real(T)
      I = (JTOTAL-JAP)/JUMP(JBLOCK,1) + 0.5
      IF(DONE.GE.MAX(3,MINTL+1).and.I.gt.4) DONE=DONE+1000

	NJDONE = NJ
      IF(DONE.GE.1000) GO TO 750
C                next PARITY :
  720 CONTINUE
      if(IAME==0.and.HERE.and.SMATL>0) then
      if(NFUS==0) then
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ
	else if(NFUS==1) then
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ,CFUSJ(1)
	else
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ,
     x		CFUSJ(1:min(7,NFUS)),sum(CFUSJ(1:NFUS))
	endif
	call flush(56)
       endif
 1446 	FORMAT(f8.1,11G12.4)
 1447 	FORMAT(f8.1,G12.4,I3)
       IF(LISTCC.GT.0.and.LISTCC.le.90) LISTCC = LISTCC - 1
      IF(DRY) LISTCC = 0
      IF(DONE.GE.1000)GO TO 750
C                next JTOTAL :
  730 CONTINUE
C                next JBLOCK :
  740 CONTINUE
	  !   if(ECMC(PEL,EXL).lt.0) go to 948
          T = TM(I)
  750 continue
      REWIND 10
      T = TM(I)
        if(OPEN18) then
		close(18,status='delete')
		OPEN18 = .false.
		written(18) = .false.
	  endif
!       if(OPEN12) then
!	close(12,status='delete')
!	OPEN12 = .false.
!  	endif
        T = TM(I)
      if(final) write(KOI,*) 'Finished all CC sets @ ',real(T)
      call flush(KOI)
C  Serial operation for this section:
         KS=10
         REWIND KS
         allocate (FUSUM(MAXCH,NSA*NJA))
         JN=-1.0
        DO 761 JF=1,JCCSET
!        write(48,*)' reading ',JF,'/',JCCSET,' from file 10',real(T)-TM0
         READ(KS,END=761) JTOTAL,NCH,MINTL,IPARIT,JUMPER
         READ(KS) (LVAL(C),JVAL(C),PART(C,1),EXCIT(C,1),C=1,NCH)
	   DO C=1,NCH
      		JPROJ(C)= JEX(1,PART(C,1),EXCIT(C,1))
      		JTARG(C)= JEX(2,PART(C,1),EXCIT(C,1))
	   ENDDO
          RS = 0.0
           IF(abs(JTOTAL-JN).gt.0.1) THEN
              RS=1.0
           ENDIF
           DO 751 C=1,NCH
           DO 751 I=1,NSA*NJA
  751      FUSUM(C,I) = 0.0
         DO 753 JIN=1,MINTL
          READ(KS,end=753) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1)
C
C     FIND CUMULATIVE SUMS OF REACTION CROSS SECTIONS
C     AND GIVE SUMMARY OF S-MATRIX ELEMENTS ON FILE 7
C     -----------------------------------------------
         DO 752 IT=0,ITCM
  752    SIGJ(IT) = 0.0
      JCOEF = (2*JTOTAL+1)* PI  * XCOEF / ABS(K(PEL,EXL))**2
      CALL SUMX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGT,SIGJ,SIGR,PART,EXCIT,
     X           JUMPER,FUSL,CSIG,JEX,ITC,K,RMASS,PEL,EXL,LVAL,
     X         NSA,NJA,JVAL,JPROJ,JTARG,JTOTAL,TOTFUS,FUSUM,CORFUS,
     X         SIGEL,SIGTOT,SMATS,ITCM,XSTABL,IF1,IF2)
  753    CONTINUE
!      write(48,*) ' read ',JF,' = ',JTOTAL,' from file 10',real(T)-TM0
      DO 755 I=1,NSA
      SIGFUS(I,2) = 0.0
  755 SIGFUS(I,3) = 0.0
      T = JCOEF*NSA*NJA*JUMPER
      IAM = 0
      DO 758  I=1,NSA
      T4 = 0.0
         DO 756 J=1,NJA
         IAM = IAM + 1
          DO 756 C=1,NCH
  756    T4 = T4 + ABS(FUSUM(C,IAM))**2
  758 SIGFUS(I,IPARIT) = T4 / NJA
      DO 759 I=1,NSA
        AMSA = I-1 -JEX(1,PEL,EXL)
        T4 = 0.0D0
        IF(JTOTAL+0.2.GT.ABS(AMSA)) 
     X   T4 = T * (RS - SIGFUS(I,2) - SIGFUS(I,3))
  759   SIGFUS(I,1) = SIGFUS(I,1) + T4
C
       JN=JTOTAL
  761  CONTINUE
       deallocate (FUSUM)
       if(MAXCH>0.and.number_calls<0.and..false.) then ! dealloc if not needed more
       if(final) write(KO,*) 'Deallocate channel arrays'
	!write(6,*) 'READWF dealloc ',allocated(READWF); call flush(6)
         deallocate (CLIST,NCLIST,NFLIST,CH,ECM,XS,INCOME,
     X    CFMAT,CGMAT,CRCRAT,CUTVAL,JVAL,JPROJ,JTARG,CHNO,CFF,
     x   LVAL,PART,EXCIT,INITL,BLOCKD,LL1,SMAT,EXCH,PSI,SAME,READWF)
         if(.not.rterms) deallocate(SRC,SRATIO,SOMA,SPAD,
     X         FIM,FIMD,SIMPLE,EMPTY,SKIP)
         if(.not.rterms) deallocate(ferw)
       
       MAXCH=0
       MAXICH=0
      endif
C
C    GIVE ACCUMULATED AND DIFFERENTIAL CROSS SECTIONS
C    ------------------------------------------------
	if(num_energies>3) then
	MCCSET = JCCSET
!	if(LEN==1) allocate (alphas(MCCSET,MAXXCH,nrbmax*mag))

	endif
      DO 764 CP=1,NCP
  764 IF(MLCALL(2,CP)) CALL NLSTAT(CP,NLOC(1,CP),NLM,EPC*.01,HNL,CENTRE)
C
      DO 765 I=2,3
      SIGEL(I) = SIGEL(I)/(SIGEL(1) + 1E-20)
      SIGTOT(I) = SIGTOT(I)/(SIGTOT(1) + 1E-20)
      SIGT(I) = SIGT(I)/(SIGT(1) + 1E-20)
      TOTFUS(I) = TOTFUS(I)/(TOTFUS(1) + 1E-20)
	DO 765 IA=1,NFUS
  765 CORFUS(I,IA) = CORFUS(I,IA)/(CORFUS(1,IA) + 1E-20)
      if(.not.final) then
      if(SIGT(1)>1e-5) WRITE(KO,1450) SIGT(1),SIGTOT(1),SIGEL(1)
 1450 FORMAT('  Sig-R,T,E =',3F11.5)
      DO 770 IC=1,NCHAN
      T = maxval(abs(SIGR(IC,1:NEX(IC))))
      if(T>1e-20) then
        WRITE(KO,14611) IC,(SIGR(IC,IA),IA=1,NEX(IC))
	endif
  770 continue
      else
      WRITE(KO,1451) 'REACTION',SIGT
 1451 FORMAT('0CUMULATIVE ',a8,' cross section                 =',
     X  F11.5,:,'  <L> =',F9.2,'  <L**2> =',F9.1)
      if(ABS(ETA(PEL,EXL))<1e-9) then
         WRITE(KO,1451) ' TOTAL  ',SIGTOT
         WRITE(KO,1451) ' ELASTIC',SIGEL
	 endif
      DO 771 IC=1,NCHAN
      T = maxval(abs(SIGR(IC,1:NEX(IC))))
      if(T>1e-3.or.T.eq.0d0) then 
        WRITE(KO,1460) IC,(SIGR(IC,IA),IA=1,NEX(IC))
	else
        WRITE(KO,1461) IC,(SIGR(IC,IA),IA=1,NEX(IC))
	endif
  771 continue
      endif
	call flush(KO)
 1460 FORMAT('0CUMULATIVE outgoing cross sections in partition',I2,
     X  ' :',7F11.5,/,(52X,7F11.5))
 1461 FORMAT('0CUMULATIVE outgoing cross sections in partition',I2,
     X  ' :',1p,7e11.3,/,(52X,7e11.3))
14611 FORMAT('  Out in',I2,' :',1p,7e11.3,/,(14X,7e11.3))
      if(final) then
      write(13,*) 'Integrated cross sections for all states'
      write(13,'(3i3,1p,e12.4,0p,f8.4)') PEL,EXL,NCHAN,ENLAB,BEPROJ
      do IC=1,NCHAN
      write(13,'(2i4)') IC,NEX(IC)
        do IA=1,NEX(IC)
         IT = ITC(IC,IA)
	if(IC==PEL.and.IA==EXL) then
          write(13,'(2(f5.1,i3,f8.4),1p,3e12.4)') 
     X    (JEX(j,IC,IA),BAND(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGT(1),
     X    SIGTOT(1),SIGEL(1)
        else
	  if(abs(ENEX(1,IC,IA))+abs(ENEX(2,IC,IA))<.01.and.CDCC/=0) then
          write(13,'(2(f5.1,i3,f8.4),1p,e12.4)')   ! special adiabatic case
     X    (JEX(j,IC,IA),BAND(j,IC,IA),ENEX0(j,IC,IA),j=1,2),SIGR(IC,IA)
	  else
          write(13,'(2(f5.1,i3,f8.4),1p,e12.4)') 
     X    (JEX(j,IC,IA),BAND(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGR(IC,IA)
          endif
	endif
	call flush(13)
	enddo
      enddo
      written(13) = .true.
      WRITE(40,14651) ENLAB,TOTFUS(1),SIGT(1),
     X               ((SIGR(IC,IA),IA=1,NEX(IC)),IC=1,NCHAN),
     x               SIGTOT(1),SIGEL(1)
      written(40) = .true.
      call flush(40)
      endif
      if(final) WRITE(KO,1465) TOTFUS,NAME(1,PEL),(SIGFUS(I,1),I=1,NSA)
 1465 FORMAT('0Cumulative ABSORBTION by Imaginary Potentials     =',
     X F11.5,'  <L> =',F9.2,'  <L**2> =',F9.1,
     X/,'   Fusion for specific ',A8,' M-states :',4F12.6,/,1X,
     * 8F12.6,/,6F12.6)
        WRITE(KO,1451) 'OUTGOING',SIGT(1)-TOTFUS(1)
	call flush(KO)
!14651 FORMAT(F8.4,16G12.4/(8X,16G12.4))
14651 FORMAT(G12.4,16G12.4/(12X,16G12.4))
!  Put total fusion & reaction cross section in file 39:
!!      WRITE(39,14651) ETOTAL,TOTFUS(1),SIGT(1)
!  Fusion & reaction cross section in file 39:
!lab      IF(NFUS.eq.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1)
!lab      IF(NFUS.NE.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1),CORFUS(1,:)
      if(final) then
      IF(NFUS.eq.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1),SIGTOT(1),
!     IF(NFUS.eq.0) WRITE(39,14651) ETOTAL,TOTFUS(1),SIGT(1),SIGTOT(1),   ! not for fission work Nov 09
     x                   SIGT(1)-TOTFUS(1),SIGEL(1)
      IF(NFUS.NE.0) WRITE(39,14651) ETOTAL,TOTFUS(1),SIGT(1),SIGTOT(1),
     x                   SIGT(1)-TOTFUS(1),SIGEL(1),CORFUS(1,:)
      written(39) = .true.
      call flush(39)
      endif
      do IA=1,NFUS
      WRITE(KO,1466) IA,KFUS,CORFUS(:,IA)
1466  FORMAT('0Cumulative ABSORBTION in state ',i3,' by Imaginary ',
     X  'Potential',I3,'   =',F11.5,'  <L> =',F9.2,'  <L**2> =',F9.1)
      enddo
      sfa = 1.
      if(abs(ETA(PEL,EXL))<50.) then
      sfa = exp(2d0*PI*ETA(PEL,EXL)) *
     x    ENLAB * MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
      if(final) WRITE(KO,14661) sfa
14661 FORMAT(' To convert to S-factors (MeV.mb = keV.b), ',
     X  ' multiply by',1P,E12.4/)
      I=0
       DO 14662 IC=1,NCHAN
       DO 14662 IA=1,NEX(IC)
       if(.not.(IC==PEL.and.IA==EXL)) then
	 I=I+1
	 SFAC(I) = SIGR(IC,IA)*sfa
         endif
14662 	continue
      if(I>0.and.final) then
	   call openif(35); call openif(75)
         WRITE(35,14651) ETOTAL,SFAC(1:I)
         WRITE(75,14651) ENLAB,SFAC(1:I)
	 call flush(35); call flush(75)
         written(35) = I>0
         written(75) = I>0
       endif
      endif
		do id=1,datasets
		if(data_type(id)==3) then

!               Adjust any datanorm search parameter!
           	datanorm=1.0
           	do ip=1,nvars
           	if(srch_kind(ip)==5.and.srch_datanorm(ip)==id) 
     x       		datanorm = datanorm * srch_value(ip)
           	enddo
!           WRITE(KO,*) ' dataset ',id,kq1,lq1,real(datanorm)

		 do ip=1,datalen(id)
 		  if(abs(enlab-data_energies(ip,id))<1e-5) then
	           if(data_ic(id)>0) then
		    T = sigr(data_ic(id),data_ia(id))
		   else
		    if(data_ia(id)==0) T = SIGT(1)
		    if(data_ia(id)==1) T = TOTFUS(1)
		    if(data_ia(id)>=2) T = CORFUS(1,data_ia(id)-1)
                    if(data_ia(id)==-1) T = SIGEL(1)   ! neutron elastic
                    if(data_ia(id)==-2) T = SIGT(1)+SIGEL(1)  ! neutron total=elastic+reaction
	 	   endif
		   if(data_idir(id)==-1) then  ! convert to absolute, first time
		       datavals(ip,id) = datavals(ip,id)/sfa  
		       dataerr(ip,id) = dataerr(ip,id)/sfa  
		     endif
	 	   EE = (T/datanorm-datavals(ip,id))/dataerr(ip,id)
                   data_chisq(id) = data_chisq(id) + EE**2
                   theoryvals(ip,id) = T
		   endif
	          enddo
		  endif
	         enddo
      SIGR(PEL,EXL) = SIGR(PEL,EXL) + SIGT(1)
	if(final.and.SMATL>0) write(56,*) '&'
	if(final.and.SMATL>0) write(156,*) '&'

C
  772 IF(abs(THINC).lt.1e-9) GO TO 785
	rewind 10
C
      NLJ = NINT(2. * LJMAX + 1.)
      LCROSS = 16
      LFAM=37
      LXSEC=0 ! 199  Not used now, but you can set to unused file no.

      JCCSET = 0
      LEG = 0
      DO 780 IC=1,NCHAN
      DO 780 IA=1,NEX(IC)
      IF(BAND(1,IC,IA) .EQ. 0) GO TO 775
C
      
      if(final) WRITE(KO,1470) (NAME(IN,IC),IN=1,2),IA,(JEX(IN,IC,IA),
     X      PSIGN(SIGN(1,BAND(IN,IC,IA))+2)    ,IN=1,2),CHSIGN(IA)
 1470 FORMAT(/' CROSS SECTIONS FOR OUTGOING ',A8,' & ',A8,' in state #',
     X I4,' with spins & parities',F4.1,1X,A2,' &',F5.1,1X,A2,';',i3/)
         IF(DRY .AND..NOT.(IC.EQ.PEL .AND. IA.EQ.EXL)) GO TO 775
         IF(.NOT.GIVEXS(IC)) GO TO 775
         IF(ECMC(IC,IA) .LE. 0.0) GO TO 775
         IF(BAND(1,IC,IA)*BAND(2,IC,IA).EQ.0) GO TO 775
         IP = SIGN(1,BAND(1,IC,IA))
         MAXPLM = max(1,nint(JEX(1,IC,IA)+JEX(2,IC,IA)+
     x                 JEX(1,PEL,EXL)+JEX(2,PEL,EXL)))
         NMULTIES=0
       if(DGAM>0) then
         IN = PP-1
         DNAME = NAME(IN,IC)
         do 1471 IB=1,IA-1
         IF(COPY(IN,IC,IB,1)>0) go to 1471
           LGAM=abs(JEX(IN,IC,IA)-JEX(IN,IC,IB))
           LGAM = max(LGAM,1)   !  no L=0 photons
           if((-1)**LGAM * SIGN(1,BAND(IN,IC,IA))*SIGN(1,BAND(IN,IC,IB))
     x        <0) LGAM=LGAM+1
           if(LGAM > JEX(IN,IC,IA)+JEX(IN,IC,IB)) go to 1471
             NMULTIES=NMULTIES+1
	     DSPINS(NMULTIES) = JEX(IN,IC,IB)
             DMULTIES(NMULTIES) = LGAM
             DLEVEL(NMULTIES) = IB
1471	 continue
        endif
      CALL CRISS(IC,IA,JEX(1,IC,IA),NSA,NJA,PEL,EXL,JEX(1,PEL,EXL),
     X           LUSED,K,ETA,RMASS,CSIG,KOORDS,MMXCH,MAXF,ENLAB,LEN,
     X           THMIN,THMAX,THINC,EXTRA(1,IC,IA),LAMPL,NEARFA,KQMAX+1,
     X           SIGR(IC,IA),XSTABL,LJMAX,NLJ,NJDONE,JBORD,JUMP,
     X           ITC(IC,IA),ITC(PEL,EXL),IEXCH,PP,PPKIND(PP),CDCC,IP,
     X           IBIN(IA),ENEX(1,IC,IA),PI,HEADNG,IOFAM(1,IC,IA),
     X		 IOFAM(2,IC,IA),LCROSS,LFAM,LXSEC,CHSIGN(IA),
     X           ETOTAL-EOFF,MASS(LIN,LAB),MASS(3-LIN,LAB),
     X           ECMC(IC,IA),MASS(LIN,IC),MASS(3-LIN,IC),LEG,MAXPLM,
     x           DSPINS,DMULTIES,NMULTIES,DLEVEL,DNAME)
	call flush(16)
C
775     continue

************************************************************************
  780 CONTINUE
          T = TM(I)
          if(final) write(KOI,*) 'Finished all xsecs @ ',real(T)
  785 KO = KOI

***      IF(VEFF.NE.0) CALL VEFFPOT
***      IF(BPM.NE.0) CALL BPMFUS
      IF(VEFF.NE.0) then
        include 'veffpot.f'
        close(15)
        endif
      IF(BPM.NE.0) then
        include 'bpmfus.f'
        endif


  948 continue
	deallocate (FNC,WNM)
C
C ----------------- NEXT INCIDENT ENERGY
       if(final.and.(abs(NLAB(1))>1.or.num_energies>0)) then
	 write(71,949) ENLAB,phase(1:min(20,nphases))
	 call flush(71)
	 endif
  949   format(f12.8,20f7.1)
  950 CONTINUE
C
  960 IF(IAME.gt.0) GO TO 970
       call flush(6)
      close(10,status='delete')
      written(10) = .false.
  970 if(final) call fkind(written,KO)
	close(40); close(71)
	if(OPN(1)) close(11)
	if(OPN(2)) close(9)
	if(OPEN12) then
	   close(12,status='delete')
	   written(12) = .false.
	   endif

      do id=1,datasets
       if(data_type(id)==6) then
        t = (srch_value(data_par(id))-datavals(1,id))/dataerr(1,id)
        data_chisq(id) = t**2
        endif
       if(data_idir(id)==-1) then
         data_idir(id) = 0
         endif
       if(data_type(id)==7.or.data_type(id)==8) then
        I = data_term(id)
         if(pralpha) write(6,*) 'Term',I,'E,W:',E_Brune(I),W_Brune(I)
         if(pralpha) write(6,*) 'E_Brune(;)',E_Brune
         if(data_type(id)==7)
     x       theoryvals(1,id) = E_Brune(data_term(id)) ! Brune energy
         if(data_type(id)==8)
     x       theoryvals(1,id) = W_Brune(data_term(id)) ! Brune total observed width
        t = (theoryvals(1,id)-datavals(1,id))/dataerr(1,id)
        data_chisq(id) = t**2
        endif
      enddo

      if(final.and.IAME==0) then
       WRITE(KO,1485) MAXQRN,MFNL,MLOC,LMAX1,MCLIST, 
     X        MMXQRN,NFNL,NF,LUSED+1,NCLREQ
 1485 FORMAT(/' PARAMETERS : MAXQRN   MFNL   MLOC  LMAX1 MCLIST'
     X  /'    ALLOWED :',5I7/'    REQUIRED:',5I7/)

	if(ENLAB>1e-5) then
 	t = hcm*K(PEL,EXL)
	NAMEV = 'OK'
	if(t > hktarg) NAMEV = 'NOT OK'
	write(ko,1486) ENLAB,t,NAMEV,hktarg
 1486 format(/' ACCURACY ANALYSIS at',f10.3,' MeV :'//
     x        '  Elastic h*k =',f6.3,' so ',a8,' compared with',f6.3/)
      DO 1467 I=1,NTHRESH
1467  IF(THRJ(I).gt..1) WRITE(KO,1468) RESM(I),THRJ(I)
1468  FORMAT('   Real(S-el) >',F5.2,' first at J =',F10.1)
      IF(JSWITCH.gt.0.1) WRITE(KO,1469) RMATR+GAP,JSWITCH
      WRITE(KO,1469) RTURN,JTOTAL
1469  FORMAT('   R-turn = ',F8.2, ' fm    at J =',F10.1)
C              give warning of limits for Coulomb excitations:
       T = 2. * ETA(PEL,EXL) * 180./PI
        TH1 = T / (JTOTAL+0.5)
        TH2 = T / (K(PEL,EXL) * max(abs(rmatch),rasym))
        WRITE(KO,1621) TH1,JTOTAL,TH2
 1621  FORMAT(/' Forward-angle excitation cut off ',
     X 'below',F8.3,' deg by max JT =',f10.1/
     x '                              and below',F8.3,' deg by max R.'/)

	do CP=1,NCP
	if(WID(CP)>0.) then
	NAMEV = 'OK'
	if(WID(CP) > RNL*1.2) NAMEV = 'NOT OK'
	WRITE(KO,220) CP,WID(CP),rnl,NAMEV,CENTR(CP),centre
220   FORMAT(/i2,
     &          ': Recommended RNL: non-local width >',F7.2,
     &              ' cf:',f7.2,' fm: ',a8/
     &        '    Recommended CENTRE: centration   ~',F7.2,
     &              ' cf:',f7.2,' fm'/)
	endif
     	enddo

	endif  ! acc
	endif  ! final

      TIME1 = SECOND() - TIME0
      IF(TIME1.GT.0.01) WRITE(KOI,1495) IAME,TIME1
 1495 FORMAT(' Total CPU',I3,' time = ',F10.2,' seconds')
      RETURN
CUNI  DEBUG SUBCHK
      END

c AMoro
c Triangle function defined by eg Joachain (2.100)
      function triangle(x,y,z)
      implicit none
      real*8 x,y,z,triangle
      triangle=(x-y-z)**2-4*y*z
      end function


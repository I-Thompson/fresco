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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PWISER(LVAL,ECM,XK,ETA,CLIST,NCLIST,NFLIST,
     *  PART,EXCIT,CHL,ACCRCY,RASYM,CFMAT,CGMAT,CRCRAT,NCH,M,MD,
     *  PWFLAG,LAMAX,STREN,LAMBDA,DERIV,FORMF,NF,NRM,GAP,
     *  ITC,SWITCH,FJSWTCH,SHOW,PRIN,GAIL,NSTEPM,ALLPART)
C     driver routine for CRCWFN from Fresco
C
	use parameters
	use io
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CHL(LMAX1,MXPEX,2),CZ,ZI,ZDL,FORMF(MAXM,NF)
      COMPLEX*16  CLIST(MAXCH,MAXCH,MCLIST)
      INTEGER LVAL(MAXCH),PART(MAXCH,3),C,J,ITC(MXP,MXX),
     *     EXCIT(MAXCH,3),LAMBDA(MLOC),CSTART,C1,C2,JF,LAM,SHOW,
     *     NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH),KMF
      LOGICAL PWFLAG(MXP),FJSWTCH,PRIN,DERIV,DONE,GAIL,ALLPART
      REAL*8 STREN(MLOC),XK(MXP,MXX),ETA(MXP,MXX),
     *    ECM(MAXCH,3),CFMAT(MAXCH,MAXCH,2),
     *    CGMAT(MAXCH,MAXCH,2),SWITCH
      REAL*8 PETA(MAXCH),PCOEF(MAXCH),PK(MAXCH),COEFO(MAXCH),
     *    COUPL(MAXCH,MAXCH,LAMAX+1),RIN(2),DC,CRCRAT(MAXCH)
      DATA DONE/.false./
C
      Z = 0.0
      if(SHOW>0) write(48,*) ' PWISER: ALLPART =',ALLPART
C
      CZ = (0D0,0D0)
      ZI = (0D0,1D0)
      if(.not.GAIL) then
      CFMAT(:,:,:)=0.0D0
      CGMAT(:,:,:)=0.0D0
      endif
C
      J=0
      CSTART = 1
C     SWITCH=1.0D3
      KMAX = LAMAX+1
      NSTEPM = 0
C
      DO 100 C=1,NCH
       J=J+1
       JST = CSTART-1
C                     JST is awkward off-set for partition.NE.1
C		      CFMAT,CGMAT have true indices,
C		      all other arrays in CRCWFN have block indices.
        IC = PART(C,1)
        IA = EXCIT(C,1)
C                                  make arrays for CRCWFN
        IF (PWFLAG(IC)) THEN
            L=LVAL(C)
            PETA(J)= ETA(IC,IA)
            PK(J) = XK(IC,IA)**2
            PCOEF(J)= 1.0D0
            COEFO(J)= -1.0D0/ECM(C,2)
	    HCM = ECM(C,3)
            RIN(2)= (M-1)*HCM
            RIN(1)= (MD-1)*HCM
            if(DERIV) RIN(1)=(NRM-1)*HCM
            if(DERIV) RIN(2)=RIN(1)
c                                 if match CRCWFN to zero
            IF (FJSWTCH) THEN
                RTURN =(PETA(1)+SQRT(PETA(1)**2+L*(L+1d0)))/SQRT(PK(1))
!!?		if(.not.DERIV) then
                RIN(2) = min(RTURN,RASYM)-GAP
                RIN(1) = RIN(2)-0.5d0
!!?		endif
!                write(KO,*) 'PWI:',RIN(2),RTURN,RASYM,GAP
c
               tkr1 = sqrt(pk(j))*rin(1)
               tkr2 = sqrt(pk(j))*rin(2)
               rnk1 = tkr1*peta(j)
               rnk2 = tkr2*peta(j)
c
               tmp1 = (L*real(L)+2.d0*rnk1-tkr1**2)
               tmp2 = (L*real(L)+2.d0*rnk2-tkr2**2)
               if (tmp1.le.0.d0.or.tmp2.le.0.d0) then
                    write(KO,*) 'Problem with CRCWFN match to zero'
                    tmp6 =1.5d0
                    goto 778 
               endif
              tmp1=sqrt(tmp1)
               tmp2=sqrt(tmp2)
c
               tmp3 = rin(2)*(L*tmp1+rnk1+L*real(L))/
     *               (rin(1)*(L*tmp2+rnk2+L*real(L)) )
               tmp5 = (tmp1*(tkr2-peta(j)) - tmp2*(tkr1-peta(j)))/
     *                 (L*real(L)+peta(j)**2)
c
               if (tmp3.le.0.d0.or.abs(tmp5).gt.1.d0) then
                    write(KO,*) 'Problem with CRCWFN match to zero'
                    tmp6 =1.5d0
                    goto 778 
               endif
              tmp3=L* log(tmp3)
               tmp5 = asin(tmp5)
c
              tmp4= sqrt( rin(2)*tmp1/(rin(1)*tmp2))
c
               tmp6= exp(tmp2-tmp1+tmp3+tmp5)*tmp4
c
778            CRCRAT(c) = tmp6 
c              write(KO,*) 'frx11: CRCRAT',c,CRCRAT(c) 
            ENDIF
            COUPL(J,J,1) =  2.D0*PETA(J)*XK(IC,IA)
            COUPL(J,J,2) =  L*(L+1d0)
        ENDIF
C
 777    FORMAT(' Rturn=',F8.2,' fm: match CRCWFN to zero',
     *  ' at',F8.2,' fm., in partition',I2)
       ICN=0
       IF (C.LT.NCH) ICN = PART(C+1,1)
C                                          detects change of partition
       IF (C.EQ.NCH.OR.(IC.NE.ICN.and..not.ALLPART)) THEN
C
         IF(FJSWTCH.and.(SHOW.ge.1.or..not.DONE) .and.CSTART==1)
     *       WRITE(KO,777) RTURN,RIN(2),PART(CSTART,1)
!	  DONE = FJSWTCH
	 KMF = 1
         DO 10 C1=CSTART,C
         IT  = ITC(PART(C1,1),EXCIT(C1,1))
         DO 10 C2=CSTART,C
	   J1=C1-JST
	   J2=C2-JST
C                                           diagonal couplings
           IF (C1.EQ.C2) THEN
               DO 11 LK=3,KMAX
 11            COUPL(J1,J2,LK) =0.0D0
              IF(PWFLAG(IC)) THEN
C                                         reorientation potl
               DO 12 NC=1,NCLIST(C1,C1)
	 	  JF = NFLIST(C1,C1,NC)
                 LAM = LAMBDA(JF)
                 IF (LAM.GT.0.AND.LAM.le.LAMAX) THEN
                     DC= CLIST(C1,C1,NC)
                   IF(ABS(DC)>1D-8.and.ABS(STREN(JF))>1D-20) then
		     COUPL(J1,J2,LAM+1)= COUPL(J1,J2,LAM+1)
     X			 + STREN(JF)*DC/(-ECM(C1,2))
		     KMF = max(KMF,LAM+1)
		     endif
                 ENDIF
 12            CONTINUE
C
      	    if(.not.GAIL) then
              CFMAT(C1,C2,1)= -ZI*(CHL(LVAL(C1)+1,IT,1))
              CFMAT(C1,C2,2)= -ZI*(CHL(LVAL(C1)+1,IT,2))
              CGMAT(C1,C2,1)= CHL(LVAL(C1)+1,IT,1)
              CGMAT(C1,C2,2)= CHL(LVAL(C1)+1,IT,2)
	    endif
C
             ELSE  if(.not.GAIL) then
              CFMAT(C1,C2,1)= -ZI*(CHL(LVAL(C1)+1,IT,2))
              CFMAT(C1,C2,2)= -ZI*(CHL(LVAL(C1)+1,IT,1))
              CGMAT(C1,C2,1)= CHL(LVAL(C1)+1,IT,2)
              CGMAT(C1,C2,2)= CHL(LVAL(C1)+1,IT,1)
              ENDIF
           ELSE
C                                          off-diagonal
               DO 14 LK=1,KMAX
 14            COUPL(J1,J2,LK) = 0.0D0
              IF (PWFLAG(IC)) THEN
C                                          coupling potls
                    ZDL = ZI**(LVAL(C1)-LVAL(C2))
               DO 18 NC=1,NCLIST(C1,C2)
	 	  JF = NFLIST(C1,C2,NC)
                  LAM=LAMBDA(JF)
                  IF (LAM.GT.0.AND.LAM.le.LAMAX) THEN
                    DC= CLIST(C1,C2,NC)*ZDL
                    IF (ABS(DC).GT.1.D-8) then
                     COUPL(J1,J2,LAM+1)= COUPL(J1,J2,LAM+1) +
     *                  STREN(JF)*DC/(-ECM(C1,2))
		     KMF = max(KMF,LAM+1)
                    ENDIF
                  ENDIF
 18            CONTINUE
              ENDIF
C
           ENDIF
 10       CONTINUE
C
         IF (PWFLAG(IC)) THEN

         IF(SHOW.gt.4)then
           DO C1=CSTART,C
           DO C2=CSTART,C
	   J1=C1-JST
	   J2=C2-JST
	write(KO,*) 'COUPL(',C1,',',C2,') =',(COUPL(J1,J2,I),I=1,KMAX)
	   enddo
	   enddo
	  endif
C
        CALL CRCWFN(C-JST,LVAL(CSTART),PETA,PK,KMAX,COUPL,PCOEF,
     *    RIN,RASYM,ACCRCY,CFMAT,CGMAT,MAXCH,PRIN,KMF,
     *    .FALSE.,.TRUE.,DERIV,.FALSE.,SWITCH,CSTART,HCM,FORMF,NSTEPD,
     * 	  MAXM,NF,M,CLIST,NCLIST,NFLIST,MAXCH,MCLIST,COEFO,SHOW)
      	  NSTEPM = max(NSTEPM,NSTEPD)
C
C   flags  DIAG CORREC CFG : set to F T F
C
            DO 40 C1=CSTART,C
           DO 40 C2=CSTART,C
         IF(SHOW.lt.2) GO TO 40
         IF(ABS(CFMAT(C1,C2,1)).GT.1E-9 .or. SHOW.ge.3)
     X    WRITE(KO,39) C1,C2,(CFMAT(C1,C2,I),CGMAT(C1,C2,I),I=1,2)
 39       FORMAT('CRCWFN to',I3,' from ch.',I3,':',1P,4G12.4)
 40        CONTINUE
C
          ENDIF
C                     zero counter for next partition
          J=0
          CSTART=C+1
        ENDIF
C
 100  CONTINUE
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CRCWFN(NCH,LVAL,ETA,ECM,KMAX,COUPL,COEF,RIN,ROUT,
     *    ACCRCY,FCC,GCC,MCH,PRIN,KMF,
     *    DIAG,CORREC,DERIV,CFG,SWITCH,CSTART,HCM,FORMF,NSTEPD,
     *    MAXM,NF,M,CLIST,NCLIST,NFLIST,MAXCH,MCLIST,COEFO,SHOW)
C
C     Coupled Real Coulomb Wavefunctions.
C     -----------------------------------
C
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ETA(MCH),ECM(MCH),COUPL(MAXCH,MAXCH,KMAX),COEF(MCH),
     *    FCC(MAXCH,MAXCH,2),GCC(MAXCH,MAXCH,2),ACCRCY,RIN(2),ROUT,
     $    RMAX,COEFO(MCH)
       real*8,save,allocatable:: RSTART(:)
       complex*16 FORMF(MAXM,NF),CLIST(MAXCH,MAXCH,MCLIST),ZDL
      INTEGER LVAL(MCH),CSTART,SHOW,
     *     NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH)
      LOGICAL DIAG,CORREC,DERIV,CFG,PRIN
C
C
	RMAX = (M-1)*HCM
        JST = CSTART-1
      IF (CFG) THEN
        LMAX1 = 1
        DO 20 I=1,NCH
 20         IF (LVAL(I).GT.LMAX1-1) LMAX1 = LVAL(I)+1
        LMAX1=LMAX1 + 4
        CALL COUL(MCH,NCH,ETA,ECM,COEF,FCC,GCC,ROUT,LMAX1,LVAL)
      ENDIF
C
!			Do scan twice: first to find NSTEPS & allocate STEPS,
!				       second to store in RSTART
	NSTEPS=0
	NSTEPD=0
	allocate (RSTART(1))
      do iscan=1,2
      CALL CUTUP(NCH,KMAX,COEF,RIN,ROUT,ACCRCY,RSTART,KMF,DERIV,MCH,
     *           iscan,NSTEPS,NSTEPD,NSTEPR2,SWITCH,SHOW,PRIN,MAXCH)
	NSTEPD = NSTEPS
	if(iscan.eq.1) then
	  deallocate (RSTART)
	  allocate (RSTART(NSTEPD+1))
	  endif
      enddo
C
      IMAX = 2*NCH*NSTEPS
      MAXA = 4*NCH**2*(2*NSTEPS-1)
      IF(SHOW.ge.1) WRITE(KO,*)'NO. OF STEPS ',NSTEPS,NSTEPR2,DERIV,IMAX
      IF(SHOW.ge.10) then
	IA = 3
	IB = 2
!	do R=20.,60.,1.0
!	do R=1.,20.,0.1
        do IR=10,200
        R=IR*0.1
	VT = VVAL(IA,IB,R)
	VM = VVALM(IA,IB,R,R+HCM)
	V1 = VVAL1(IA,IB,R)
	V2 = VVAL2(IA,IB,R)
	write(110,*) R,VT,VM,V1,V2
	enddo
	do I=1,M
	write(111,333) (I-1)*HCM,FORMF(I,1)
333	format(f9.3,2g15.5)
	enddo
	endif
C
      CALL MATRIX(NCH,ECM,KMAX,COEF,FCC,GCC,MCH,
     *     DIAG,CORREC,DERIV,NSTEPS,NSTEPR2,RSTART,
     *     MAXA,IMAX,SWITCH,MAXCH,CSTART,SHOW)
C
	deallocate (RSTART)
      RETURN
      CONTAINS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	FUNCTION VVAL(IA,IB,R)
        IMPLICIT REAL*8(A-H,O-Z)
	INTEGER C1,C2
C					Coupling 
	if(R.ge.RMAX) then
        VT = 0.0D0
	 RR = 1d0
        DO 23 NK=1,KMAX
         RR = RR/R
23        VT = VT + COUPL(IA,IB,NK)*RR
	  VT = VT/COEF(IA)
	else
	   C1 = IA+JST
	   C2 = IB+JST
	   VT = 0.0d0
	   if(C1.eq.C2) VT = LVAL(IA)*(LVAL(IA)+1)/R**2/COEF(IA)
	   I  = nint(R/HCM)+1
               ZDL = (0.,1.)**(LVAL(IA)-LVAL(IB))
               DO 18 NC=1,NCLIST(C1,C2)
	 	  JF = NFLIST(C1,C2,NC)
                    DC= CLIST(C1,C2,NC)*ZDL
		    VT = VT + FORMF(I,JF)*DC*COEFO(IA)
 18            CONTINUE
	endif
	VVAL = VT
	RETURN
	END FUNCTION VVAL

	FUNCTION VVALM(IA,IB,R0,R1)
        IMPLICIT REAL*8(A-H,O-Z)
	INTEGER C1,C2
C					Mean Coupling 
	if(R.ge.RMAX) then
            VBAR = COUPL(IA,IB,1)*LOG(R1/R0)
            DO 15 NK= 2,KMAX
               VBAR = VBAR + COUPL(IA,IB,NK)/real(NK-1)*
     *                   (R0**(1-NK) - R1**(1-NK))
 15         CONTINUE
            VBAR =  VBAR/((R1-R0)*COEF(IA))
	else
	   R = (R0+R1)*0.5
	   VBAR = VVAL(IA,IB,R)
	endif
	VVALM = VBAR
	RETURN
	END FUNCTION VVALM

	FUNCTION VVAL1(IA,IB,R)
        IMPLICIT REAL*8(A-H,O-Z)
	INTEGER C1,C2
C					First derivative
	if(R.ge.RMAX) then
            DDV = 0.0D0
            DO 50 K=1,KMAX
50             DDV = DDV - K * COUPL(IA,IB,K)
     *           /(COEF(IA)*R**(K+1))
	else
	   C1 = IA+JST
	   C2 = IB+JST
	   DDV = 0.0d0
	   if(C1.eq.C2) DDV = -2.*LVAL(IA)*(LVAL(IA)+1)/R**3/COEF(IA)
	   I  = min(M-2, nint(R/HCM)+1 )
           ZDL = (0.,1.)**(LVAL(IA)-LVAL(IB))
               DO 18 NC=1,NCLIST(C1,C2)
	 	  JF = NFLIST(C1,C2,NC)
                    DC= CLIST(C1,C2,NC)*ZDL
		    DDV = DDV + (FORMF(I+1,JF)-FORMF(I-1,JF))
     X			*DC*COEFO(IA) / (2.0*HCM)
 18            CONTINUE
	endif
	VVAL1 = DDV
	RETURN
	END FUNCTION VVAL1

	FUNCTION VVAL2(IA,IB,R)
        IMPLICIT REAL*8(A-H,O-Z)
	INTEGER C1,C2
C					Second derivative
	if(R.ge.RMAX) then
            DDV = 0.0D0
            DO 30 K=1,KMAX
 30            DDV = DDV +  K*(K+1) * COUPL(IA,IB,K)
     *           /(COEF(IA)*R**(K+2))
	else
	   C1 = IA+JST
	   C2 = IB+JST
	   DDV = 0.0d0
	   if(C1.eq.C2) DDV = 6.*LVAL(IA)*(LVAL(IA)+1)/R**4/COEF(IA)
	   I  = min(M-2, nint(R/HCM)+1 )
           ZDL = (0.,1.)**(LVAL(IA)-LVAL(IB))
               DO 18 NC=1,NCLIST(C1,C2)
	 	  JF = NFLIST(C1,C2,NC)
                    DC= CLIST(C1,C2,NC)*ZDL
		    DDV = (FORMF(I+1,JF)-2.*FORMF(I,JF)+FORMF(I-1,JF))
     X			*DC*COEFO(IA) / HCM**2 + DDV
 18            CONTINUE
	endif
	VVAL2 = DDV
	RETURN
	END FUNCTION VVAL2

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CUTUP(NCH,KMAX,COEF,RIN,ROUT,ACCRCY,RSTART,KMF,
     *    DERIV,MCH,iscan,NSTEPS,NSTEPD,NSTEPR2,SWITCH,SHOW,PRIN,MAXCH)
C     devide integration range into piecewise segments
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 COEF(MCH),ACCRCY,RIN(2),ROUT
       real*8 RSTART(NSTEPD+1)
      LOGICAL DERIV,FLAGR2,PRIN
      INTEGER SHOW
C
      FLAGR2 = .NOT. DERIV
      IF (ABS(RIN(1)-RIN(2)).LT.1D-6)  THEN
         FLAGR2=.FALSE.
         NSTEPR2=0
      ENDIF
C
C---------------------------- DIVIDE INTEGRATION RANGE------------
C

      R = RIN(1)
      HO = 0.1D0*R
      NSTEPS =0
      if(iscan.eq.2) RSTART(NSTEPS+1)=R
C
 10   NSTEPS = NSTEPS+1
      DDVM = 1.0D-30
      IF (R.LT.SWITCH) THEN
C        Airy functions step depends on 2nd derivative
         DO 20 IC=1,NCH
	    DDV = VVAL2(IC,IC,R)
 20         IF (ABS(DDV).GT.DDVM) DDVM = ABS(DDV)
C
         H = 4.0D0* SQRT(ACCRCY/DDVM)
      ELSE
C        sines steps depend on 1st derivative
         DO 40 IC=1,NCH
	    DDV = VVAL1(IC,IC,R)
40          IF (ABS(DDV).GT.DDVM) DDVM = ABS(DDV)
C        arbitrary factor 16
         H = 1.6D1*ACCRCY/DDVM
      ENDIF
C
C     ensure next step is less than 1.4 * last
      IF (H.GT.HO) H= HO
      HO = H*1.4D0
C
      IF (FLAGR2) THEN
          IF (R+H.GE.RIN(2).AND.RIN(2).GT.RIN(1)) THEN
              H =RIN(2)-R
              NSTEPR2 = NSTEPS+1
              FLAGR2 = .FALSE.
          ENDIF
      ENDIF
C
C   ENSURE AT LEAST TWO STEPS:
      IF(NSTEPS.EQ.1.AND.R+H.GT.ROUT) H=(ROUT-RIN(1))*0.5D0
C
      R = R+H
      IF (R.LT.ROUT) THEN
C         ensure that we don't get a tiny last step
          IF (ROUT-R.LT.0.6D0*H) R= ROUT -0.6D0*H
          if(iscan.eq.2) RSTART(NSTEPS+1) = R
          GOTO 10
      ENDIF
	if(iscan.eq.1) then
	  if(PRIN.and.KMF>1) then
	  WRITE(KO,98) NSTEPS,RIN(2),ROUT,KMF
98    FORMAT(' CRCWFN: ',I5,' steps from',F10.3,' to',F10.3,' fm.',
     X		' for powers <=',I2)
	  call flush(6)
	  endif
	  return
	endif
      RSTART(NSTEPS+1) = ROUT
      IF(SHOW.ge.2) WRITE(KO,99) NSTEPS,RIN,ROUT,SWITCH,
     X		(RSTART(I),I=1,NSTEPS+1)
99    FORMAT(' CUTUP: ',I5,' steps from',2F10.3,' to',F10.3,
     X      ';',f10.1,:,', at',/,(1X,12F10.3))
      RETURN
      END SUBROUTINE CUTUP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE MATRIX(NCH,ECM,KMAX,COEF,FCC,GCC,
     *      MCH,DIAG,CORREC,DERIV,NSTEPS,NSTEPR2,RSTART,
     *      MAXA,IMAX,SWITCH,MAXCH,CSTART,SHOW)
C     formulate and solve boundary condition matrix
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (S1=21.0D0, S2=770.0D0, S3=50666.0D0, PI=3.1415926536D0)
      REAL*8 ECM(MCH),COEF(MCH),
     *    FCC(MAXCH,MAXCH,2),GCC(MAXCH,MAXCH,2),ALPHA,RSTART(NSTEPS+1),
     *    AIR0(NCH,4),VBAR(NCH),BETA(NCH),A(NCH,NCH),
     *    B(NCH,NCH),C(NCH,NCH),RHS(IMAX),BC(12),AIR1(NCH,2)
      real*8,save,allocatable:: AMAT(:)
      INTEGER J0(IMAX),CSTART,SHOW
      LOGICAL DIAG,CORREC,DERIV,AIRY
C
      CALL JARRAY(J0,NCH,IMAX)
      allocate (AMAT(MAXA))
      do 1 I=1,MAXA
1	AMAT(I) = 0d0
      AIRY = .TRUE.
C
C************************************** LARGE LOOP OVER EACH RADIAL STEP
      DO 100 NS=1,NSTEPS
         R0 = RSTART(NS)
         R1 = RSTART(NS+1)
         RH = R1-R0
         RB = R0 + 0.5D0*RH
         IF (R0.GE.SWITCH) AIRY = .FALSE.
C------------------------------------- LSQ PART 1
         RLS = RH/20.0D0
         SUM = S1
         SUMXX = S2*RLS**2
         SUM4X = S3*RLS**4
         TED= 1.0D0/(SUM*SUM4X- SUMXX**2)
         D11 = SUM4X*TED
         D13 = -SUMXX*TED
         D22 = 1.0D0/SUMXX
         D33 = SUM*TED
C      ted=1/determinant
C---------------------------------- DVBAR---------
C
         IF (AIRY) THEN
C           dvbar is gradient at mid-point
            DVBAR = 0.0D0
            DO 10 NC=1,NCH
10          DVBAR = DVBAR + VVAL1(NC,NC,RB)/NCH
            IF(DVBAR.LT.0.0) THEN
               ALPHA = -1.0D0*(ABS(DVBAR)**(1.0D0/3.0D0))
            ELSE
               ALPHA = DVBAR**(1.0D0/3.0D0)
            ENDIF
         ENDIF
C
C---------------------------------------VBAR------
C
         DO 30 NC=1,NCH
C           vbar is average potential
            VBAR(NC) =  VVALM(NC,NC,R0,R1)
            IF (AIRY) THEN
              BETA(NC) = (VBAR(NC)-ECM(NC)/COEF(NC))/DVBAR - RB
            ELSE
             TEMP = ECM(NC)/COEF(NC) - VBAR(NC)
              IF (TEMP.LT.0.0) THEN
                 WRITE(KO,*)'E < V at SWITCHING RADIUS with R>SWITCH'
                 WRITE(KO,*)'Use Airy functions (ie Increase SWITCH)'
                 STOP
              ENDIF
              BETA(NC) = SQRT(TEMP)
              ALPHA =1.0D0
            ENDIF
!	write(KO,*) 'NS,R0,R1,ALPHA,BETA(1) =',
!     x      NS,R0,R1,ALPHA,BETA(1)
C
C--------------------------------------  SUMMING V-VO
C
            DO 25 NC2= 1,NCH
               SUMF = 0.0D0
               SUMFX = 0.0D0
               SUMFXX = 0.0D0
C      super-fast 21 point least-squares fit to U-Uo
               DO 22 NL=-10,10
                  RD = NL*RLS
	          VT = VVAL(NC,NC2,RB+RD)
                  IF (NC.EQ.NC2) THEN
                    VT = VT -VBAR(NC)
                    IF (AIRY) VT = VT -RD*DVBAR
                  ENDIF
                  SUMF = SUMF+ VT
                  SUMFX = SUMFX+ VT*RD
                  SUMFXX = SUMFXX + VT*RD**2
 22            CONTINUE
C-------------------------------------------- LSQ PART 2 -------
               AT = D11*SUMF + D13*SUMFXX
               BT = D22*SUMFX
               CT = D13*SUMF + D33*SUMFXX
               IF (AIRY) THEN
                  A(NC,NC2) = AT - BT*RB +CT*RB*RB
                  B(NC,NC2) = BT - (CT+CT)*RB
                  C(NC,NC2) = CT
               ELSE
                  A(NC,NC2) = AT
                  B(NC,NC2) = BT
                  C(NC,NC2) = CT
               ENDIF
C
C---------------------------------------------------------------
 25         CONTINUE
 30      CONTINUE
C
         IA=(NS-2)*2*NCH
         JA= IA+2*NCH
         JB= JA + NCH
         ID= IA + NCH
         IC= ID + NCH
         IDC= IC + NCH
C-------------------------------------calculate and store Airy's
         DO 40 NC=1,NCH
C
            IF (AIRY) THEN
               DX0 = ALPHA*(R0+BETA(NC))
               DX1 = ALPHA*(R1+BETA(NC))
               CALL DAIRY(DX0,AI0,AD0,BI0,BD0)
               CALL DAIRY(DX1,AI1,AD1,BI1,BD1)
            ELSE
               N2PI = INT(R0/(2.D0*PI))
               R0 = R0 - N2PI*2.0D0*PI
               R1 = R1 - N2PI*2.0D0*PI
               RB = (R1+R0)/2.D0
               DX0 = BETA(NC)*R0
               DX1 = BETA(NC)*R1
               AI0 = COS(DX0)
               AI1 = COS(DX1)
               BI0 = SIN(DX0)
               BI1 = SIN(DX1)
               AD0 = -BETA(NC)*BI0
               AD1 = -BETA(NC)*BI1
               BD0 =  BETA(NC)*AI0
               BD1 =  BETA(NC)*AI1
            ENDIF
C
            IF (IA.LT.0) THEN
C              values at Rmin not needed in the matrix
               AIR0(NC,1) = AI0
               AIR0(NC,2) = BI0
               AIR0(NC,3) = AD0
               AIR0(NC,4) = BD0
!	write(KO,*) 'AIR0:',NC,AI0,BI0,AD0,BD0,'@',ALPHA,BETA(NC)
               JA = 0
               JB = NCH
            ELSE
               AMAT(J0(IA+NC)+JA+NC) =-AI0
               AMAT(J0(IA+NC)+JB+NC) =-BI0
               AMAT(J0(ID+NC)+JA+NC) =-AD0
               AMAT(J0(ID+NC)+JB+NC) =-BD0
            ENDIF
            AMAT(J0(IC+NC)+JA+NC) = AI1
            AMAT(J0(IC+NC)+JB+NC) = BI1
            AMAT(J0(IDC+NC)+JA+NC) = AD1
            AMAT(J0(IDC+NC)+JB+NC) = BD1
C
            IF (.NOT.DERIV) THEN
            IF (NS.EQ.NSTEPR2) THEN
               AIR1(NC,1) = AI0
               AIR1(NC,2) = BI0
            ENDIF
            ENDIF
 40      CONTINUE
C
         IF (AIRY) PAL = PI/ALPHA
C----------------------------------  CORRECTION LOOP
         DO 50 NC=1,NCH
C
            IF (IA.LT.0) THEN
               AI0 = AIR0(NC,1)
               BI0 = AIR0(NC,2)
               AD0 = AIR0(NC,3)
               BD0 = AIR0(NC,4)
            ELSE
               AI0 =-AMAT(J0(IA+NC)+JA+NC)
               BI0 =-AMAT(J0(IA+NC)+JB+NC)
               AD0 =-AMAT(J0(ID+NC)+JA+NC)
               BD0 =-AMAT(J0(ID+NC)+JB+NC)
            ENDIF
            AI1 = AMAT(J0(IC+NC)+JA+NC)
            BI1 = AMAT(J0(IC+NC)+JB+NC)
            AD1 = AMAT(J0(IDC+NC)+JA+NC)
            BD1 = AMAT(J0(IDC+NC)+JB+NC)
C
            BET1 = BETA(NC)
           IF (CORREC) THEN
C---------------------------------------- DO OFF DIAGONAL'S FIRST
              DO 45 NC2=NC+1,NCH
                 BET2 = BETA(NC2)
                 IF (IA.LT.0) THEN
                    AI2 = AIR0(NC2,1)
                    BI2 = AIR0(NC2,2)
                    AD2 = AIR0(NC2,3)
                    BD2 = AIR0(NC2,4)
                 ELSE
                    AI2 =-AMAT(J0(IA+NC2)+JA+NC2)
                    BI2 =-AMAT(J0(IA+NC2)+JB+NC2)
                    AD2 =-AMAT(J0(ID+NC2)+JA+NC2)
                    BD2 =-AMAT(J0(ID+NC2)+JB+NC2)
                 ENDIF
                 AI3 = AMAT(J0(IC+NC2)+JA+NC2)
                 BI3 = AMAT(J0(IC+NC2)+JB+NC2)
                 AD3 = AMAT(J0(IDC+NC2)+JA+NC2)
                 BD3 = AMAT(J0(IDC+NC2)+JB+NC2)
C
             IF (AIRY) THEN
               CALL OFDAG(R0,R1,ALPHA,BET1,BET2,BC,AI0,BI0,AD0,BD0,
     *           AI1,BI1,AD1,BD1,AI2,BI2,AD2,BD2,AI3,BI3,AD3,BD3,SHOW)
             ELSE
               CALL SCOFF(R0,R1,BET1,BET2,BC,AI0,BI0,AI1,BI1,
     *           AI2,BI2,AI3,BI3)
             ENDIF
C
                 W = A(NC,NC2)
                 WR= B(NC,NC2)
                 WRR = C(NC,NC2)
                 IF (.NOT.AIRY) THEN
                    W = W - WR*RB +WRR*RB*RB
                    WR = WR - 2.D0*WRR*RB
                    PAL = 1.D0/BET1
                 ENDIF
      WP371  =    PAL*(W*BC(3)+WR*BC(7)+WRR*BC(11))
      WP159  =    PAL*(W*BC(1)+WR*BC(5)+WRR*BC(9))
      WP261  =    PAL*(W*BC(2)+WR*BC(6)+WRR*BC(10))
      WP481  =    PAL*(W*BC(4)+WR*BC(8)+WRR*BC(12))
C
C-------------------------------------correct (NC,NC2)
      AMAT(J0(IC+NC)+JA+NC2)= -AI1*WP371+ BI1*WP159
      AMAT(J0(IDC+NC)+JA+NC2)=(-AD1*WP371+ BD1*WP159)*ALPHA
      AMAT(J0(IC+NC)+JB+NC2)=  BI1*WP261- AI1*WP481
      AMAT(J0(IDC+NC)+JB+NC2)= (BD1*WP261- AD1*WP481)*ALPHA
C
                 W = A(NC2,NC)
                 WR= B(NC2,NC)
                 WRR = C(NC2,NC)
                 IF (.NOT.AIRY) THEN
                    W = W - WR*RB +WRR*RB*RB
                    WR = WR - 2.D0*WRR*RB
                    PAL = 1.D0/BET2
                 ENDIF
C--------------------------------------correct (NC2,NC)
      WP371  =    PAL*(W*BC(3)+WR*BC(7)+WRR*BC(11))
      WP159  =    PAL*(W*BC(1)+WR*BC(5)+WRR*BC(9))
      WP261  =    PAL*(W*BC(2)+WR*BC(6)+WRR*BC(10))
      WP481  =    PAL*(W*BC(4)+WR*BC(8)+WRR*BC(12))
C
C-------------------------------------correct (NC,NC2)
      AMAT(J0(IC+NC2)+JA+NC)= -AI3*WP261+ BI3*WP159
      AMAT(J0(IDC+NC2)+JA+NC)=(-AD3*WP261+ BD3*WP159)*ALPHA
      AMAT(J0(IC+NC2)+JB+NC)=  BI3*WP371- AI3*WP481
      AMAT(J0(IDC+NC2)+JB+NC)=(BD3*WP371- AD3*WP481)*ALPHA
C
C
 45           CONTINUE
C-------------------------------- diagonal corrections
           ENDIF
            IF (.NOT.CORREC) THEN
                DO 46 JI=1,12
 46                BC(JI)=0.0D0
                PAL = 0.0D0
            ELSE
               IF (AIRY) THEN
                 CALL CDIAG(R0,R1,ALPHA,BET1,BC,AI0,BI0,AD0,BD0,
     *            AI1,BI1,AD1,BD1)
               ELSE
                 CALL SCORR(R0,R1,BET1,BC,AI0,BI0,AI1,BI1)
               ENDIF
C
            ENDIF
            W = A(NC,NC)
            WR= B(NC,NC)
            WRR = C(NC,NC)
            IF (.NOT.AIRY) THEN
               W = W - WR*RB +WRR*RB*RB
               WR = WR - 2.D0*WRR*RB
               PAL = 1.D0/BET1
            ENDIF
C
            WP258 = PAL*(W*BC(2)+WR*BC(5)+WRR*BC(8))
            WP147 = PAL*(W*BC(1)+WR*BC(4)+WRR*BC(7))
            WP369 = PAL*(W*BC(3)+WR*BC(6)+WRR*BC(9))
C
      AMAT(J0(IC+NC)+JA+NC)= AI1*(1.0D0-WP258)+BI1*WP147
      AMAT(J0(IDC+NC)+JA+NC)=(AD1*(1.0D0-WP258) +BD1*WP147)*ALPHA
      AMAT(J0(IC+NC)+JB+NC)= BI1*(1.0D0+WP258) -AI1*WP369
      AMAT(J0(IDC+NC)+JB+NC)=(BD1*(1.0D0+WP258) -AD1*WP369)*ALPHA
C
C----------------------------------- mult. deriv.s by alpha
          IF (AIRY) THEN
            IF (IA.LT.0) THEN
               AIR0(NC,3) = AIR0(NC,3)*ALPHA
               AIR0(NC,4) = AIR0(NC,4)*ALPHA
            ELSE
               AMAT(J0(ID+NC)+JA+NC) = AMAT(J0(ID+NC)+JA+NC)*ALPHA
               AMAT(J0(ID+NC)+JB+NC) = AMAT(J0(ID+NC)+JB+NC)*ALPHA
            ENDIF
          ENDIF
C
 50      CONTINUE
C--------------------------end of correction loop
 100  CONTINUE
C********************************************* end of Nsteps loop
           IF(SHOW.ge.5) WRITE(KO,101) (AMAT(I),I=1,MAXA)
 101       FORMAT(' AMAT:',1P,6E11.3)
C
      CALL LUDCST(AMAT,MAXA,IMAX,IMAX,J0,NCH)
C
C
C                     JST is awkward off-set for partition.NE.1
      JST = CSTART-1
      IST = IMAX-2*NCH
C----------------------------------- BACKSUBSTITUTE
      DO 140 NC=1,NCH
         DO 150 I=1,IMAX
 150        RHS(I) = 0.0D0
          IF(SHOW.ge.6) WRITE(KO,49) '0',NC,(AIR0(NC,I),I=1,4)
          IF(SHOW.ge.6) WRITE(KO,49) '1',NC,(AIR1(NC,I),I=1,2)
 49       FORMAT('AIR',A1,'(1-4) @',I3,':',1P,4E11.3)
         RHS(IST+NC) = FCC(NC+JST,NC+JST,1)
         RHS(IST+NC+NCH) = FCC(NC+JST,NC+JST,2)
           IF(SHOW.ge.5) WRITE(KO,151) NC,(RHS(I),I=1,IMAX)
 151       FORMAT('RHS @',I3,':',1p,6E11.3/(9X,1P,6E11.3))
         CALL LUBKST(AMAT,MAXA,IMAX,IMAX,J0,NCH,RHS)
           IF(SHOW.ge.5) WRITE(KO,152) NC,(RHS(I),I=1,IMAX)
 152       FORMAT('SOL @',I3,':',1p,6E11.3/(9X,1p,6E11.3))
           IF(SHOW.ge.5) WRITE(KO,152) NC,(RHS(I),I=1,IMAX)
         DO 160 NC2=1,NCH
           FCC(NC2+JST,NC+JST,1)= RHS(NC2)*AIR0(NC2,1)+
     *                            RHS(NCH+NC2)*AIR0(NC2,2)
           IF (DERIV) THEN
                FCC(NC2+JST,NC+JST,2)= RHS(NC2)*AIR0(NC2,3)
     *                        +RHS(NCH+NC2)*AIR0(NC2,4)
           ELSE
                IA = 2*(NSTEPR2-1)*NCH
C               IF (IA.GT.0)THEN
		  if(IA+NC2>IMAX) then
		  write(6,*) IA,NC2,IMAX,NSTEPR2,NCH,NSTEPS
		  stop
		  endif
                  FCC(NC2+JST,NC+JST,2)= RHS(IA+NC2)*AIR1(NC2,1)
     *                          +RHS(IA+NCH+NC2)*AIR1(NC2,2)
C               ELSE
C                 FCC(NC2+JST,NC+JST,2)=FCC(NC2+JST,NC+JST,1)
C               ENDIF
           ENDIF
        IF(SHOW.ge.2.and.ABS(FCC(NC2,NC,1)).GT.1E-9.or.SHOW.ge.3)
     X   WRITE(KO,39) 'F',NC2,NC,IA,(FCC(NC2,NC,I),I=1,2)
     X            ,RHS(IA+NC2),AIR1(NC2,1) ,RHS(IA+NCH+NC2),AIR1(NC2,2)
 39       FORMAT('MATRIX ',A1,' to',I3,' from ch.',I3,I5,':',1P,6E11.3)
 160     CONTINUE
C---------------------------------done the F's now do the G's
         DO 170 I=1,IMAX
 170        RHS(I) = 0.0D0
         RHS(IST+NC) = GCC(NC+JST,NC+JST,1)
         RHS(IST+NC+NCH) = GCC(NC+JST,NC+JST,2)
           IF(SHOW.ge.5) WRITE(KO,151) NC,(RHS(I),I=1,IMAX)
         CALL LUBKST(AMAT,MAXA,IMAX,IMAX,J0,NCH,RHS)
           IF(SHOW.ge.5) WRITE(KO,152) NC,(RHS(I),I=1,IMAX)
         DO 180 NC2=1,NCH
            GCC(NC2+JST,NC+JST,1)= RHS(NC2)*AIR0(NC2,1)+
     *                             RHS(NCH+NC2)*AIR0(NC2,2)
            IF (DERIV) THEN
                GCC(NC2+JST,NC+JST,2)= RHS(NC2)*AIR0(NC2,3)
     *                                +RHS(NCH+NC2)*AIR0(NC2,4)
            ELSE
C               IF (IA.GT.0) THEN
                   GCC(NC2+JST,NC+JST,2)= RHS(IA+NC2)*AIR1(NC2,1)
     *                           +RHS(IA+NCH+NC2)*AIR1(NC2,2)
C               ELSE
C                  GCC(NC2+JST,NC+JST,2)= GCC(NC2+JST,NC+JST,1)
C               ENDIF
            ENDIF
        IF(SHOW.ge.2.and.ABS(GCC(NC2,NC,1)).GT.1E-9.or.SHOW.ge.3)
     X   WRITE(KO,39) 'G',NC2,NC,IA,(GCC(NC2,NC,I),I=1,2)
 180     CONTINUE
C
 140  CONTINUE
C--------------------------------End of Story!
      deallocate (AMAT)
      RETURN
      END SUBROUTINE MATRIX
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE COUL(MCH,NCH,ETA,ECM,COEF,FCC,GCC,ROUT,LMAX1,LVAL)
Cget uncoupled Coulomb wfns - better to do in main program.
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ETA(MCH),ECM(MCH),COEF(MCH),ROUT,
     *    FCC(MCH,MCH,2),GCC(MCH,MCH,2),
     *    F(LMAX1),G(LMAX1),FP(LMAX1),GP(LMAX1)
      INTEGER LVAL(MCH)
C
      XLM =  0.0D0
      XLX =  LMAX1-1
      DO 100 NC =1,NCH
        AK = SQRT(ECM(NC)/COEF(NC))
        RHO = ROUT*AK
        KFN=0
        MODE=1
        IFAIL=0
        CALL COULFG(RHO,ETA(NC),XLM,XLX,F,G,FP,GP,MODE,KFN,IFAIL,M1)
        I = LVAL(NC)+1
        FCC(NC,NC,1) =F(I)
        FCC(NC,NC,2) =FP(I)*AK
        GCC(NC,NC,1) =G(I)
        GCC(NC,NC,2) =GP(I)*AK
100     CONTINUE
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CDIAG(R0,R1,ALP,BET,BC,AI0,BI0,AD0,BD0,AI1,BI1,
     *      AD1,BD1)
C     diagonal correction integralsfor Airy fns
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (PI=3.141592653589793D0)
C
      REAL*8 ALP,BET,BET2,BET3,PHA2,PHA3,AI0,AD0,BI0,BD0,
     +       AI1,AD1,BI1,BD1,C11,C12,C13,C01,C02,C03,CB,CC,
     +       PHA,BC,RB0,RB1
      REAL*8 R0,R02,R03,R1,R12,R13
      DIMENSION BC(12)
C
C     correction terms
      BET2 = BET*BET
      BET3 = BET2*BET
      PHA = 1.0D0/ALP
      PHA2 = PHA*PHA
      PHA3 = PHA2*PHA
      RB0 = R0+BET
      RB1 = R1+BET
      R02 = R0*R0
      R03 = R0*R02
      R12 = R1*R1
      R13 = R1*R12
      AI12 = AI1*AI1
      AI02 = AI0*AI0
      AD12 = AD1*AD1
      AD02 = AD0*AD0
      AID1 = AI1*AD1
      AID0 = AI0*AD0
      AIBI1= AI1*BI1
      AIBI0= AI0*BI0
      AIBD1= AI1*BD1
      AIBD0= AI0*BD0
      ADBI1= BI1*AD1
      ADBI0= BI0*AD0
      ADBD1 = AD1*BD1
      ADBD0 = AD0*BD0
      BI12 = BI1*BI1
      BI02 = BI0*BI0
      BID1 = BI1*BD1
      BID0 = BI0*BD0
      BD12 = BD1*BD1
      BD02 = BD0*BD0
C
      C11=(R1+BET)
      C13=- PHA
      C01=(R0+BET)
      C03= C13
C-----------------  int(AA), int(AB), int(BB)
      BC(1) = C11*AI12 + C13*AD12 - C01*AI02 - C03*AD02
      BC(2) = C11*AIBI1 + C13*ADBD1 - C01*AIBI0 - C03*ADBD0
      BC(3) = C11*BI12 + C13*BD12 - C01*BI02 - C03*BD02
C
      CB = 1.0D0/3.0D0
C     C11= CB*(RB1*RB1-3.0D0*RB1*BET)
      C11= CB*(R12-BET*R1-BET2-BET2)
      C12= CB*PHA2/2.0D0
      C13= CB*(BET+BET-R1)*PHA
      C01= CB*(R02-BET*R0-BET2-BET2)
      C02= C12
      C03= CB*(BET+BET-R0)*PHA
C
C------------------ int(RAA), int(RAB), int(RBB)
      BC(4) = (C13*AD12-C03*AD02) + (C11*AI12-C01*AI02)
     *      + 2.0D0 *(C12*AID1-C02*AID0)
      BC(5) = (C13*ADBD1-C03*ADBD0) + (C11*AIBI1-C01*AIBI0)
     *      + C12*(ADBI1+AIBD1) - C02*(ADBI0+AIBD0)
      BC(6) = (C13*BD12-C03*BD02) + (C11*BI12-C01*BI02)
     *      + 2.0D0 *(C12*BID1-C02*BID0)
C
      CC = 1.0D0/15.0D0
      B3P3 = 8.0D0*BET3 - 3.0D0*PHA3
      C11= CC*(3.0D0*R13 -R12*BET +4.0D0*R1*BET2 + B3P3)
      C12= CC*(3.0D0*R1 -2.0D0*BET)*PHA2
      C13= CC*(4.0D0*R1*BET -3.0D0*R12 -8.0D0*BET2)*PHA
      C01= CC*(3.0D0*R03 -R02*BET +4.0D0*R0*BET2 + B3P3)
      C02= CC*(3.0D0*R0 -2.0D0*BET)*PHA2
      C03= CC*(4.0D0*R0*BET -3.0D0*R02 -8.0D0*BET2)*PHA
C------------------ int(RRAA), int(RRAB), int(RRBB)
C
      BC(7)= (C13*AD12-C03*AD02) +(C11*AI12-C01*AI02)
     *     + 2.0D0* (C12*AID1-C02*AID0)
      BC(8)= (C13*ADBD1-C03*ADBD0) +(C11*AIBI1-C01*AIBI0)
     *     + C12*(AIBD1+ADBI1) - C02*(AIBD0+ADBI0)
      BC(9)= (C13*BD12-C03*BD02) +(C11*BI12-C01*BI02)
     *     + 2.0D0* (C12*BID1-C02*BID0)
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE OFDAG(R0,R1,ALP,BET1,BET2,BC,AI0,BI0,AD0,BD0,
     *           AI1,BI1,AD1,BD1,AI2,BI2,AD2,BD2,AI3,BI3,AD3,BD3,SHOW)
C     off diagonal couplingsfor Airy fns
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (PI=3.141592653589793D0)
C
      REAL*8 ALP,BET2,PHA2,PHA3,AI0,AD0,BI0,BD0,
     +       AI1,AD1,BI1,BD1,C11,C12,C13,C01,C02,C03,CB,CC,
     +       PHA,BC
      REAL*8 R0,R1
      INTEGER SHOW
      DIMENSION BC(12)
       LOGICAL redo
       redo = .false.
       degen = 1d-7
       degen = 1d-5
C
C     correction terms
      BPB = BET1+BET2
      BMB  = BET1-BET2
C
C     if degenerate
C
c      IF (ABS(BMB/BET1).LT.1.0D-3) THEN
 10    IF (ABS(BMB/BET1).LT.2*degen.or.redo) THEN
        CALL CDIAG(R0,R1,ALP,BET1,BC,AI0,BI0,AD0,BD0,
     *       AI1,BI1,AD1,BD1)
          BC(12)=BC(9)
          BC(11)=BC(8)
          BC(10)=BC(8)
          BC(9)=BC(7)
          BC(8)=BC(6)
          BC(7)=BC(5)
          BC(6)=BC(5)
          BC(5)=BC(4)
          BC(4)=BC(3)
          BC(3)=BC(2)
C         WRITE(KO,*) BC(9)
          RETURN
      ENDIF
C
      BMB2 = BMB*BMB
      DBMB = 1.0D0/BMB
      DBMB2 = DBMB*DBMB
      PHA = 1.0D0/ALP
      PHA2 = PHA*PHA
      PHA3 = PHA2*PHA
      P2DB = PHA2*DBMB
      P3DB2 = PHA3*DBMB2
      ADI13 = AD1*AI3 - AI1*AD3
      ADI02 = AD0*AI2 - AI0*AD2
      AB13  = AD1*BI3 - AI1*BD3
      AB02  = AD0*BI2 - AI0*BD2
      BA13  = BD1*AI3 - BI1*AD3
      BA02  = BD0*AI2 - BI0*AD2
      BDI13 = BD1*BI3 - BI1*BD3
      BDI02 = BD0*BI2 - BI0*BD2
C------------------------- int(AA), int(AB), int(BA), int(BB)
      C12= P2DB
      BC(1) = C12*(ADI13-ADI02)
      BC(2) = C12*(AB13-AB02)
      BC(3) = C12*(BA13-BA02)
      BC(4) = C12*(BDI13-BDI02)
C
      P5D3 = P3DB2*P2DB
      C11 = -(BPB+R1+R1)*P3DB2
      C01 = -(BPB+R0+R0)*P3DB2
      C12 = R1*P2DB + P5D3 +P5D3
      C02 = R0*P2DB + P5D3 +P5D3
      C13 = 2.0D0*PHA *P3DB2
      C03 = C13
C------------------------- int(RAA), int(RAB), int(RBA), int(RBB)
      BC(5) = C11*AI1*AI3 + C12*ADI13 + C13*AD1*AD3 -
     *        C01*AI0*AI2 - C02*ADI02 - C03*AD0*AD2
      BC(6) = C11*AI1*BI3 + C12*AB13  + C13*AD1*BD3 -
     *        C01*AI0*BI2 - C02*AB02  - C03*AD0*BD2
      BC(7) = C11*BI1*AI3 + C12*BA13  + C13*BD1*AD3 -
     *        C01*BI0*AI2 - C02*BA02  - C03*BD0*AD2
      BC(8) = C11*BI1*BI3 + C12*BDI13 + C13*BD1*BD3 -
     *        C01*BI0*BI2 - C02*BDI02 - C03*BD0*BD2
C
      ALP3 = ALP**3
      A3B2 = ALP3*BMB2
      P6D4 = P3DB2*P3DB2
      DB5 = DBMB2*DBMB2*DBMB *PHA3*PHA3*PHA2
      CC = 24.0D0 + 2.0D0*ALP3*BPB*BMB2
      C11 = -P6D4 * (12.0D0*BPB + CC*R1 +
     *       4.0D0*A3B2*R1*R1)
      C01 = -P6D4 * (12.0D0*BPB + CC*R0 +
     *       4.0D0*A3B2*R0*R0)
      C12 = DB5 * (CC+
     *      12.0D0*A3B2*R1 + A3B2*A3B2*R1*R1)
      C02 = DB5 * (CC+
     *      12.0D0*A3B2*R0 + A3B2*A3B2*R0*R0)
      C13 = 4.0D0*PHA*P6D4*(6.0D0+A3B2*R1)
      C03 = 4.0D0*PHA*P6D4*(6.0D0+A3B2*R0)
      CB = 4.0D0*P5D3
c
c     test for loss of precision: if so do diagonal coupling routine
c
       errloss = abs((c13-c03)/c13)
      if (errloss.lt.degen) then
           IF(SHOW.ge.2) write(KO,*)
     X          'Rounding error loss ',real(errloss),' < ',real(degen),
     X              ' so using diagonal integral'
          redo = .true.
           go to 10
      endif
C
C------------------------- int(RRAA), int(RRAB), int(RRBA), int(RRBB)
      BC(9) = C11*AI1*AI3 + C12*ADI13 + C13*AD1*AD3 -
     *        C01*AI0*AI2 - C02*ADI02 - C03*AD0*AD2 +
     *        CB*(BET2*(AD1*AI3-AD0*AI2) - BET1*(AI1*AD3-AI0*AD2))
      BC(10)= C11*AI1*BI3 + C12*AB13  + C13*AD1*BD3 -
     *        C01*AI0*BI2 - C02*AB02  - C03*AD0*BD2 +
     *        CB*(BET2*(AD1*BI3-AD0*BI2) - BET1*(AI1*BD3-AI0*BD2))
      BC(11)= C11*BI1*AI3 + C12*BA13  + C13*BD1*AD3 -
     *        C01*BI0*AI2 - C02*BA02  - C03*BD0*AD2 +
     *        CB*(BET2*(BD1*AI3-BD0*AI2) - BET1*(BI1*AD3-BI0*AD2))
      BC(12)= C11*BI1*BI3 + C12*BDI13 + C13*BD1*BD3 -
     *        C01*BI0*BI2 - C02*BDI02 - C03*BD0*BD2 +
     *        CB*(BET2*(BD1*BI3-BD0*BI2) - BET1*(BI1*BD3-BI0*BD2))
C
C     WRITE(KO,717)BC(9),BET1,BMB,BMB/BET1,R1
C717  FORMAT(5(D12.4,3X))
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
       SUBROUTINE JARRAY(J0,NCH,IMAX)
       INTEGER J0(IMAX)
C      New array technique for staircase matrix.
C
       DO 10 I=1,IMAX
          K= (I-1)/(2*NCH)
          J0(I) = (I-1)*4*NCH - 2*K*NCH
 10       CONTINUE
C
       I=IMAX-2*NCH+2
       JSUB =2*NCH
C
 20    J0(I) = J0(I) -JSUB
       I=I+1
       JSUB = JSUB + 2*NCH
       IF (I.LE.IMAX) GOTO 20
C
       END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      SUBROUTINE LUDCST(A,NP1,NP2,N,J0,NCH)
C
C Given an N * N matrix A, with physical dimension NP, this routine
C replaces it with the LU decomposition of a rowwise permutation of
C itself. Staircase modifications as 1-D array.
C NO PIVOTING!
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      PARAMETER (TINY=1.0D-35)
      DIMENSION A(NP1),J0(NP2)
      DO 19 J=1,N
         JMD = ((J-1)/(2*NCH)-1)*2*NCH
         DO 14 I=MAX(1,1+JMD),J-1
         SUM = A(J0(I)+J)
           IMD = ((I-1)/(2*NCH))*2*NCH
           DO 13 K=MAX(1,IMD+1,JMD+1),
     *             MIN(I-1,IMD+4*NCH,JMD+4*NCH)
13         SUM = SUM - A(J0(I)+K)*A(J0(K)+J)
14       A(J0(I)+J) = SUM
C
         DO 16 I=J,MIN(N,JMD+4*NCH)
            SUM = A(J0(I)+J)
            IMD = ((I-1)/(2*NCH))*2*NCH
            DO 15 K=MAX(1,IMD+1,JMD+1),
     *              MIN(J-1,IMD+4*NCH,JMD+4*NCH)
15            SUM = SUM - A(J0(I)+K)*A(J0(K)+J)
            A(J0(I)+J) = SUM
C
16          CONTINUE
C
         IF(ABS(A(J0(J)+J)).LT.TINY) A(J0(J)+J) = TINY
         IF(J.NE.N) THEN
            DUM = 1.0D0 /A(J0(J)+J)
            DO 18 I=J+1,MIN(N,JMD+4*NCH)
18             A(J0(I)+J) = A(J0(I)+J)*DUM
         ENDIF
19     CONTINUE
      RETURN
      END
      SUBROUTINE LUBKST(A,NP1,NP2,N,J0,NCH,B)
C
C Solves the set of N linear equations A.X = B.
C Here, A is input, not as the matrix A, but rather its LU decomposition
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION A(NP1),B(NP2),J0(NP2)
      II = 0
      DO 12 I=1,N
        SUM = B(I)
        IF(II.NE.0) THEN
          IMD = ((I-1)/(2*NCH))*2*NCH
          DO 11 J=IMD+1,I-1
11           SUM = SUM - A(J0(I)+J)*B(J)
          ELSE IF(ABS(SUM).GT.0.0) THEN
             II = I
          ENDIF
12        B(I) = SUM
      DO 14 I=N,1,-1
        SUM = B(I)
        IF(I.LT.N) THEN
           IMD = ((I-1)/(2*NCH))*2*NCH
           DO 13 J = I+1,MIN(N,IMD+4*NCH)
13            SUM = SUM - A(J0(I)+J)*B(J)
           ENDIF
14      B(I) = SUM/A(J0(I)+I)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SCORR(R0,R1,BET1,BC,AI0,BI0,AI1,BI1)
C     diagonal correction integrals for sines + cosines
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 R0,R1,BET1,BC(12),AI0,AI1,BI0,BI1
C
      SN20 = 2.D0*AI0*BI0
      SN21 = 2.D0*AI1*BI1
      CS20 = AI0*AI0 - BI0*BI0
      CS21 = AI1*AI1 - BI1*BI1
      R02 = R0 *R0
      R12 = R1 *R1
      BTE = 1.D0/BET1
      SMS = 0.5D0*BTE*(SN21-SN20)
      CMC = 0.5D0*BTE*(CS21-CS20)
      BT4 = 0.25D0*BTE
C
      BC(1) = 0.5D0*(R1-R0 + BTE*(AI1*BI1 - AI0*BI0))
      BC(2) = 0.5D0*BTE * (AI0*AI0 - AI1*AI1)
      BC(3) = 0.5D0*(R1-R0 - BTE*(AI1*BI1 - AI0*BI0))
C
      CR3 = 0.25D0*(R12-R02)
      CR4 = BT4*(R1*SN21-R0*SN20 +CMC)
      BC(4) = CR3 + CR4
      BC(6) = CR3 - CR4
      BC(5) = BT4*(R0*CS20-R1*CS21 +SMS)
C
      CR3 = (R1*R12-R0*R02)/6.0D0
      CR4 =  BT4*(R12*SN21-R02*SN20 + BTE*
     *    (R1*CS21 - R0*CS20  - SMS))
      BC(7) = CR3 + CR4
      BC(9) = CR3 - CR4
      BC(8)= BT4*(R02*CS20-R12*CS21 +BTE*
     *    (R1*SN21 - R0*SN20 + CMC))
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SCOFF(R0,R1,BET1,BET2,BC,AI0,BI0,AI1,BI1,
     *     AI2,BI2,AI3,BI3)
C        off diagonal sines+cosines correction integrals
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 BC(12)
C
      APB = BET1+BET2
      AMB = BET1-BET2
      BMA = BET2-BET1
C
C     if degenerate then SCORR
C
       degen = 1d-6
      IF (ABS(AMB/BET1).LT.degen) THEN
        CALL SCORR(R0,R1,BET1,BC,AI0,BI0,AI1,BI1)
        BC(12)= BC(9)
        BC(11)= BC(8)
        BC(10)= BC(8)
        BC(9) = BC(7)
        BC(8) = BC(6)
        BC(7) = BC(5)
        BC(6) = BC(5)
        BC(5) = BC(4)
        BC(4) = BC(3)
        BC(3) = BC(2)
C        WRITE(KO,*)BC(9)
        RETURN
      ENDIF
C
C
C         A=COS   B=SIN  1,3 = R1   0,2 = R0
C
      SAPB1 = (AI1*BI3 + BI1*AI3)/APB
      SAMB1 = (BI1*AI3 - AI1*BI3)/AMB
      CAPB1 = (AI1*AI3 - BI1*BI3)/APB
      CAMB1 = (AI1*AI3 + BI1*BI3)/AMB
      SBMA1 = SAMB1
      CBMA1 = -CAMB1
C
C             SAPB1 = SIN(APB*R1)/APB
C             SAMB1 = SIN(AMB*R1)/AMB
C             CAPB1 = COS(APB*R1)/APB
C             CAMB1 = COS(AMB*R1)/AMB
C
      SAPB0 = (AI0*BI2 + BI0*AI2)/APB
      SAMB0 = (BI0*AI2 - AI0*BI2)/AMB
      CAPB0 = (AI0*AI2 - BI0*BI2)/APB
      CAMB0 = (AI0*AI2 + BI0*BI2)/AMB
      SBMA0 = SAMB0
      CBMA0 = -CAMB0
C
C             SAPB0 = SIN(APB*R0)/APB
C             SAMB0 = SIN(AMB*R0)/AMB
C             CAPB0 = COS(APB*R0)/APB
C             CAMB0 = COS(AMB*R0)/AMB
C
      SP1 = 0.5D0*(SAPB1+SAMB1)
      SP0 = 0.5D0*(SAPB0+SAMB0)
      SM1 = 0.5D0*(SAMB1-SAPB1)
      SM0 = 0.5D0*(SAMB0-SAPB0)
      CA1 = 0.5D0*(CAPB1+CAMB1)
      CA0 = 0.5D0*(CAPB0+CAMB0)
      CB1 = 0.5D0*(CBMA1+CAPB1)
      CB0 = 0.5D0*(CBMA0+CAPB0)
 
      BC(1) = SP1 - SP0
      BC(2) = CB0 - CB1
      BC(3) = CA0 - CA1
      BC(4) = SM1 - SM0
C
      CP = (CAPB1 - CAPB0)/APB
      CA = (CAMB1 - CAMB0)/AMB
      CB = (CBMA1 - CBMA0)/BMA
      SP = (SAPB1 - SAPB0)/APB
      SA = (SAMB1 - SAMB0)/AMB
      SB = (SBMA1 - SBMA0)/BMA
C
      BC(5) = R1*SP1 - R0*SP0 + 0.5D0*(CA+CP)
      BC(6) = R0*CB0 - R1*CB1 + 0.5D0*(SB+SP)
      BC(7) = R0*CA0 - R1*CA1 + 0.5D0*(SA+SP)
      BC(8) = R1*SM1 - R0*SM0 + 0.5D0*(CA-CP)
C
      CP = CP/APB
      SP = SP/APB
      SA = SA/AMB
      RCAM = (R1*CAMB1-R0*CAMB0)/AMB
      RCAP = (R1*CAPB1-R0*CAPB0)/APB
      RSAP = (R1*SAPB1-R0*SAPB0)/APB
      R02 = R0*R0
      R12 = R1*R1
C
      BC(9) = R12*SP1 - R02*SP0 +RCAM + RCAP -SA -SP
      BC(10)= R02*CB0 - R12*CB1 +RSAP +CB/BMA + CP +
     *   (R1*SBMA1-R0*SBMA0)/BMA
      BC(11)= R02*CA0 - R12*CA1 +RSAP +CA/AMB + CP +
     *   (R1*SAMB1-R0*SAMB0)/AMB
      BC(12) = R12*SM1 - R02*SM0 +RCAM - RCAP -SA +SP
C
      RETURN
      END
C
      SUBROUTINE DAIRY(DX,AI,AIP,BI,BIP)
C  FOR DOUBLE PRECISION ARGUMENTS, THIS ROUTINE CALCULATES THE AIRY
C     FUNCTION AI(X) AND ITS DERIVATIVE AIP(X).  IT ALSO FINDS
C     THE OTHER REAL LINEARLY INDEPENDENT SOLUTION BI(X) AND
C     ITS DERIVATIVE BIP(X).
C     THE DEFINITIONS AND NORMALIZATIONS ARE AS IN NBS HANDBOOK
C     OF MATHEMATICAL FUNCTIONS,P.446
C     THE METHODS USED ARE POWER SERIES EXPANSION FOR SMALL X
C     AND GAUSSIAN INTEGRATION FOR LARGE X
      DIMENSION X(16),W(16),XSQ(16)
      DOUBLE PRECISION DX,AI,AIP,BI,BIP
      DOUBLE PRECISION    XS ,XCUBE,AISUM,AIPSUM
      DOUBLE PRECISION DF,DFP,DG,DGP
      DOUBLE PRECISION FJM2,FJM1,FJ,FJP1,FJP2,FACTOR
      DOUBLE PRECISION C1,C2,ROOT3
      DOUBLE PRECISION DZETA,DARG,DROOTX
      DOUBLE PRECISION ROOT4X,S,CO,RATIO,EFAC,ZETASQ
      DOUBLE PRECISION SUMR,SUMI,SUMRP,SUMIP,TERMR
      DOUBLE PRECISION DZERO,DA,DB,ONE
      DOUBLE PRECISION X,W,XSQ
      DOUBLE PRECISION TEMP,   RTPI,RTPI2
      DOUBLE PRECISION TERMA,TERMB
      LOGICAL NEEDBI
      DATA DZERO,ONE /0.0D0,1.0D0/
      DATA ROOT3/1.732050807568877D0/
      DATA C1,C2 /.355028053887817D0, .258819403792807D0/
      DATA RTPI /.2820947917738781D0/
      DATA RTPI2/.5641895835477562D0/
C *** INSERTED ARB 6.80
C  POSITIONS AND WEIGHTS FOR 10-TERM SUM FOR AIRY FUNCTIONS
       DATA W( 1) /  3.1542515762964787D-14/
       DATA W( 2) /  6.6394210819584921D-11/
       DATA W( 3) /  1.7583889061345669D-08/
       DATA W( 4) /  1.3712392370435815D-06/
       DATA W( 5) /  4.4350966639284350D-05/
       DATA W( 6) /  7.1555010917718255D-04/
       DATA W( 7) /  6.4889566103335381D-03/
       DATA W( 8) /  3.6440415875773282D-02/
       DATA W( 9) /  1.4399792418590999D-01/
       DATA W(10) /  8.1231141336261486D-01/
       DATA X( 1) /  1.4083081072180964D+01/
       DATA X( 2) /  1.0214885479197331D+01/
       DATA X( 3) /  7.4416018450450930D+00/
       DATA X( 4) /  5.3070943061781927D+00/
       DATA X( 5) /  3.6340135029132462D+00/
       DATA X( 6) /  2.3310652303052450D+00/
       DATA X( 7) /  1.3447970824609268D+00/
       DATA X( 8) /  6.4188858369567296D-01/
       DATA X( 9) /  2.0100345998121046D-01/
       DATA X(10) /  8.0594359172052833D-03/
      DATA XSQ( 1) /0.19833317248562170D 03/
      DATA XSQ( 2) /0.10434388535311650D 03/
      DATA XSQ( 3) /0.55377438020178170D 02/
      DATA XSQ( 4) /0.28165249974668990D 02/
      DATA XSQ( 5) /0.13206054139355800D 02/
      DATA XSQ( 6) /0.54338651079380440D 01/
      DATA XSQ( 7) /0.18084791929954200D 01/
      DATA XSQ( 8) /0.41202095387883690D 00/
      DATA XSQ( 9) /0.40402390924418070D-01/
      DATA XSQ(10) /0.64954507303538390D-04/
C  POSITIONS AND WEIGHTS FOR  4-TERM SUM FOR AIRY FUNCTIONS
       DATA W(11) /  4.7763903057577263D-05/
       DATA W(12) /  4.9914306432910959D-03/
       DATA W(13) /  8.6169846993840312D-02/
       DATA W(14) /  9.0879095845981102D-01/
       DATA X(11) /  3.9198329554455091D+00/
       DATA X(12) /  1.6915619004823504D+00/
       DATA X(13) /  5.0275532467263018D-01/
       DATA X(14) /  1.9247060562015692D-02/
      DATA XSQ(11) /0.15365090398596670D 02/
      DATA XSQ(12) /0.28613816631634610D 01/
      DATA XSQ(13) /0.25276291648668180D 00/
      DATA XSQ(14) /0.37044934027789980D-03/
C  POSITIONS AND WEIGHTS FOR  2-TERM SUM FOR AIRY FUNCTIONS
       DATA W(15) /  9.6807280595773604D-01/
       DATA W(16) /  3.1927194042263958D-02/
       DATA X(15) /  3.6800601866153044D-02/
       DATA X(16) /  1.0592469382112378D+00/
      DATA XSQ(15) /0.13542842977111070D-02/
      DATA XSQ(16) /0.11220040761098810D 01/
      IF(DX.LT.-5.0D0) GO TO 100
      NEEDBI=.FALSE.
      IF(DX.GT.3.7D0) GO TO 200
C     THIS ROUTE FOR SMALLX, USING POWER SERIES.
C     INITIALIZE
10    XS  = DX*DX
      XCUBE = XS *DX
      XS  = XS *0.5D0
      DF = C1
      DFP = C1*XS
      DG = C2*DX
      DGP = C2
      AISUM = DF - DG
      AIPSUM = DFP - DGP
      BI = DF + DG
      BIP = DFP + DGP
      FJM2=-2.0D0
20    FJM2=FJM2+3.0D0
      FJM1=FJM2+ONE
      FJ=FJM1+ONE
      FJP1=FJ+ONE
      FJP2=FJP1+ONE
      RATIO = XCUBE/FJ
      DF = DF*RATIO/FJM1
      DFP = DFP*RATIO/FJP2
      DG = DG*RATIO/FJP1
      DGP = DGP*RATIO/FJM2
      BI = BI + (DF+DG)
      BIP = BIP + (DFP+DGP)
      IF(NEEDBI) GO TO 80
      AISUM = AISUM + (DF-DG)
      AIPSUM = AIPSUM + (DFP-DGP)
C     CONVERGENCE TEST
80    IF(DABS(DF).GT.1.0D-16) GO TO 20
C     CONVERGENCE. COMPUTE FUNCTIONS
      BI = ROOT3*BI
      BIP = ROOT3*BIP
C  THIS RETURNS IF X IS BETWEEN 3.7 AND 8.0, SINCE IN SUCH CASES MORE
C  ACCURATE VALUES OF AI AND AIP HAVE ALREADY BEEN FOUND BY GAUSSIAN
C  INTEGRATION
      IF(NEEDBI)RETURN
      AI = AISUM
      AIP = AIPSUM
      RETURN
C  GAUSSIAN INTEGRATION FOR LARGE NEGATIVE X
100   DROOTX = DSQRT(-DX)
      ROOT4X = DSQRT(DROOTX)
      DZETA = -.6666666666666667D0*DX*DROOTX
      DARG = DZETA - .7853981633974483D0
      SUMR = DZERO
      SUMI = DZERO
      SUMRP = DZERO
      SUMIP = DZERO
C  TEST TO SEE HOW MANY TERMS ARE NEEDED IN GAUSSIAN INTEGRATION
      IF(DX.LT.(-200.D0)) GO TO 140
      IF(DX.LT.(-15.D0)) GO TO 130
C  THIS CASE FOR DX BETWEEN -5.0 AND -15.0
      LIMLO=1
      LIMHI=10
      GO TO 149
C  THIS CASE FOR DX BETWEEN -15.0 AND -200.
130   LIMLO=11
      LIMHI=14
      GO TO 149
C  THIS CASE FOR DX.LT.-200.
140   LIMLO=15
      LIMHI=16
149   ZETASQ=DZETA**2
      DO 150 K=LIMLO,LIMHI
      TERMR=W(K)/((ZETASQ+XSQ(K))**2)
      SUMR = SUMR + TERMR
      TERMR=TERMR*X(K)
      SUMI=SUMI+TERMR
      TERMR=TERMR*X(K)
      SUMRP=SUMRP+TERMR
150   SUMIP=SUMIP+TERMR*X(K)
      SUMR=(SUMR*ZETASQ+SUMRP)*ZETASQ
      TEMP=SUMI*ZETASQ
      SUMI=(TEMP+SUMIP)*DZETA
      SUMRP=SUMRP*DZETA
      SUMIP=SUMIP-TEMP
C  FORM AIRY FUNCTIONS
      S = DSIN(DARG)
      CO = DCOS(DARG)
      RATIO = RTPI2/ROOT4X
      AI = RATIO*(CO*SUMR + S*SUMI)
      BI = RATIO*(CO*SUMI - S*SUMR)
      SUMRP=SUMRP+SUMRP
      RATIO = -.25D0/DX
      FACTOR = -RTPI2*ROOT4X
      AIP = RATIO*AI - DROOTX*BI + FACTOR*(CO*SUMRP+S*SUMIP)
      BIP = RATIO*BI + DROOTX*AI + FACTOR*(CO*SUMIP-S*SUMRP)
      RETURN
C   GAUSSIAN INTEGRATION FOR LARGE POSITIVE X
200   DROOTX = DSQRT(DX)
      DZETA = .6666666666666667D0*DX*DROOTX
      EFAC = DEXP(-DZETA)
      ROOT4X = DSQRT(DROOTX)
      AI = DZERO
      BI = DZERO
      AIP = DZERO
      BIP = DZERO
      IF(DX.LT.8.0D0) NEEDBI=.TRUE.
C  TEST TO SEE HOW MANY TERMS ARE NEEDED IN GAUSSIAN INTEGRATION
      IF(DX.GT.15.0D0) GO TO 230
C  THIS CASE FOR DX BETWEEN 3.7 AND 15.
      LIMLO=1
      LIMHI=10
      GO TO 249
C  THIS CASE FOR DX GREATER THAN 15.
230   LIMLO=11
      LIMHI=14
249   DO 250 K=LIMLO,LIMHI
      DA=DZETA+X(K)
      TERMA = W(K)/DA
      AI = AI + TERMA
      AIP=AIP+TERMA*X(K)/DA
      IF(NEEDBI) GO TO 250
      DB=DZETA-X(K)
      TERMB = W(K)/DB
      BI = BI + TERMB
      BIP=BIP+TERMB*X(K)/DB
250   CONTINUE
C  FORM FUNCTIONS
      FACTOR=RTPI*DZETA/ROOT4X
      RATIO = 0.25D0/DX
      AI=AI*EFAC*FACTOR
      AIP=-(DROOTX+RATIO)*AI+RTPI*ROOT4X*EFAC*AIP
C  THIS IS SATISFIED ONLY FOR X BETWEEN 3.7 AND 8.0  IN THESE CASES
C  THE BI AND BIP ABOUT TO BE COMPUTED ARE NOT SUFFICIENTLY ACCURATE.
C  THUS RETURN TO POWER SERIES FOR BI AND BIP.
      IF(NEEDBI) GO TO 10
      FACTOR=FACTOR+FACTOR
      BI=BI*FACTOR/EFAC
      BIP=(DROOTX-RATIO)*BI-RTPI2*ROOT4X*BIP/EFAC
      RETURN
      END

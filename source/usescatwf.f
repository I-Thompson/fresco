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

       TIME0J = TM(I)
       if(SMATL.ge.3.and.NCHAN.ge.3) then
       write(KO,1433) 
 1433       format(' With velocity factors, final S-matrix elements:')
!	  rtrace = 0.
	  XSC = 0.
          DO 540 C=1,NCH
            IC = PART(C,1)
            IA = EXCIT(C,1)
!@@
            IF(RMASS(IC).lt.1e-5) then
                T = sqrt(RMASS(PEL)*amu/ (HBC*K(PEL,EXL)))
            ELSE IF(RMASS(PEL).lt.1e-5) then
                T = 1./sqrt(RMASS(IC)*amu/ (HBC*K(IC,IA)))
            ELSE
                T = sqrt(K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL)))
            ENDIF
!	    see justification at frxx3.f
!
!           T = sqrt(K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL)))
!@@
          SOMA(C,1) = SMAT(C) * T
!	  if(C.ne.initl(JIN).and.abs(SMAT(C))>1d-50) rtrace = abs(SMAT(C))
	  T4 = abs(SOMA(C,1))
	  if(C.ne.initl(JIN).and.T4>1d-50) XSC = T4
540       CONTINUE
             WRITE(KO,1430) (SOMA(C,1),C=1,NCH)
!	  if(abs(rtrace)>1d-50) write(155+JIN,*) ECM(initl(JIN),1),rtrace
!	  if(abs(XSC)>1d-50) write(165+JIN,*) ECM(initl(JIN),1),XSC
       endif
C
      if(NFUS>0) FUSL(2:1+NFUS) = 0.0
C-------------------------------------------CORE FUSION POTENTIAL
      DO 556 IT=1,NFUS
      DO 556 C1=1,NCH
      if(ITC(PART(C1,1),EXCIT(C1,1))==IT) then
        DO 555 JF=1,NF
          IF(PTYPE(1,JF).NE.KFUS) GO TO 555
          IF(PTYPE(2,JF).EQ.0)    GO TO 555
          IF(PTYPE(2,JF).GT.2)    GO TO 555
          IF(PTYPE(3,JF).NE.0)    GO TO 555
          T4 = 0.0
          DO 554 I=1,N
554       T4 = T4 + AIMAG(FORMF(I,JF)) * ABS(PSI(I,C1))**2
          T4 = T4 * HP(PEL) / (K(PEL,EXL) * ECM(EL,1))
          FUSL(1+IT) = FUSL(1+IT)-T4 * XCOEF * (2*JTOTAL+1.) * 4.*PI
555       CONTINUE
      ENDIF
556   continue
C
      WRITE(KS) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1)
C
      IF(VEFF.NE.0) THEN
C-------------------------------------------EFFECTIVE POTENTIALS
            DO 560 I=1,NLN
  560       VPOT(I,1) = PSI((I-1)*MR+1,EL)
            CALL CHECK(IEX,MAXICH,25)
        if(NL>0.and.NRBASES/=0) write(6,*)
     x          ' NONLOCAL COUPLINGS NOT INCLUDED IN VEFF!!'
!       write(48,*) 'N,(N-1)*HP(1) =',N,(N-1)*HP(1)
!       write(48,*) 'NLN1,(NLN1-1)*HP(1)*MR =',NLN1,(NLN1-1)*HP(1)*MR
!       write(48,*) 'NLN,(NLN-1)*HP(1)*MR =',NLN,(NLN-1)*HP(1)*MR
        write(48,*) 'IEX =',IEX
        DO 590 C2=1,IEX
        DO 590 NC=1,NCLIST(EL,C2)
         IF = NFLIST(EL,C2,NC)
!           write(48,*) 'EL,C2,NC =',EL,C2,NC,' IF =',IF
C     IF(EL.GT.IEX .OR. PTYPE(3,IF).LE.0) GO TO 590
      IF(.NOT.(EL.LE.IEX .AND.
     X   (PTYPE(3,IF).GT.0.OR.PTYPE(3,IF).EQ.0.AND.EL.NE.C2))) GOTO 590
           C6 = CLIST(EL,C2,NC)
!           write(48,*) 'EL,C2,NC =',EL,C2,NC,' S=',C6, ' IF =',IF
           IF(ABS(C6).LT.1E-10) GO TO 590
           DO 580 I=1,N
  580      SRC(I,1) = SRC(I,1) + C6 * FORMF(I,IF) * PSI(I,C2)
  590   CONTINUE
      REWIND 15
C   Choose T as J-dependent factor for Weighted Equivalent Potential
         T = (2*JTOTAL+1.) * (1. - ABS(SMAT(EL))**2)
         IF(ABS(VEFF).EQ.2 .AND. ABS(SMAT(EL)) .LT. 0.1) T = 0.0
         IF(ABS(VEFF).EQ.3) THEN
          T = 0.0
          DO 595 C=1,NCH
  595       IF(C.NE.EL) T = T + ABS(SMAT(C))**2
          T = T * (2.*JTOTAL + 1.)
         ENDIF
         IF(ABS(VEFF).EQ.4) T = 1.0
           write(48,*) 'Weight factor =',T
C
          CALL CHECK(NSA,MAXCH,5)
         DO IMA=1,NSA
           READ(15) (FNC(I,IMA),I=1,NLN1)
         ENDDO
         DO IMA=1,NSA
           READ(15) (WNM(I,IMA,2),I=1,NLN1)
         ENDDO
         IMA = NINT(JVAL(EL) - LVAL(EL) + JEX(1,PEL,EXL)) + 1
         DO I=1,NLN
          J = (I-1)*MR + 1
          FNC(I,IMA) = FNC(I,IMA) + CONJG(PSI(J,EL))*SRC(J,1) * T
          WNM(I,IMA,2) = WNM(I,IMA,2) +   ABS(VPOT(I,1))**2     * T
         ENDDO
      REWIND 15
         DO IMA=1,NSA
          WRITE(15) (FNC(I,IMA),I=1,NLN1)
         ENDDO
         DO IMA=1,NSA
          WRITE(15) (WNM(I,IMA,2),I=1,NLN1)
         ENDDO
      ENDIF
C
      DO 655 C=1,NICH
         IC = PART(C,1)
         IA = EXCIT(C,1)
         IT = ITC(IC,IA)
      AMDSQS = ABS(SMAT(C))**2
      IF(mod(WAVES,2).ne.0) THEN
         DO 652 I=CUTVAL(C),N
            VPOT(I,1) =PSI(I,C)
            IF(abs(WAVES).GE.4) THEN
                C7 = PSI(I,EL) + 1E-30
                if(WAVES.le.0) C7 = VPOT(I,3) + 1E-30
                GO TO  651
                ENDIF
            IF(WAVES.GE.0) GO TO 652
        CALL COULFG((I-1)*ECM(C,3)*K(IC,IA),ETA(IC,IA),Z,LVAL(C)+1D0,
     X                  CFG,CFG(1,2),CFG(1,3),CFG(1,4),2,0,J,L1)
            IF(J.GT.0) THEN
               VPOT(I,1) = 0.
               GO TO 652
            ENDIF
            L1 = LVAL(C)+1
            T = 0.0
            IF(C.EQ.EL) T = 1.0
               CF = CFG(L1,1)
               CG = CFG(L1,2)
            C6 = (CG + CI*CF) * (0.,0.5)
            C7 = (CG - CI*CF) * (0.,-.5)
            C7 = -(T* C7 + SMAT(C)*C6)
C           C7 = CF
 
  651       VPOT(I,1) = VPOT(I,1)/C7
            if(C.eq.1) VPOT(I,2) = 0.0
            VPOT(I,2) = VPOT(I,2) + VPOT(I,1)
  652     CONTINUE
       ENDIF
      IF(mod(WAVES,2).NE.0.AND.ABS(AMDSQS).GT.1E-16) THEN
      IF(abs(WAVES).le.3) THEN
        WRITE(KO,1440) C,SMAT(C),AMDSQS,(VPOT(I,1),I=CUTVAL(C),M)
 1440 FORMAT(//' For channel no.',I3,', the S-matrix element is',
     X F15.5,' +i*',F10.5,   ' (square modulus =',F15.8,')',/
     X ' with the radial wave function'//
C    @   (5(F12.5,' +i*',F9.5,',') ) )
     X      1P,            (5(E12.3,' +i*',E9.2,',') ) )
       ELSE
        IF(C.NE.EL.or.WAVES.le.0)
     X    WRITE(KO,1441) C,SMAT(C),AMDSQS,(C,ENLAB,(I-1)*HP(IC),
     X LVAL(C),ABS(VPOT(I,2))**2,VPOT(I,1),I=MR+1,M,MR)
 1441 FORMAT(//' For channel no.',I3,', the S-matrix element is',
     X F15.5,' +i*',F10.5,   ' (square modulus =',F15.8,')',/
     X ' with ratio to elastic wavefunction:'//
C    @   (5(F12.5,' +i*',F9.5,',') ) )
     &   (' MX:',I2,2F6.2,I4,2F12.5,' +i*',F9.5) )
      ENDIF
      ENDIF
  655 CONTINUE
C
      IF(WDISK.NE.0) THEN
	  call openif(17); written(17) = .true.
	  NA = (N-1)/2*2 + 1 ! print these
      IF(WDISK>0) WRITE(17,657) NA,HP(PEL),ENLAB,JTOTAL,PARITY,
     x               MASS(:,PEL)
      IF(WDISK<0) WRITE(17) NA,HP(PEL),ENLAB,JTOTAL,PARITY,MASS(:,PEL)
657   FORMAT(I4,2F8.4,F8.1,I3,2f12.6,2f8.3)
      DO 690 IT=1,ITCM
          IF(ABS(MOD(WDISK,2)).EQ.1.AND.IT.NE.ITC(PEL,EXL)) GOTO 690
	  EL = initl(JIN)
          DO 670 C=1,NICH
            IF(ITC(PART(C,1),EXCIT(C,1)).NE.IT) GOTO 670
            IF(ABS(MOD(WDISK,2)).EQ.1.AND.C.ne.EL) GOTO 670
            IF(WDISK.GT.0)
     X     WRITE(17,660) IT,LVAL(C),JVAL(C),JTOTAL,LVAL(EL),JVAL(EL),
     X                   SMAT(C),ETA(PART(C,1),EXCIT(C,1))
            IF(WDISK.LT.0)
     X     WRITE(17) IT,LVAL(C),JVAL(C),JTOTAL,LVAL(EL),JVAL(EL),
     X                   SMAT(C),ETA(PART(C,1),EXCIT(C,1))
  660         FORMAT(2I4,2F6.1,I4,F6.1,2F15.10,f12.8)
           IF(WDISK.GT.0.AND.WDISK.LE.2) THEN
              WRITE(17,'(1p,6e12.4)') (PSI(I,C),I=1,NA)
           ELSE IF(WDISK.GE.3) THEN
              WRITE(17,681) (PSI(I,C),I=1,NA)
           ELSE
             DO 665 I=2,NA
             C6 = PSI(I,C)
  665        WRITE(17) C6
           ENDIF
  670      CONTINUE
  681      FORMAT(1P,6E13.6)
  690  CONTINUE
        IF(WDISK.GT.0) WRITE(17,660) -1,0,0.,0.,0,0.,0.,0.,0.
        IF(WDISK.LT.0) WRITE(17) -1,0,0.,0.,0,0.,0.,0.  ,0.
       ENDIF
C
C     PRINT OUT REACTION CROSS SECTIONS
C     ---------------------------------
      DO 692 IT=0,ITCM
 692  SIGJ(IT) = 0.0
      CALL DISPX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGJ,PART,EXCIT,
     X           JEX,ITC,K,RMASS,PEL,EXL,LVAL,FUSL,OUTJ,
     X           JVAL,JPROJ,JTARG,JTOTAL,PARITY,
     X           SMATS,CHANS+DONE,JTOTAL.GE.JTMAX-3.,SCALE,
     X           ITCM,IF1,IF2)
      if(IAME==0) then ! do it now!
       FUSJ = FUSJ + FUSL(1)
       TOTJ = TOTJ + OUTJ
       FUSJPI = FUSJPI + FUSL(1)
       if(NFUS>0) CFUSJ(1:NFUS) = CFUSJ(1:NFUS) + FUSL(2:1+NFUS)
      endif
      IF(REAL(SMAT(EL))>0.9) then
      SMALLJ =  SMALLCOUP>0.
      DO 693 IC=1,NCHAN
        NA = NEX(IC)
        DO 693 IA=1,NEX(IC)
        IT = ITC(IC,IA)
      CHSIZES(IT) = SIGJ(IT) / JCOEF
      IF(IT/=ITC(PEL,EXL).and.K(IC,IA)>1e-5) then
      IF(CHSIZES(IT)<SMALLCHAN.and.CHPRES(IT)>0) 
     X		SMALLS(IT) = SMALLS(IT) + 1
      IF(CHSIZES(IT)>SMALLCOUP) SMALLJ = .false.
        T = ETA(IC,IA)
        T4 =(T+SQRT(T**2 + JTOTAL*(JTOTAL+1d0)))/K(IC,IA)
*	write(139,*) ic,ia,jtotal,real(rturn),real(rturn+gap),real(t4)
*	call flush(139)
!		Drop a channel if its  turning point is too far outside
!		the turning point of the elastic channel.
	if(T4>RTURN+GAP.and.SMALLCHAN>0.) SMALLS(IT) = SMALLS(IT) + 1
      endif
693   CONTINUE
!	write(138,'(f7.1)')  JTOTAL
!	write(138,'(1p,8e10.2)')  CHSIZES(1:ITCM)
!	call flush(138)
      ELSE
       SMALLJ =  .false.
      ENDIF

C
      IF(.not.FAIL) then
      IF(SMATS.EQ.3.AND.ABS(DBLE(SMAT(EL)) - 0.5).LT.0.45) SMATL = 4
      IF(SMATS.EQ.3.AND.ABS(DBLE(SMAT(EL)) - 0.5).GT.0.45) SMATL = 2
      ENDIF
C
      IF(XS(EL).LT.ABSEND.AND.JTOTAL.GT.ABS(JTMIN).and.
     x   (IEXCH==0.or.XS(EL)>1e-5).and.JTOTAL>jleast) DONES = DONES + 1
       DO 695 I=1,NTHRESH
  695  IF(DBLE(SMAT(EL)).gt.RESM(I).and.THRJ(I).le..1) THRJ(I)=JTOTAL

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

****BPMFUS**************************************************************
!	SUBROUTINE BPMFUS
!      	implicit real*8(a-h,o-z)
***        COMPLEX*16 FFC
C                      Calculate fusion cross sections using BPM &WKB::
C                     (Ignoring any spins of projectile and/or target,
C                      by using only scalar parts of optical potential,
C                      and using the spin-average of the equivalent
C                          local potential FNC from just above)
      	LL = JTOTAL+0.5
         VPOT(1:N,1) = 0.0
         DO 860 IF=1,NF
	  I = PTYPE(2,IF)
          if((I.le.2.or.(I.ge.10.and.I.le.13) )
     X		.and. PTYPE(3,IF).eq.0) then
         VPOT(1:N,1) = VPOT(1:N,1)+ELCOEF(IF)*FORMF(1:N,IF)
	  endif
  860   CONTINUE
C***********************************************************************
C    Use this R1 to distinguish interior from exterior regions
C                 (full fusion & none, respectively)
C***********************************************************************
      R1 = 1.0 * (MASS(1,PEL)**0.333 + MASS(2,PEL)**0.333)
C***********************************************************************
      C = 1
      IF(VEFF.NE.0) C=2
      DO 940 C1=1,C
        DO 861 I=1,3
  861   CORFUS(I,1) = 0.0
      IF(C1.EQ.2) THEN
C                     Case 2: add local equivalent potential (VOPT).
        DO 862 I=1,M
         R = (I-1)*HP(PEL) / (HP(PEL)*MR)
         S = FFC(R,VOPT,NLN)
         VPOT(I,2) = VPOT(I,1) + S
862      VPOT(I,1) = S
C     Print out Polarisation Potential for latter use:
         IF(BPM.LT.0) WRITE(7,680) (VPOT(I,1),I=1,M)
  680      FORMAT(1P,6E12.4)
         IF(BPM.LT.0) written(7) = .true.
       ENDIF
        WRITE(KO,863)
863     FORMAT(/' BPM: Using the real scalar part of the ',
     X          ' Elastic Channel Potential')
        IF(C1.EQ.2) WRITE(KO,*) 'Including CRC local equivalent potentia
     Xl.'
        WRITE(KO,864) ECMC(PEL,EXL),R1
864     FORMAT(/' Incident energy (cm) =',F8.3,'. Fusion inside radius o
     Xf',F8.3,' fm.'/)
        IMA = R1/HP(PEL) + 1.5
C        IF(ABS(BPM).GE.3)  WRITE(6,680) (VPOT(I,C1),I=1,M)
        DO 900 L=0,LL
           T = L*(L+1) / (FMSCAL * RMASS(PEL))
           JJ = 0
           EKL1 = 0.0
           RTURN1 = 0.0
           RTURN2 = 0.0
           W = 'Turn'
           BAR = -1E4
           POS = 0.0
           RSO = 0.0
           DO 870 I=M,IMA,-1
            R = (I-1)*HP(PEL)
            S = VPOT(I,C1) + T/R**2 - ECMC(PEL,EXL)
            RS = DBLE(S)
      		IF(L.le.-11.and.BPM>1) 
     x	 WRITE(KO,*) 'I,RS,RSO,RTURN2 =',I,RS,RSO,RTURN2,EKL1
C   RS is positive during sub-barrier penetration.
            IF(I.LT.M) THEN
              IF(RS*RSO.LT.0.0) THEN
                JJ = JJ+1
                IF(JJ.EQ.1) RTURN2 = R
                IF(JJ.EQ.2) RTURN1 = R
!		write(6,*) ' Found t.p. no. ',JJ,' @ ',R,' for L =',L
		call flush(6)
                IF(JJ.GT.2) GO TO 871
              ENDIF
            IF(RS.GT.0.) EKL1 = EKL1 + SQRT(ABS(RS))
            ENDIF
            IF(BAR.LT.RS) THEN
               BAR = RS
               DER = (RS - 2.*RSO + RSOO) / HP(PEL)**2
               POS = R
             ENDIF
            RSOO = RSO
            RSO  = RS
 870      CONTINUE
 871  EKL = EKL1 * HP(PEL) * SQRT(FMSCAL*RMASS(PEL))
      IF(BAR.LE.0.) THEN
         HOM = SQRT(ABS(DER)*2./(FMSCAL*RMASS(PEL)))
C
C    I DO NOT KNOW why HOM has to be multiplied by 2 in the next line!!!
       HOM = HOM*2
C     IF(BAR.GT.0) WRITE(KO,*) 'EKL : WKB = ',EKL,', H-W = ',PI*BAR/HOM
          RTURN1 = POS
          RTURN2 = HOM
          W = 'Peak'
         EKL = PI * BAR / HOM
       ENDIF
      TLE = 0.0
      if(EKL.lt.10.) TLE = 1.0/(1 + EXP(2.0*EKL))
      SIGL = 10.0 * PI/ABS(K(PEL,EXL))**2 * (2*L+1d0) * TLE
      IF(ABS(BPM).GE.2.AND.TLE.GT.1E-5)
     x      WRITE(KO,890) L,HOM,BAR,W,RTURN1,RTURN2,TLE,SIGL
 890  FORMAT(' L =',I3,', Barrier =',2F6.2,', ',A4,' @',2F6.2,', Tl =',
     X  F7.4,': Xsec-F =',F8.3,' mb')
      CORFUS(1,1) = CORFUS(1,1) + SIGL
      CORFUS(2,1) = CORFUS(2,1) + SIGL * L
      CORFUS(3,1) = CORFUS(3,1) + SIGL * L**2
 900  CONTINUE
      DO 901 I=2,3
 901  CORFUS(I,1) = CORFUS(I,1)/(CORFUS(1,1) + 1E-20)
      WRITE(KO,902) C1,CORFUS(:,1)
 902  FORMAT('0Cumulative FUSION cross section by BPM (case',I2,')',
     X       '   =',E11.3,'  <L> =',F9.2,'  <L**2> =',F9.1)
      IF(C1.EQ.1) BARE = CORFUS(1,1)
  940 CONTINUE
      IF(C.EQ.1) WRITE(KO,945) ENLAB,BARE,TOTFUS(1)
      IF(C.EQ.2) WRITE(KO,945) ENLAB,BARE,TOTFUS(1),CORFUS(1,1)
  945 FORMAT('0Fusion at',F8.3,' MeV: Bare =',1P,E11.3,', Actual =',
     X E11.3,:,', BPM-equiv =',E11.3,' mb')
      write(391,'(G12.4,1p,E12.5)') ENLAB,BARE
!	END SUBROUTINE BPMFUS

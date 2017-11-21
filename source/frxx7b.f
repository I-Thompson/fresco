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

       subroutine readpt(ki,namep,massp,zp,nex,pwf,namet,masst,zt,qval)
       implicit real*8(a-h,o-z)
	character*8 namep,namet
	real*8 massp,masst
	logical pwf
	integer readstates
	namelist/partition/namep,massp,zp,nex,pwf,namet,masst,zt,qval,
     x                      readstates
	namep='        ';massp=0; zp=0; nex=0; pwf=.false.
	namet='        ';masst=0; zt=0; qval=0
	read(ki,nml=partition)
	return
	end
       subroutine readst(ki,jp,copyp,bandp,ep,kkp,tp,cpot,
     X		jt,copyt,bandt,et,kkt,tt,fexch,ignore,infam,outfam)
       implicit real*8(a-h,o-z)
	real*8 jp,kkp,jt,kkt
	integer copyp,bandp,copyt,bandt,cpot,ptyp,ptyt,infam,outfam
	logical keep,fexch,ignore
	namelist/states/ jp,  copyp,bandp,ptyp,ep,  kkp,tp, cpot,
     X		jt,copyt,bandt,ptyt,et,kkt,tt,  fexch,ignore,keep,
     X          infam,outfam
	jp=0;   copyp=0; bandp =0; ep=0;   kkp=0; tp=0;  cpot=0; 
     	jt=0;   copyt=0; bandt =0; et=0;   kkt=0; tt=0;
	fexch=.false.; ignore=.false.; infam=0; outfam=0
	read(ki,nml=states)
	return
	end

      SUBROUTINE FFNL(FNC,NLL,NLO,SCALE,LTR,PTR,TTR, ICTO,ICFROM,IB,IA)
      use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 FNC(NLL,NLO),CSCAL,SCALE
C             EXPONENTIAL FORM FACTOR, range RANGE, with R**FK factor,
C           and GAUSSIAN NON-LOCALITY, range B
      READ (4,*) RANGE,FK,B,CSCAL
      WRITE(KO,10) RANGE,FK,B,CSCAL
10    FORMAT('0Constructing NON-LOCAL FORM FACTOR, an EXPONENTIAL with '
     $ ,' RANGE ',F8.4,' fm., R**',F4.2,' factor, ',/,
     $' and NON-LOCAL WIDTH ',F8.4,' fm.  Complex scale factor =',2F8.4)
      DO 20 I=1,NLL
      DO 20 J=1,NLO
       R = REAL(FNC(I,J))
       F = 0.0
       IF(R.GT.0) F = R ** FK * EXP(-R/RANGE) * SCALE * CSCAL
       IF(I.GE.NLL-2) F = 0.0
        RF = AIMAG(FNC(I,J))
        DNL= RF - R
20      FNC(I,J) = F * EXP(-(DNL/B)**2)
      RETURN
      END
      SUBROUTINE FFNL2(FNC,NLL,NLO,SCALE,LTR,PTR,TTR, ICTO,ICFROM,IB,IA)
      use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 FNC(NLL,NLO),CSCAL,F,SCALE
C        WS FORM FACTOR,
C         Strength CSCAL, radius R0*A**1/3, diff. AV, with R**FK factor,
C            and GAUSSIAN NON-LOCALITY, range B
      READ (4,*) CSCAL,R0,A,AV,FK,B
      WRITE(KO,10) CSCAL,R0,A,AV,FK,B
10    FORMAT('0Constructing NON-LOCAL FORM FACTOR, a WOODS-SAXON of '
     $ ,'size ',2F8.4,', r0 =',F8.4,' (A=',F6.2,'), AV =',F7.4,
     $ /,'0with R**',F4.2,' factor and NON-LOCAL WIDTH ',F8.4,' fm.')
      RAD = R0 * A**0.3333333
      CSCAL = CSCAL * SCALE / (SQRT(3.14159)*B)
      DO 20 I=1,NLL
      DO 20 J=1,NLO
       R = REAL(FNC(I,J))
       F = 0.0
       E = 1./(1. + EXP((R-RAD)/AV))
       IF(R.GT.0) F = R ** FK * E * SCALE * CSCAL
       IF(I.GE.NLL-2) F = 0.0
        RF = AIMAG(FNC(I,J))
        DNL= RF - R
20      FNC(I,J) = F * EXP(-(DNL/B)**2)
      RETURN
      END
      SUBROUTINE FNLSET(FNC,NLL,NLO,HTARG,DNL,
     X                  SCALE,LTR,PTR,TTR, ICTO,ICFROM,IB,IA)
      use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 FNC(NLL,NLO),CSCAL,SCALE
      REAL*8 DNL(NLO)
C             EXPONENTIAL FORM FACTOR, range RANGE, with R**FK factor,
C           and GAUSSIAN NON-LOCALITY, range B
      READ (4,*) RANGE,FK,B,CSCAL
      WRITE(KO,10) RANGE,FK,B,CSCAL
10    FORMAT('0Constructing NON-LOCAL FORM FACTOR, an EXPONENTIAL with '
     $ ,' RANGE ',F8.4,' fm., R**',F4.2,' factor, ',/,
     $' and NON-LOCAL WIDTH ',F8.4,' fm.  Complex scale factor =',2F8.4)
      DO 20 I=1,NLL
       R = (I-1)*HTARG
       F = 0.0
       IF(R.GT.0) F = R ** FK * EXP(-R/RANGE) * SCALE * CSCAL
       IF(I.GE.NLL-1) F = 0.0
      DO 20 J=1,NLO
        RP = R + DNL(J)
        IF(RP.LE.0.0) GO TO 20
        FNC(I,J) = F * EXP(-(DNL(J)/B)**2)
20      CONTINUE
      RETURN
      END

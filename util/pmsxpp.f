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

       implicit real*8(a-h,o-z)
       dimension XSEC(7,7),T0(0:6),PM(0:6)
       call logfac(20)
        read(5,*) KQ1PR
C       write(6,*) KQ1PR
       j = (KQ1PR-1)/2
       xj = j
       z  = 0d0
      write(0,*) 'M states for spin =',j
      a = 1.0/SQRT(KQ1PR*1.)
       do 2 m=0,j
       do 2 k=0,j*2,2
2     write(0,*)' Coefficient of t(',k,'0) in p(',m,') is ',
     X     a*(-1)**(j-m) * cleb6(xj,m+z,xj,-m+z,k+z,z)
1     read(5,*,err=20,end=20) th,xs,
     X     ((XSEC(KQ1,LQ1),LQ1=2-MOD(KQ1,2),KQ1),KQ1=2,KQ1PR)
C     write(6,*) th,xs,
C    X     ((XSEC(KQ1,LQ1),LQ1=2-MOD(KQ1,2),KQ1),KQ1=2,KQ1PR)
       do 5 k=2,j*2,2
 5     t0(k) = XSEC(k+1,1)
       t0(0) = 1.
      a = 1.0/SQRT(KQ1PR*1.)
       do 10 m=0,j
         if(m.eq.1) a=2.*a
       pm(m) = 0.0
       do 10 k=0,j*2,2
10     pm(m) = pm(m)+a*(-1)**(j-m) * cleb6(xj,m+z,xj,-m+z,k+z,z)*t0(k)
      print *,real(th),(real(pm(m)),m=0,j)
      go to 1
  
20    STOP
      end
      SUBROUTINE LOGFAC(L)
C  FACT(I)= LN((I-1)!)
      REAL*8 FACT
      COMMON/LFAC/FACT(101)
C  THE VALUE OF L CORRESPONDS TO THE DIMENSION OF (FACT(
      FACT(1)=0.
      DO 100 J=2,L
  100 FACT(J)=FACT(J-1)+LOG(J-1D0)
      RETURN
      END
      FUNCTION CLEB6(A,AL,B,BE,C,M)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL* 8 M,FACT
      COMMON / LFAC / FACT(11)
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      GA = -M
C
      IF(AL+BE+GA) 11,10,11
C11    WIG3J = 0.0
11    CLEB6 = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GO TO 11
      IF(A-ABS(AL)) 11,14,14
14    IF(B-ABS(BE)) 11,15,15
15    IF(C-ABS(GA)) 11,13,13
13    IA = C-B+AL
      IB = C-A-BE
      IF(IA) 20,21,21
21    IF(IB) 20,23,23
23    MIN = 0
      GO TO 24
20    MIN = -IA
      IF(MIN+IB) 25,24,24
25    MIN = -IB
24    IC = A-AL
      ID = B+BE
      IE = A+B-C
      NIN = MIN
      T=PHASE(MIN)
      S = T
30    MIN = MIN+1
      IZ = IC-MIN+1
      IF(IZ) 29,29,26
26    IY = ID-MIN+1
      IF(IY) 29,29,27
27    IX = IE-MIN+1
      IF(IX) 29,29,28
28    TA = IX*IY*IZ
      TB = MIN*(IA+MIN)*(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    I = B-A+GA
      IF(S.EQ.0.0) GO TO 11
      IH = A+AL
      II = B-BE
      IJ = C+GA
      IK = C-GA
      IL = A+C-B
      IM = B+C-A
      IN = A+B+C+1.0
      XDDD = 0.5*(FACT(IH+1)+FACT(IC+1)+FACT(ID+1)
     $           +FACT(II+1)+FACT(IJ+1)+FACT(IK+1)
     $           +FACT(IE+1)+FACT(IM+1)+FACT(IL+1)-FACT(IN+1))
     $      - ( FACT(IC-NIN+1)+FACT(IA+NIN+1)+FACT(ID-NIN+1)+
     $          FACT(IB+NIN+1)+FACT(NIN+1)+FACT(IE-NIN+1))
C     WIG3J = (-1.0)**I *  EXP(XDDD) * S
      CLEB6 = SQRT(2*C+1.)*EXP(XDDD) * S
      RETURN
      END

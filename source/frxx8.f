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

*****FRXX8*************************************************************
      SUBROUTINE COULPO(C,RC,AC,H,T,NM)
      REAL*8 C(NM),RC,AC,H,T,I1,U,B,Z0,Z1,Z2,Y0,Y1,Y2,A
      INTEGER J,N
      IF(RC.NE.0)GO TO 2
      I1=T/H
      DO 10 N=1,NM
   10 C(N)=I1/N
      RETURN
    2 IF(AC.NE.0)GO TO 1
      J=RC/H+1.0E-10
      I1=T/H
      IF(J)11,11,12
   12 DO 13 N=1,J
   13 C(N)=T*(3-(N*H/RC)**2)/(2*RC)
   11 IF(J.GE.NM)GO TO 14
      JJJ=J+1
      DO 15 N=JJJ,NM
   15 C(N)=I1/N
   14 RETURN
    1 B=EXP(RC/AC)
      U=EXP(-H/AC)
      J=(20.*AC+RC)/H+1.0E-10
      IF(J.LT.NM)J=NM
      Y1=0.
      Y2=1.5*H/RC
      C(1)=Y2
      B=B*U
      Z0=0.
      Z1=B/(1.+B)
      DO 16 N=2,J
      B=B*U
      Z2=N*B/(1+B)
      Y0=Y1
      Y1=Y2
      Y2=2*Y1-Y0+10*Z1+Z0+Z2
      IF(NM.GE.N)C(N)=Y2
      Z0=Z1
   16 Z1=Z2
      A=((Y1-Y0+Y1-Y2)*(J-0.5)+Y2-Y1)/H
      B=Y2-A*J*H
      A=A*T/B
      B=T/(B*H)
      DO 17 N=1,NM
   17 C(N)=B*C(N)/N-A
      RETURN
      END
      FUNCTION FF(Y,HI,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,Y,HI,P1,P2,Q,X,FF
      DATA X/.16666666666667D0/
      P=Y*HI
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FF=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FF=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FF=F(N)
      RETURN
      END
      FUNCTION FF4(Y,HI,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,Y,HI,P1,P2,Q,X,FF4
      DATA X/.16666666666667E0/
      P=Y*HI
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FF4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FF4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FF4=F(N)
      RETURN
      END
      FUNCTION FFC4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 F(N),FFC4
      REAL*8 P,Y,P1,P2,Q,X
      DATA X/.16666666666667E0/
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC4=F(N)
      RETURN
      END
      FUNCTION FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      DATA X/.16666666666667E0/
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      END
      SUBROUTINE SPLINT(PP,N,I,NP,A)
C
C     FIND COEFFICIENTS A(J),J=1,NP FOR SPLINE INTERPOLATION OF A
C          FUNCTION AT 'P' GIVEN F(I+J-1),J=1,NP
C            F(I+1) HAS THE VALUE AT P=I
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(4)
      DATA X/.16666666666667D0/
      P = PP
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      IF(P.EQ.0.) GO TO 10
      P1=P-1.
      P2=P-2.
      Q=P+1.
      C1 = P * P1 * X
      C2 = Q * P2 * .5
C     FF=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      NP = 4
      A(1) = -P2 * C1
      A(4) =   Q * C1
      A(2) =  P1 * C2
      A(3) =  -P * C2
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 NP = 1
      I = 1
      A(1) = 1
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 NP = 1
      I = N
      A(1) = 1.
      RETURN
   10 NP = 1
      I = I+1
      A(1) = 1.0
      RETURN
      END
      SUBROUTINE DISPLY(A,M,N,MA,DMAX)
	use io
      REAL*8   A(MA,N),EMAX
      REAL*8 DMAX
      CHARACTER CHARS(11),CHARSM(11),LINE(132),PLUS,CAPI,MINUS,SIGNS
      DATA  CHARS      / ' ','1',' ','3',' ','5',' ','7',' ','9','T' /
      DATA  CHARSM     / ' ','a',' ','c',' ','e',' ','g',' ','i','J' /
      DATA PLUS,MINUS,CAPI / '+','-','I' /
      MP= MIN(M,130)
      EMAX = 0
      DO 10 I=1,N
      DO 10 J=1,MP
10    EMAX = MAX(EMAX, ABS(A(J,I)) )
      DMAX = EMAX
      IF(EMAX.EQ.0) GOTO 100
      SCALE = 10/EMAX
C20    FORMAT(// 20('*') //)
      IF(ABS(LOG10(EMAX)).LE.4) THEN
      WRITE(KO,12) EMAX
       ELSE
      WRITE(KO,13) EMAX
       ENDIF
12    FORMAT(' Full scale =',F12.5/ )
13    FORMAT(' Full scale =',1P,E12.4/ )
      DO 50 I=1,N
      DO 30 J=1,MP
         ANO= SCALE * ABS(A(J,I))
         NO = ANO
         IS = SIGN(1D0,A(J,I))
          SIGNS = PLUS
           LINE(J) = CHARS(NO+1)
          if(IS.lt.0) then
          	SIGNS = MINUS
                LINE(J) = CHARSM(NO+1)
          endif
         IF(NO.LT.1 .AND. ANO .GE. 0.10) LINE(J) = SIGNS
30      CONTINUE
      WRITE(KO,35) (LINE(K),K=1,MP),CAPI
35    FORMAT(' I', 131A1)
50    CONTINUE
C     WRITE(KO,20)
      RETURN
100   WRITE(KO,105)
105   FORMAT(/' Everywhere zero???'/)
      RETURN
      END
      SUBROUTINE PLM(X,N,M,NA,PL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PL(NA,M+1),L,X
      N1 = N+1
      M1 = M+1
      DO 10 J=1,M1
      DO 10 I=1,N1
10    PL(I,J)=0.
      PL(1,1) = 1.
      if(NA>1) PL(2,1) = X
      SX = SQRT(1.-X*X)
      if(M>0) PL(2,2) = SX
	  FACT=1.
	  PMM=1.
	  DO 15 J=2,min(M1,N1)
	    mm = J-1
	    PMM = PMM*FACT*SX
		FACT=FACT+2.
		PL(J,J) = PMM
		if(J+1.le.N1) PL(J+1,J) = X*(2*mm+1.) * PL(J,J)
15      CONTINUE
	  
	  DO 20 J=1,M1
	   mm = J-1
	  DO 20 I=J+2,N1
	   ll = I-1
      PL(I,J)=((2.*ll-1.)*X*PL(I-1,J) - (ll+mm-1.)*PL(I-2,J))/(ll-mm)	  
20    CONTINUE
      RETURN
      END
      FUNCTION YLMC(L,M)
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(IN):: L,M
      PHASE(I) = (-1)**I
C
C      CALCULATE THE COEFFICIENT OF P(L,M)*E(I*M*PHI) IN Y(L,M)
C
         LF1 = L + 1
         MA = ABS(M)
                             R =  FACT(LF1-MA)-FACT(LF1+MA)
      R = SQRT((2*L+1)/(4*PI)*EXP(R))
     &      * PHASE(M)
      IF(M.LT.0) R = R * PHASE(MA)
      YLMC = R
      RETURN
      END
      SUBROUTINE LOGFAC(L)
	use factorials
C  FACT(I)= LN((I-1)!)
C  THE VALUE OF L CORRESPONDS TO THE DIMENSION OF (FACT(
      FACT(1)=0.
      DO 100 J=2,L
  100 FACT(J)=FACT(J-1)+LOG(J-1D0)
      RETURN
      END
      FUNCTION WIG3J(A,B,C,AL,BE,GA)
	use factorials
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8,intent(IN):: A,B,C,AL,BE,GA
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(AL+BE+GA) 11,10,11
11    WIG3J = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GOTO 11
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
28    TA = real(IX)*real(IY)*real(IZ)
      TB = MIN*real(IA+MIN)*real(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    IF(S.EQ.0.0) GO TO 11
      I = B-A+GA
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
      WIG3J = PHASE(I)  *  EXP(XDDD) * S
      RETURN
      END
C     FUNCTION WIG6J
C
C     THIS CALCULATES WIGNER 6-J COEFFICIENTS AS DEFINED IN BRINK AND
C     SATCHLER. IT REQUIRES TWO SUBPROGRAMS, FACT AND COMP... SEE
C     LISTING OF 3-J SYMBOL. IT USES THE CLOSED FORM FOR W-COEFICIENTS
C     DUE TO RACAH ALSO IN BRINK AND SATCHLER. THE SAME RULES APPLY FOR
C     THE INPUT PARAMETERS AS STATED IN THE 3-J PROGRAM.
C
      FUNCTION WIG6J(A,B,C,D,E,F)
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,B,C,D,E,F
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(FAIL3(A,B,C)) GO TO 10
      IF(FAIL3(A,E,F)) GO TO 10
      IF(FAIL3(B,D,F)) GO TO 10
      IF(FAIL3(C,D,E)) GO TO 10
      GO TO 14
   10 WIG6J=0.
      RETURN
   14 IC=A+B-C
      ID=E+D-C
      IE=A+E-F
      IG=B+D-F
      IA=C+F-A-D
      IB=C+F-B-E
      IH=A+B+E+D+1.
      M=MIN(IH,IC,ID,IE,IG)
      IF(M)10,17,17
   17 MUP=MIN(IA,IB,0)
      T=PHASE(M)
      N = M
      S=T
   18 M=M-1
      IF(M+MUP)24,16,16
   16 TA=real(IA+M+1)*real(IB+M+1)*real(IH-M)*real(M+1)
      TB=real(IC-M)*real(ID-M)*real(IE-M)*real(IG-M)
      T=-T*TA/TB
      S=S+T
      GO TO 18
   24 IT=A+B+C+1.
      IU=A+E+F+1.
      IV=B+D+F+1.
      IW=C+D+E+1.
      XD  =  .5*(FACT(1+IC)+FACT(1+IE+IB)+FACT(1+IA+IG)
     1+FACT(1+IE)+FACT(1+IB+IC)+FACT(1+IA+ID)+FACT(1+IG)+FACT(1+IC+IA)
     1+FACT(1+ID+IB)
     2+FACT(1+ID)+FACT(1+IA+IE)+FACT(1+IB+IG)-FACT(1+IT)-FACT(1+IU)
     3-FACT(1+IV)-FACT(1+IW))
     4+FACT(1+IH-N)-FACT(1+N)-FACT(1+IA+N)-FACT(1+IB+N)-FACT(1+IC-N)
     5-FACT(1+ID-N)-FACT(1+IE-N)-FACT(1+IG-N)
      WIG6J=PHASE(IH-1)  *S*EXP(XD)
      RETURN
 99    WIG6J = 1.0
      RETURN
      END
C     FUNCTION WIG9J(A,B,C,D,E,F,G,H,Z)
C
C     THIS CALCULATES 9-J SYMBOLS, OR X-COEFICIENTS AS DEFINED IN BRINK
C     AND SATCHLER. IT USES THE FORMULA FOR 9-J SYMBOLS IN TERMS OF 6-J
C     SYMBOLS. IT THEREFORE NEEDS THE 6-J SUBPROGRAM, AND THE CONDITION
C     ON THE INPUT PARAMETERS IS THE SAME AS FOR THE 3-J PROGRAM#      0
C
       FUNCTION WIG9J(A,B,C,D,E,F,G,H,Z)
       IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,B,C,D,E,F,G,H,Z
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(FAIL3(A,B,C)) GO TO 20
      IF(FAIL3(D,E,F)) GO TO 20
      IF(FAIL3(C,F,Z)) GO TO 20
      IF(FAIL3(A,D,G)) GO TO 20
      IF(FAIL3(B,E,H)) GO TO 20
      IF(FAIL3(G,H,Z)) GO TO 20
      GO TO 26
   20 WIG9J=0.
      RETURN
   26 S=0.
      XA=ABS(A-Z)
      XB=ABS(D-H)
      XC=ABS(B-F)
      X=XA
      IF(X-XB)10,11,11
   10 X=XB
   11 IF(X-XC)12,13,13
   12 X=XC
   13 IF(X-A-Z)14,14,15
   14 IF(X-D-H)16,16,15
   16 IF(X-B-F)17,17,15
   17 S=S+(2.*X+1.)*WIG6J(A,Z,X,H,D,G)*WIG6J(B,F,X,D,H,E)*WIG6J(A,Z,X,F,
     1B,C)
      X=X+1.
      GO TO 13
   15 IF(S-0.)18,20,18
   18 K=2.0*(A+B+D+F+H+Z)
      WIG9J=PHASE(K)*S
      RETURN
      END
      FUNCTION CLEB6(A,AL,B,BE,C,M)
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,AL,B,BE,C,M
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
28    TA = real(IX)*real(IY)*real(IZ)
      TB = MIN*real(IA+MIN)*real(IB+MIN)
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
      FUNCTION RACAH(A,B,C,D,E,F)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8,intent(IN):: A,B,C,D,E,F
      PHASE(I) = (-1)**I
      Z = ABS(A+B+C+D)
      I = Z + 0.5
      RACAH = PHASE(I) * WIG6J(A,B,E,D,C,F)
      RETURN
      END
      SUBROUTINE GAUSS(N,NR,A,SING,DET,EPS,SHOW)
C
C    SOLVE BY GAUSSIAN ELIMINATION SUM(J): A(I,J).P(J) = A(I,N+1)
C             WHERE P(J) IS LEFT IN MAT(J,N+1)
C
	use io
       IMPLICIT REAL*8(A-H,O-Z)
	parameter(M=1)
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
         WRITE(KO,3) K,DET/LOG(10.),A(K,K)
    3  FORMAT(//' THE MATRIX IS SINGULAR AT',I3,', Log10 determinant is
     &  ',2E16.8,' WITH ',2e16.8/)
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
      FUNCTION FFC(PP,F,N)
      COMPLEX*16 FFC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC=F(N)
      RETURN
      END
      FUNCTION FFR(PP,F,N)
      COMPLEX*16 F(N)
      REAL*8 PP
      REAL*8 FFR
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR = ( - P2*DBLE(F(I)) + Q*DBLE(F(I + 3)))*(P*P1*X)
     &        + (P1*DBLE(F(I+1)) - P*DBLE(F(I+2)))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR=F(N)
      RETURN
      END
      SUBROUTINE EXPAND(FNL,NLD,NLN,NLO,F,J,MLT)
      REAL*8 FNL(NLD,NLO),F(NLN)
      P = (J-1)/FLOAT(MLT)
      JJ = P
      IF(JJ.EQ.0) JJ = 1
      IF(JJ.GT.NLO-3) JJ = NLO - 3
      P = P - JJ
      IF(ABS(P).LT.1E-5) THEN
         DO 10 I=1,NLN
10       F(I) = FNL(I,JJ+1)
      ELSE
         P1 = P - 1.
         P2 = P - 2.
         Q  = P + 1.
         X = P * P1 / 6.0
         Y = Q * P2 * 0.5
         DO 20 I=1,NLN
20       F(I) = (-P2 * FNL(I,JJ) + Q * FNL(I,JJ+3)) * X
     &             +(P1 * FNL(I,JJ+1)    - P * FNL(I,JJ+2)  ) * Y
      ENDIF
      RETURN
      END
      SUBROUTINE ECPAND(FNL,NLD,NLN,NLO,F,J,MLT)
      COMPLEX*16 FNL(NLD,NLO),F(NLN)
      P = (J-1)/FLOAT(MLT)
      JJ = P
      IF(JJ.EQ.0) JJ = 1
      IF(JJ.GT.NLO-3) JJ = NLO - 3
      P = P - JJ
      IF(ABS(P).LT.1E-5) THEN
         DO 10 I=1,NLN
10       F(I) = FNL(I,JJ+1)
      ELSE
         P1 = P - 1.
         P2 = P - 2.
         Q  = P + 1.
         X = P * P1 / 6.0
         Y = Q * P2 * 0.5
         DO 20 I=1,NLN
20       F(I) = (-P2 * FNL(I,JJ) + Q * FNL(I,JJ+3)) * X
     &             +(P1 * FNL(I,JJ+1)    - P * FNL(I,JJ+2)  ) * Y
      ENDIF
      RETURN
      END
      SUBROUTINE QLM(X,N,M,NA,QL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 QL(NA,9),L
      N1 = N+1
      M1 = M+1
      DO 10 J=1,M1
      DO 10 I=1,N1
10    QL(I,J)=0.
      IF(ABS(X).EQ.1.) RETURN
      FL = LOG((1.+X)/(1.-X)) * 0.5
      QL(1,1) = FL
      QL(2,1) = X * FL - 1.0
      SX = SQRT(1.-X*X)
      QL(2,2) = -(X * QL(2,1) - QL(1,1)) / SX
      DO 20 I=3,N1
      L = I-1
      QL(I,1)=((2.*L-1.)*X*QL(I-1,1) - (L-1.)*QL(I-2,1))/L
      JM = MIN0(I,M1)
      DO 20 J=2,JM
      M = J-1
      MM= M-1
      QL(I,J)=-((L-MM)*X*QL(I,J-1) - (L+MM)*QL(I-1,J-1))/SX
20    CONTINUE
      RETURN
      END
      COMPLEX*16 FUNCTION FFCI(PP,F,N,L)
      INTEGER, INTENT(IN):: N,L
      COMPLEX*16, INTENT(IN):: F(N)
      REAL*8, INTENT(IN):: PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFCI=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFCI=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFCI=F(N) * ((N-1.)/PP)**(L+1)
      RETURN
      END

C    SUBROUTINE GETMASS TO INTERPRET A NUCLEUS NAME eg 28SI or SI28
C    AND RETURN REAL*8 VALUES OF THE MASSES AND ERROR.
C    INPUT:NAME (READ AS CHARACTER*8)
C    OUTPUT:RM    MASS  (AMU)................DOUBLE PRECISION
C           AM    MASS DEFECT     (KEV)......DOUBLE PRECISION
C           EM    ERROR IN MASS   (KEV)......DOUBLE PRECISION
C           IA    MASS  
C           IZ    ATOMIC NUMBER                               
C           IFLAG  (0)  O.K.
C                  (1)  REQUIRED SYMBOL NOT IN FILE
C                  (2)  MASSES ARE INTERPOLATED VALUES
C                  (-1) NO MASS NUMBER FOUND
C		   (other) IOSTAT from opening MASFIL
C
C
      SUBROUTINE GETMASS(RM,AM,EM,IA,IZ,IFLAG,NAME,MASFIL)
      CHARACTER*10 NOS
      CHARACTER*26 LETTS
      CHARACTER*8 LINEIN,NAME
      CHARACTER*2 ATOM,NUCNAM(110)
      INTEGER MAXIZ,IAMIN(0:4),IAMAX(0:4),IREC(0:4),INF,IINF
      CHARACTER FORMAT(3)*4
      REAL*8 RM,EM,AM,RMASS,ERRM
      CHARACTER * (*) MASFIL
      REAL CKZ,DCZ,DELCKZ,DELPA,DELSN,DELWAH,JM,MEA,MJ,
     +     MN,PA,SJ,SN,TEA,WAH
      DATA FORMAT /'(I1)','(I2)','(I3)'/
      DATA NUCNAM/'NN','H ','HE',
     +'LI','BE','B ','C ','N ','O ','F ','NE',
     +'NA','MG','AL','SI','P ','S ','CL','AR',
     +'K ','CA','SC','TI','V ','CR','MN','FE','CO',
     +'NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     +'RB','SR','Y ','ZR','NB','MO','TC','RU','RH',
     +'PD','AG','CD','IN','SN','SB','TE','I ','XE',
     +'CS','BA','LA','CE','PR','ND','PM','SM',
     +'EU','GD','TB','DY','HO','ER','TM','YB',
     +'LU','HF','TA','W ','RE','OS','IR','PT',
     +'AU','HG','TL','PB','BI','PO','AT','RN',
     +'FR','RA','AC','TH','PA','U ','NP','PU',
     +'AM','CM','BK','CF','ES','FM','MD','NO',
     +'LR','RF','HA','NH','NS','UO','UE'/
      DATA IFIRST/0/
      DATA MAXIZ/109/ 
      DATA AMU/931501.6/  
      DATA NOS/'0123456789'/
      DATA LETTS/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA ITYPE1/11/
C		Remove leading blanks
      LINEIN = ADJUSTL(NAME)
C		Convert to upper case
      CALL UCASE(LINEIN)
C		Replace common names
      CALL NICKS(LINEIN)
      RM=0.0
      AM=0.0
      EM=0.0
      IA=0
      IZ=0
      IFLAG=0
	if(LINEIN(1:3)=='GAM') return
      IF (IFIRST.NE.0) GOTO 10
      OPEN(7,STATUS='OLD',ACCESS='DIRECT',RECL=64,FILE=MASFIL,
     X	 ERR=5,IOSTAT=IFLAG)
    5	if(IFLAG/=0) then
	  write(6,*) 'Could not open mass table : ',MASFIL
	  write(6,*) 'Get from ftp://ftp.ph.surrey.ac.uk/',
     X		'pub/fresco/m88lrb if necessary.'
	  stop
	  endif
!    5	if(IFLAG/=0) continue
      IFIRST=1
   10 IFLAG=0
C
C     *READING ALONG TO FIND NUCLEUS*
C
	N1 = SCAN(LINEIN,NOS)
	N2 = SCAN(LINEIN,NOS,.true.)
	N1 = max(N1,N2-2)
	N = N2-N1+1
	if(N1==0) N=0
!			N1 is position of first digit
!			N2 is position of last digit
	NV = VERIFY(LINEIN(N1:N2),NOS)

	I1 = SCAN(LINEIN,LETTS)
	I2 = SCAN(LINEIN,LETTS,.true.)
	I = I2-I1+1
	if(I1==0.or.I>2) I=0
	IV = VERIFY(LINEIN(I1:I2),LETTS)
!			I1 is position of first letter
!			I2 is position of last letter
!	write(6,*) 'N1,N2,I1,I2 =',N1,N2,I1,I2,'V:',NV,IV
      ATOM=LINEIN(I1:I2)
	if(N==0..or.I==0.or.NV/=0.or.IV/=0) then
	  IFLAG=-1
	  return
	  endif
C
C     *OBTAINING A AND Z FOR NUCLEUS*
C
      READ(LINEIN(N1:N2),FORMAT(N)) IA

      IMAX=MAXIZ+1
      DO 110 I1=1,IMAX
      IF (ATOM.EQ.NUCNAM(I1)) GOTO 120
  110 CONTINUE
      GOTO 270
  120 IZ=I1-1
C
C     *OBTAINING DATA FROM MASS FILE*
C
  140 IF (IZ.LT.0.OR.IA.LT.0) GOTO 270
      IF (IZ.GT.MAXIZ) GOTO 270
!	return
      KEYPTR=IZ/5                
      READ(7,REC=KEYPTR+2)(IREC(I),IAMIN(I),IAMAX(I),I=0,4)
      IPTR=IZ-(KEYPTR*5)
      IF (IA.LT.IAMIN(IPTR).OR.IA.GT.IAMAX(IPTR)) GOTO 270
      KEYPTR=IA+IREC(IPTR)-IAMIN(IPTR)
      READ(7,REC=KEYPTR)PA,DELPA,DCZ,MN,MEA,CKZ,DELCKZ,SN,DELSN,
     +     TEA,SJ,JM,MJ,WAH,DELWAH,INF
C
C     decode the INF word, which contains information on which
C     values are present for each nucleus..so dave love says
C
      IINF=INF/2**(ITYPE1*2-2)
      IINF=MOD(IINF,4)
      IF (IINF.NE.0) THEN
        RMASS=DBLE(WAH)*1000.0D0
        ERRM=DBLE(DELWAH)*1000.0D0
      ELSE
        RMASS=DBLE(MN)*1000.0D0
        IFLAG=2
      ENDIF
C
C     *CALCULATION   DATA PUT INTO ARRAYS*
C
      IF (ERRM.LT.0.0) IFLAG=2
      RM=RMASS/DBLE(AMU)+DBLE(IA)
      AM=RMASS
      EM=ERRM 
	write(6,250) LINEIN,RM,IZ
250	format('   From nuclide mass table: ',A8,' has mass',f10.5,
     X		' and charge',i3)
      RETURN
C
C     *ERROR MESSAGES*
C
  270 IFLAG=1
      RM=999.0 
      EM=999.0   
      IA=0 
      IZ=0
      RETURN
      END
      SUBROUTINE UCASE(LINEIN)
*     
*     CONVERT LOWER CASE TO UPPER CASE IF APPROPRIATE
*
      CHARACTER*8 LINEIN
      DO 20 J=1,8
      IF (LGE(LINEIN(J:J),'a').AND.LLT(LINEIN(J:J),'z')) THEN
       LINEIN(J:J)=CHAR(ICHAR(LINEIN(J:J))-ICHAR('a')+ICHAR('A'))
      ENDIF
 20   CONTINUE
      RETURN
      END
      SUBROUTINE NICKS(N)
*     
*     CONVERT COMMON NUCLEAR NAMES (at least 3 letters required)
*
      CHARACTER*8 N
	if(N(1:3)=='NEU') N = '1NN'
	if(N(1:3)=='PRO') N = '1H'
	if(N(1:3)=='DEU') N = '2H'
	if(N(1:3)=='TRI') N = '3H'
	if(N(1:3)=='HEL') N = '3HE'
	if(N(1:3)=='ALP') N = '4HE'
	if(N(1:3)=='ALF') N = '4HE'
	return
	end
      SUBROUTINE HDIAG(H,N,NSMAX,IEGEN,U,NR)
C     MIHDI3,FORTRAN 22 DIAGONALIZATION OF A REAL SYMMETRIC MATRIX BY
C        THE JACOBI METHOD.
C     MAY 19, 1959
C     CALLING SEQUENCE FOR DIAGONALIZATION
C            CALL  HDIAG(H, N, IEGEN, U, NR)
C     IEGEN MUST BE SET UNEQUAL TO ZERO IF ONLY EIGENVALUES ARE
C             TO BE COMPUTED.
C
C     N IS THE ORDER OF THE MATRIX.  H.
C            WHERE H IS THE ARRAY TO BE DIAGONALIZED.
C     IEGEN MUST BE SET EQUAL TO ZERO IF EIGENVALUES AND EIGENVECTORS
C           ARE TO BE COMPUTED.
C
C     U IS THE UNITARY MATRIX USED FOR FORMATION OF THE EIGENVECTORS.
C
C     NR IS THE NUMBER OF ROTATIONS.
C
C     A DIMENSION STATEMENT MUST BE INSERTED IN THE SUBROUTINE.
C     DIMENSION H(N,N), U(N,N), X(N), IQ(N)
C
C
C     THE SUBROUTINE OPERATES ONLY ON THE ELEMENTS OF H THAT ARE TO THE
C             RIGHT OF THE MAIN DIAGONAL. THUS, ONLY A TRIANGULAR
C             SECTION NEED BE STORED IN THE ARRAY H.
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(NSMAX,NSMAX), U(NSMAX,NSMAX), X(NSMAX), IQ(NSMAX)
      ONE = 1d0
      TWO = 2d0
      IF(IEGEN) 15,10,15
   10 DO 14 I = 1,N
      DO 14 J = 1,N
      IF(I - J) 12,11,12
   11 U(I,J) = 1.0
      GO TO 14
   12 U(I,J) = 0.
   14 CONTINUE
C
   15 NR=0
      IF(N - 1) 1000,1000,17
C
C     SCAN FOR LARGEST OFF DIAGONAL ELEMENT IN EACH ROW
C     X(I) CONTAINS LARGEST ELEMENT IN ITH ROW
C     IQ(I) HOLDS SECOND SUBSCRIPT DEFINING POSITION OF ELEMENT
C
   17 NMI1 = N - 1
      DO 30 I = 1,NMI1
      X(I) = 0.
      IPL1 = I + 1
      DO 30 J = IPL1,N
      IF(X(I) - ABS (H(J,I))) 20,20,30
   20 X(I) = ABS (H(J,I))
      IQ(I) = J
   30 CONTINUE
C
C     SET INDICATOR FOR SHUT-OFF.RAP=2**-27,NR = NO. OF ROTATIONS
      RAP=7.450580596d-9
      HDTEST=1.0d100
C
C     FIND MAXIMUM OF X(I) S FOR PIVOT ELEMENT AND
C     TEST FOR END OF PROBLEM
C
   40 DO 70 I = 1,NMI1
      IF(I - 1) 60,60,45
   45 IF(XMAX-X(I)) 60,70,70
   60 XMAX = X(I)
      JPIV = IQ(I)
      IPIV = I
   70 CONTINUE
C
C     IS MAX. X(I) EQUAL TO ZERO, IF LESS THAN HDTEST, REVISE HDTEST
      IF(XMAX) 1000,1000,80
   80 IF (HDTEST) 90,90,85
   85 IF(XMAX - HDTEST) 90,90,148
   90 HDIMIN = ABS (H(1,1))
      DO 110 I = 2,N
      IF(HDIMIN - ABS (H(I,I))) 110,110,100
  100 HDIMIN = ABS (H(I,I))
  110 CONTINUE
C
      HDTEST = HDIMIN*RAP
C
C     RETURN IF MAX.H(J,I)LESS THAN(2**-27)ABSF(H(K,K)-MIN)
      IF (HDTEST - XMAX) 148,1000,1000
  148 NR = NR + 1
C
C     COMPUTE TANGENT, SINE AND COSINE, H(I,I),H(J,J)
  150 TANG = SIGN (TWO,(H(IPIV,IPIV)-H(JPIV,JPIV)))*H(JPIV,IPIV)/(ABS (H
     1(IPIV,IPIV)-H(JPIV,JPIV))+SQRT ((H(IPIV,IPIV)-H(JPIV,JPIV))**2+4.0
     2*H(JPIV,IPIV)**2))
      COSINE = 1.0/SQRT (1.0+TANG**2)
      SINE = TANG*COSINE
      HII = H(IPIV,IPIV)
      H(IPIV,IPIV) = COSINE**2*(HII+TANG*(2.*H(JPIV,IPIV)+TANG*H(JPIV,JP
     1IV)))
      H(JPIV,JPIV)=COSINE**2*(H(JPIV,JPIV)-TANG*(2.*H(JPIV,IPIV)-TANG*H
     1II))
      H(JPIV,IPIV)= 0.
C
C     PSEUDO RANK THE EIGENVALUES
C     ADJUST SINE AND COS FOR COMPUTATION OF H(IK) AND U(IK)
      IF(H(IPIV,IPIV) -H(JPIV,JPIV)) 152,153,153
  152 HTEMP = H(IPIV,IPIV)
      H(IPIV,IPIV) = H(JPIV,JPIV)
      H(JPIV,JPIV) = HTEMP
      HTEMP = SIGN (ONE, -SINE)*COSINE
      COSINE = ABS (SINE)
      SINE = HTEMP
  153 CONTINUE
C
C     INSPECT THE IQS BETWEEN I + 1 AND N-1 TO DETERMINE
C     WHETHER A NEW MAXIMUM VALUE SHOULD BE COMPUTED SINCE
C     THE PRESENT MAXIMUM IS IN THE I OR J ROW.
C
      DO 350 I = 1,NMI1
      IF(I-IPIV) 210,350,200
  200 IF(I-JPIV) 210,350,210
  210 IF(IQ(I)-IPIV) 230,240,230
  230 IF(IQ(I)-JPIV)350,240,350
  240 K = IQ(I)
  250 HTEMP=H(K,I)
      H(K,I) = 0.
      IPL1 = I + 1
      X(I) = 0.
C
C     SEARCH IN DEPLETED ROW FOR NEW MAXIMUM
C
      DO 320 J = IPL1,N
      IF(X(I)-ABS (H(J,I))) 300,300,320
  300 X(I) = ABS (H(J,I))
      IQ(I) = J
  320 CONTINUE
      H(K,I) = HTEMP
  350 CONTINUE
C
      X(IPIV) = 0.
      X(JPIV) = 0.
C
C     CHANGE THE OTHER ELEMENTS OF H
C
      DO 530 I = 1,N
C
      IF(I-IPIV) 370,530,420
  370 HTEMP = H(IPIV,I)
      H(IPIV,I) = COSINE*HTEMP + SINE * H(JPIV,I)
      IF(X(I) - ABS (H(IPIV,I)))380,390,390
  380 X(I) = ABS (H(IPIV,I))
      IQ(I) = IPIV
  390 H(JPIV,I) = -SINE*HTEMP + COSINE*H(JPIV,I)
      IF(X(I) - ABS (H(JPIV,I))) 400,530,530
  400 X(I) = ABS (H(JPIV,I))
      IQ(I) = JPIV
      GO TO 530
C
  420 IF(I-JPIV) 430,530,480
  430 HTEMP = H(I,IPIV)
      H(I,IPIV) = COSINE*HTEMP + SINE*H(JPIV,I)
      IF(X(IPIV) - ABS (H(I,IPIV)) ) 440,450,450
  440 X(IPIV) = ABS (H(I,IPIV))
      IQ(IPIV) = I
  450 H(JPIV,I) = -SINE*HTEMP + COSINE*H(JPIV,I)
      IF(X(I) - ABS ( H(JPIV,I)) ) 400,530,530
C
  480 HTEMP = H(I,IPIV)
      H(I,IPIV) = COSINE*HTEMP + SINE*H(I,JPIV)
      IF(X(IPIV) - ABS ( H(I,IPIV)) ) 490,500,500
  490 X(IPIV) = ABS (H(I,IPIV))
      IQ(IPIV) = I
  500 H(I,JPIV) = -SINE*HTEMP + COSINE*H(I,JPIV)
      IF(X(JPIV) - ABS ( H(I,JPIV))) 510,530,530
  510 X(JPIV) = ABS (H(I,JPIV))
      IQ(JPIV) = I
  530 CONTINUE
C
C     TEST FOR COMPUTATION OF EIGENVECTORS
C
      IF(IEGEN) 40,540,40
  540 DO 550 I = 1,N
      HTEMP = U(IPIV,I)
      U(IPIV,I) = COSINE*HTEMP + SINE*U(JPIV,I)
  550 U(JPIV,I) = -SINE*HTEMP + COSINE*U(JPIV,I)
      GO TO 40
 1000 RETURN
      END

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


*****POTENT*************************************************************
      SUBROUTINE POTENT(FORMF,NF,PTYPE,N,HCM,TRENEG,MR,MATCH,HASO,
     X              HP,HPOT,STREN,LAMBDA,NCHPMAX,VARYL,LDEP,LSHAPE,
     X              MATRIX,MEK,NIX, JEX,MASS,CPOT,NEX,NCHAN,COPY,STOP19,
     x              FORMDEF)
	use parameters
	use io
	use searchpar
	use fresco1, only: pluto,npluto
	!  N here is actually MAXM
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(INFORM=4,MMX=2)
      COMPLEX*16 FORMF(MAXM,MLOC),V,CO(MMX*MAXM),FFC4
      REAL*8 HCM, JEX(6,MXP,MXX),MASS(4,MXP),Z8,JA,JB,KG
     X       ,H,HP(MXP),HPOT(MLOC),STREN(MLOC),VARYL(2,MLOC),CLEB6
      REAL*8 RL(MMX*N),P(0:8),DEF(0:8),CCPT(2),MEK(MPAIR),AI(4),
     x       IM(MMX*N),FORMDEF(1:7,MLOC)
      LOGICAL PT(2),WRIT25,itt,LDEPP,HASO(KPMAX),STOP19,nosub,DUPOT,SHNG
      INTEGER KPOP(KPMAX)
      INTEGER PTYPE(8,MLOC),SHAPE,SHAPEI,TRENEG,TYPE,TYPEI,LOC(0:7),
     X        MATRIX(6,MPAIR), NEX(MXP),CPOT(MXP,MXX),COPY(2,MXP,MXX),
     X        LAMBDA(MLOC),NCHPMAX(MXP),KKA,LSHAPE(MLOC),LDEP(MLOC)
      EQUIVALENCE (P1,DEF(1)),(P2,DEF(2)),(P3,DEF(3)),(P4,DEF(4)),
     &            (P5,DEF(5)),(P6,DEF(6)),(P0,DEF(7)),(P,DEF)
      CHARACTER*11 WHO(31),SHP(0:14),METYPE(2),SHPNAM
      CHARACTER*5 NNKIND(12),PTNAME(2)
      CHARACTER cit,JL
      CHARACTER*7 datavals
      CHARACTER*80 COMMENT
      DATA NNKIND /'V1S0','V3S1','V3S3D','V3D1','V1P1','V3P0','V3P1',
     &             'V3P2','V3P3F','V3F2','V1D2','V3D2' /
      DATA WHO /'Coulomb','Volume','Surface','Projtl S.O.','Target S.O.'
     X         ,'Projtl Tr','Target Tr','Prj+Targ Tr',
     X       	'Spin.spin','UNDEFINED',
     X          'Proj. defrm','Target defm','Proj. table','Targ. table',
     X          'Proj.2order','Targ.2order',
     X       	'Mut.2order','All-orders',2*'UNDEFINED',
     X          'NN: SSC(C)','NN: User"s',8*'UNDEFINED','L(L+1) Vol'/
      DATA SHP /'Fourier-Bsl','Woods-Saxon','WS squared','Gaussian',
     X      'Yukawa','Exponential','RSC T=0','RSC T=1',
     X      'Read:real','Read:imag','Read:complx',
     X      'derivatives','quadrat,k>0','quadrat&k=0',
     X      'quad-pl&k=0' /,
     X      PTNAME / 'proj.','targ.' /, PI / 3.1415926 /,
     X      METYPE / 'Coulomb', 'Nuclear' /
	integer line(mvars)
	namelist/potl/ kpi,typei,itt,nosub,shapei,p,LSHAPEI,XLVARY,
     x        ALVARY,JL
	namelist/step/ib,ia,k,str
C
C  READ(I3,I2,A1,I2,7F8.4) KP,TYPE,IT,SHAPE,P1,P2,P3,P4,P5,P6,P7
C
C  TYPE  IT,   P1,   P2,   P3,   P4,   P5,   P6,P0
C       SHAPE
C    0   -     A#1   A#2   r0C   ac                Size & Coulomb
C
C              ------real------  ---imaginary----  -spin tensor-
C    1  SHAPE  Vr    Vrr0  Vra   Vi    Vir0  Via   Central Volume
C    2  SHAPE  Wr    Wrr0  Wra   Wi    Wir0  Wia   Central Derivative
C    3  SHAPE  Vso   Vsor0 Vsoa  Vsoi  Vsir0 Vsia  Projectile Spin-orbit
C    4  SHAPE  Uso   Usor0 Usoa  Usoi  Usir0 Usia  Target spin-orbit
C    5  SHAPE  VTr   VTrr0 VTra  VTri  VTir0 VTia  Projectile Tr tensor
C    6  SHAPE  UTr   UTrr0 UTra  UTri  Utir0 Utia  Target Tr tensor
C    7  SHAPE  T2r   T2rr0 T2ra  T2ri  T2ir0 T2ia  Joint spin-spin tensor
C    8  SHAPE  S.S   SSrr0 SSra  SSri  SSir0 SSia  Joint spin-spin scalar
C
C   10  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Projectile Deform
C   11  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Target Deformed
C   12  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Projectile table
C   13  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Target table
C   14  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Projectile 2-order table
C   15  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Target 2-order table
C   16  SHAPE  Def1  Def2  Def3  Def4  Def5  Def6,Def0 Mutual second-order
C   17                                              Proj-target all-order matrix
C               (Inelastic monopole = Def0) 
C
C   20  MAX#    Super-soft core of de Tourreil & Sprung (SSC (C))
C   21  MAX#    User-supplied N-N potential via subroutine NNPOT
C   30  SHAPE  VLL   VLLr0 VLLa  WLL   WLLr0 WLLa  Central L(L+1) 
C
C If -8 < TYPE < 0, then merge with the previous form,
C                   to make a summed form of type abs(TYPE).
C
C  KP=0 ends reading in potentials,
C  KP<0 indicates this is the last data card, use abs(KP) for this one.
C
C  all the r0 radii are multiplied by CC,
C   which is updated by a TYPE=0 card  to CC =  A#1**1/3 + A#2**1/3
C         & by other cards if P0.gt.0  to CC =  P0**1/3
C
C    IT = blank  : only off-diagonal potentialsiterated.
C       ne blank : potential component iterated, diagonal or off-diag.
C
C    SHAPE = 0 :  Woods-Saxon
C            1 :  Woods-Saxon squared
C            2 :  Gaussian
C            3 :  Yukawa
C            4 :  exponential
C            5 :  Reid soft core for T=0
C            6 :  Reid soft core for T=1
C            7 :  read:real             Scale real part by P1,
C            8 :  read:imaginary        and imag part by P2,
C            9 :  read:complex          and use P3 as r0 for radius defn
C      -7,-8,-9:  as for 7,8 & 9,  but rewind input file first
C           10 :  deformation by derivative of previous form factor.
C           11 :  deformation & projection of previous form factor, k>0
C           12 :  deformation & projection of previous form f., all k.
C           13 :  deformation & projection of previous form f., all k.,
C       		with no volume conservation term.
C           10-19: write SHAPE-10 to file 25 (TYPE 1 to 7 only)
C           20-24:  no shape: read in L-dependent shape from file=SHAPE,
C                           RECORD = 1 for + and 2 for - parity
C           25-29:  no shape: read in L-dependent shape from file=SHAPE,
C                           RECORD = [JT]+1 FOR EACH TOTAL CC SPIN JT
C           30-39:  use SHAPE-30 * L-dependent form factors from Cards 10.5,
C           40:  use potential P1 for + parity, P2 for - parity 
C           41:  use potential P(L+1) for partial waves L<=5,  P0 for all L>5
C           42:  use potential P(L+1) for CC sets J<=5,  P0 for all J>5 (L=int(J))
C
C           -1:  Fourier-Bessel shape
C        <= -2:  <new shapes go here>

C      for TYPE = 1-8 or 30 (incl), must have  -1 <= SHAPE <= 9
C          TYPE = 0,16              SHAPE not used
C                      (both unless SHAPE is 20 or above)
C          TYPE = 10 - 13 incl,  must have  7 .le. SHAPE .le. 13
C          TYPE = 14 - 17 incl,  must have  7 .le. SHAPE .le. 10
C          TYPE = 20 or 21       must have  1 .le. SHAPE .le. 12
C
C   for JF=1,NF the array PTYPE(i,JF) has information about potentials:
C
C       PTYPE(1,JF)  =  KP  :  number for referring to this potential
C             2      =  TYPE:    (as above, from 0 to 30)
C             3      =  k   :  multipole deformation (7=inelastic monopole)
C                                (but if TYPE=20,21, k=particular NN potl)
C             3      =  ICK :  for TYPE=17 (partition number) on first line
C             4      =  L   :  JF number of undeformed potential shape
C                                  (if -1, ignore this JF)
C             5      =  ITYP:  TYPE of original JF form, for TYPE ge 10
C                       SHAPE: file to read L-dep shapes for TYPE lt 10
C             6      =  0,1 :  1 to force potential to be iterated.
C             7,8    unused
C
C   for IX=1,NIX the array MATRIX(i,IX) has info about matrix elements:
C
C      MATRIX(1,IX)  =  IB  :  state fed by coupling
C             2      =  IA  :  state feeding the state IB
C             3      =  k   :  multipole deformation
C             4      =  NF  :  JF number of deformation multipole
C             5      =  NFD :  JF number of original potential
C             6      =  NUC :  1 for projectile, 2 for target
C         MEK(IX)    =  Coefficient of coupling NF to give M(Ek) desired
C
      TH = 1.0/3.0
C     conls is the constant factor for spin-orbit forces
      CONLS = 2.000
      DEF(0) = 1.0
      Z = 0.0
      RR = 0.0
      if(final) WRITE(KO,5)
5     FORMAT(//' ',132('*'),
     &  //' The following POTENTIALS are defined :'//
     & '  KP#       TYPE      IT    SHAPE      at  ',
     & '   V1        r1        a1     ',
     & '   V2        r2        a2        A-in     A-used'/)
      NF1 = NF+1
      DO 6 M=NF1,MLOC
      LAMBDA(M)=0
      LDEP(M) = 0
      DO 6 I=1,N
      FORMF(I,M) = (0d0,0d0)
6     CONTINUE
      KPOP(:) = -1
	do i=1,npluto
	 if(pluto(i)>0) KPOP(pluto(i))=0
	enddo
      NFD = 0
C     NIX = 1
      KPLAST = 0
      IR0 = 0
      CC = 0.0
	line(:) = 0
10    continue
!	READ(KI,12) KPI,TYPEI,itt,SHAPEI,(DEF(K),K=1,7)
	kpi=0; typei=0; itt=.false.;shapei=0; def(:)=0; nosub=.false.
	READ(KI,nml=potl)
	if(.not.nosub) then
	 cit = ' '
	 if(itt) cit = '1'
	else
	 cit = '2'
	 if(itt) cit = '3'
	endif
	LDEPP = SHAPEI>=30 .and. SHAPEI<40
	if(LDEPP)  SHAPEI = SHAPEI-30
	if(LDEPP.and..not.(JL.eq.'J' .or. JL.eq.'L')) JL='J'
	LDEPP=LDEPP .or. JL.eq.'J' .or. JL.eq.'L'
!12    FORMAT(I3,I2,L1,I2,7F8.4)
      IF(KPI.EQ.0) GO TO 200
      KP = ABS(KPI)
	if(KPOP(KP)==0) then
	 write(900+KP,880) HCM,(N-1)*HCM
880	 format('&pluto',/,'  frescocall=T   h=',f7.3,' rmax=',f7.2)
	 KPOP(KP) = 1
	endif
!		Adjust any search parameter!
	do ip=1,nvars
	if(srch_kind(ip)==1.and.abs(KP)==srch_kp(ip)) then
	  line(ip) = line(ip) + 1
	  if(line(ip)==srch_pline(ip)) then
	    if(abs(srch_value(ip)-nul)>.001) then
              DEF(srch_col(ip)) = srch_value(ip)
	    else
	     srch_value(ip) =  DEF(srch_col(ip))
	    endif
          endif
	endif
	enddo
      IF(KP.NE.KPLAST) THEN
C                      find step size for this potential:
       ICK = 0
       DUPOT = .false.
	NEXK = 0
       DO  1202 IC=1,NCHAN
       DO  1202 IA=1,NEX(IC)
         IF(CPOT(IC,IA).NE.KP) GO TO 1202
            IF(INH.GT.0 .AND. ICK.NE.0 .AND. ICK.NE.IC) THEN
             WRITE(KO,*) '***** ERROR, AS INH =',INH, ', SO SEPARATE ',
     X            'POTENTIALS NEEDED FOR PARTITIONS ',ICK,' &',IC
C                STOP 'POTENTIALS'
                 CALL ABEND(4)
            ENDIF
	  DUPOT = ICK>0.and.ICK.ne.IC	! if potential used in another partition
         ICK = IC
 1202  CONTINUE
       H = HCM
       IF(ICK.NE.0) H = HP(ICK)
       H4 = H
	if(ICK>0) NEXK = NEX(ICK)
      ENDIF
         IF(KPLAST.NE.0 .AND. KP.NE.KPLAST.and.final) WRITE(KO,174)
         KPLAST = KP
      TYPE = ABS(TYPEI)
 	SHNG = SHAPEI<0
	if(7<=abs(SHAPEI).and.abs(SHAPEI)<=9) SHAPEI = abs(SHAPEI)   ! just use SHNG then, for rewind option!
      SHAPE = SHAPEI   ! same meaning otherwise & now on
!	HASO(KP) = HASO(KP).or.TYPEI==3
!	write(KO,*) 'HASO',KP , HASO(KP),TYPEI==3
         MF = NF+1
      IF(TYPE.GE.10.and.TYPE.le.16) THEN
               ITYP = PTYPE(2,NFD)
               IF(SHAPE.EQ.0.OR.SHAPE.GT.13) SHAPE = 10
	       if(RR<1e-5) RR = 1.
               IF(TYPE>=14.and.TYPE<=16) SHAPE = min(10,SHAPE)
C                        remove next line when Coul quadrature written.
               IF(ITYP.EQ.0) SHAPE = MIN(SHAPE,10)
            DO 13 K=0,7
               LOC(K) = 0
            IF(K.NE.0 .AND. ABS(DEF(K)).LT.1E-9) GO TO 13
            IF(K.EQ.0 .AND. SHAPE.NE.12.and.SHAPE.NE.13) GO TO 13
               NF= NF + 1
               CALL CHECK(NF,MLOC,24)
            HPOT(NF) = H
            PTYPE(1,NF) = PTYPE(1,NFD)
            PTYPE(2,NF) = TYPE
            PTYPE(3,NF) = K
            PTYPE(4,NF) = NFD
            PTYPE(5,NF) = ITYP
            PTYPE(6,NF) = 0
            IF(itt) PTYPE(6,NF) = 1
            IF(nosub) PTYPE(6,NF) = PTYPE(6,NF)+2
               IF(K.EQ.0) PTYPE(4,NFD) = -NF
            PTYPE(7,NF) = 0
            PTYPE(8,NF) = 0
            LOC(K) = NF
	      LDEP(NF) = 0
	      if(JL.eq.'J') LDEP(NF) = 1
	      if(JL.eq.'L') LDEP(NF) = 2
	      LSHAPE(NF) = LSHAPEI
	      VARYL(1,NF) = XLVARY
	      VARYL(2,NF) = ALVARY
13        CONTINUE
       FORMDEF(1:7,NF) = DEF(1:7)
       IF(ITYP.GT.0) THEN
         WRITE(KO,175) KP,TYPE,WHO(TYPE+1),cit,SHAPE,SHP(SHAPE+1),MF,
     X             (DEF(K),K=1,7),RR
       ELSE
         WRITE(KO,176) KP,TYPE,WHO(TYPE+1),cit,SHAPE,SHP(SHAPE+1),MF,
     X             (DEF(K),K=1,7),RR
       ENDIF
C
       IF(TYPE.LT.12.or.TYPE.ge.18) GO TO 1350
1301        continue
!         READ(KI,*) J,(MATRIX(I,NIX),I=2,3),MEK(NIX)
	 ib=0; ia=0; k=0; str=0
         READ(KI,nml=step)
	   J = ib
         	IF(J.EQ.0) GO TO 1350
           	CALL CHECK(NIX,MPAIR,16)
!	Adjust any search parameter!
		do ip=1,nvars
		if(srch_kind(ip)==1.and.abs(KP)==srch_kp(ip)) then
		  line(ip) = line(ip) + 1
		  if(line(ip)==srch_pline(ip)) then
		    if(abs(srch_value(ip)-nul)>.001) then
       		       str = srch_value(ip)
		    else
		      srch_value(ip) =  str
		    endif
       		   endif
		endif
		enddo
            MATRIX(6,NIX) = 1
            if(k<0.or.ia<0) MATRIX(6,NIX) = 2
	   MATRIX(2,NIX) = abs(ia)
            KK = k
	    if(KK==0) k=7
	   MATRIX(3,NIX) = abs(k)
	   MEK(NIX) = str
	  
!1302     FORMAT(4X,3I4,F8.4)
            MATRIX(1,NIX) = ABS(J)
            K = MATRIX(3,NIX)
            M = LOC(K)
            MATRIX(4,NIX) = M
            MATRIX(5,NIX) = NFD
         WRITE(KO,1305) METYPE(MIN(ITYP+1,2)),(MATRIX(I,NIX),I=1,2),
     X    KK,(MATRIX(I,NIX),I=4,5),PTNAME(MATRIX(6,NIX))(1:1),MEK(NIX)
             IF(M.EQ.0) THEN
               WRITE(KO,*) 'MULTIPOLE ',KK,' NOT DEFINED!!'
                GO TO 1340
               ENDIF
1305        FORMAT(10X,A7,' reduced matrix element ',I3,' <-',I3,
     X      ' of k =',I2,' (at',I4,' from',i4,1x,a1,')  ==',F10.4)
            MEK(NIX) = MEK(NIX) / DEF(K)
1340        NIX = NIX + 1
            IF(J.GT.0) GO TO 1301
1350        CONTINUE
C
            I = MOD(TYPE,2) + 1
            R = RR * CCPT(I) / CC
         DO 14 M=MF,NF
            K = PTYPE(3,M)
	    KK=K
	    if(K==7) KK=0
!		write(6,*) ' @@ look at M =',M,' with K,KK,T= ',K,KK,TYPE
         IF(K.EQ.0) THEN
            WRITE(KO,171) M,KK
        ELSE IF(ITYP.GT.0 .AND. TYPE.LE.11) THEN
C                  plain nuclear deformations
            WRITE(KO,171) M,KK,DEF(K)/RR
            IF(PT(I) .AND. DEF(K).NE.Z)
     &         WRITE(KO,1715) PTNAME(I),DEF(K)/R,R
        ELSE IF(ITYP.eq.0 .AND. TYPE.LE.11.and..TRUE.) THEN
C                Coul deformations using Mn(Ek)
!      DO 1400 IC=1,NCHAN
	 IC = ICK
	 IF(IC.eq.0) GO TO 1400
         IF(CPOT(IC,1).NE.KP) GO TO 1400
         Z1Z2 = MASS(2+I,IC)
           BETA = DEF(K) * 4 * PI/(3 * RR **KK * Z1Z2)
           WRITE(KO,171) M,KK,BETA
            IF(PT(I) .AND. DEF(K).NE.Z) THEN
              BETA = DEF(K)  * 4*PI / (3 * R**KK * Z1Z2)
               WRITE(KO,1715) PTNAME(I),BETA,R
            ENDIF
1400     CONTINUE
        ELSE
C               need to know state spins, to determine beta:
            Z8 = Z
!      DO 149 IC=1,NCHAN
	if(ICK==0) goto 1491
	IC = ICK
      DO 149 IA=1,NEX(IC)
         IF(CPOT(IC,IA).NE.KP) GO TO 149
      DO 147 IB=1,NEX(IC)
	DO 146 IS=1,2
           R1 = 1.0
        IF(TYPE.GE.12.and.TYPE.le.16) THEN
C                                find table entry and I
           R1 = Z
           DO 145 J=1,NIX-1
            IF(MATRIX(1,J).NE.IB .OR. MATRIX(2,J).NE.IA .OR.
     X         MATRIX(3,J).NE.K  .OR. MATRIX(4,J).NE.M ) GO TO 145
           R1 = R1 + MEK(J)
	     if(TYPE>=14.and.TYPE<=16) then
		I=MATRIX(6,J)
		R1 = MEK(J)
!		write(6,*) ' @@ Set I =',I,R1
	  	if(I==IS) GO TO 1451
**************  If two J solutions (e.g, proj,targ), do both separately!
		go to 145
	     endif
145        CONTINUE
1451       CONTINUE
         ENDIF
	  if(I/=IS) GO TO 146
         IF(COPY(I,IC,IA).NE.0.or.COPY(I,IC,IB).NE.0) GO TO 146
      Z1Z2 = 1.0
      IF(ITYP.EQ.0) Z1Z2 = MASS(2+I,IC)
      JA=JEX(I,IC,IA)		
      KG=JEX(I+2,IC,IA)
      JB=JEX(I,IC,IB)
      PH = (-1)**NINT( (JA-JB + ABS(JA-JB))*0.5 )
      RC = Z1Z2 * PH * SQRT(2.*JA+1.) * CLEB6(JA,KG,KK+Z8,Z8,JB,KG)
!		write(6,*) ' @@ JA,JB,KG,KK,RC,I=',JA,JB,KG,KK,RC,I
      IF(ABS(RC).LT.1E-5) GO TO 146
      IF(TYPE.le.11) RC = Z1Z2 
!		write(6,*) ' @@ for IB,IA =',IB,IA,' MEK=',R1
!		write(6,*) ' @@ DEF,RR,RC,ITYP =',DEF(K),RR,RC,ITYP
       BETA               = DEF(K) * R1 / (RR * RC)
       IF(ITYP.EQ.0) BETA = DEF(K) * R1 * 4 * PI/(3 * RR ** KK  *  RC)
         WRITE(KO,1705) IB,IA,KK,BETA,IC,M
1705     FORMAT(14X,'Step',I2,' <-',I2,' of k =',I3,:,
     &          ' (Potl. beta =',F8.5,')',I3,' @ ',I3)
        R = RR * CCPT(I) / CC
          IF(PT(I).AND.BETA.NE.0.0) THEN
                           BETA = DEF(K) * R1 / (R * RC)
             IF(ITYP.EQ.0) BETA = DEF(K) * R1  * 4*PI / (3 * R**KK * RC)
             WRITE(KO,1715) PTNAME(I),BETA,R
            ENDIF
146       CONTINUE
147       CONTINUE
149       CONTINUE
1491      CONTINUE
        ENDIF
14      CONTINUE
        GO TO 195
      ELSE IF(TYPE.eq.17) then
C					All-order matrix couplings
        ITYP = PTYPE(2,NFD)
        WRITE(KO,177) KP,TYPE,WHO(TYPE+1),cit
	if(ICK==0) then
	   write(KO,*) ' Cannot find partition using TYPE=17!'
	   stop
	  endif
	I = NCHPMAX(ICK)*(NCHPMAX(ICK)+1)/2
	write(KO,640) KP,ICK,NCHPMAX(ICK),I
640	format('     Potential ',I3,' used in partition ',i3,
     x	  ', so ',i4,' channels and ',i5,' coupling forms required')
	if(DUPOT) then
	   write(6,*) 'ALL-ORDER COUPLINGS CANNOT BE USED IN MORE THAN',
     x		' ONE PARTITION! Stop.'
	   stop
	endif
	NF = NF + I
        CALL CHECK(NF,MLOC,24)
	STOP19 = .true.
C					Read in couplings as deformation lengths
650	 ib=0; ia=0; k=0; str=0
         READ(KI,nml=step)
	   J = ib; ib=abs(ib)
	   if(J==0) go to 660
           CALL CHECK(NIX,MPAIR,16)
	   MATRIX(1,NIX) = ib
	   MATRIX(2,NIX) = abs(ia)
		KK=abs(k)
		if(KK==0) k=7
	   MATRIX(3,NIX) = KK
	   MATRIX(4,NIX) = 0
	   MATRIX(5,NIX) = NFD
	   MATRIX(6,NIX) = 1
	   if(ia<0.or.k<0) MATRIX(6,NIX) = 2
	   MEK(NIX) = str
           WRITE(KO,655) METYPE(MIN(ITYP+1,2)),ia,ib,KK,
     X		PTNAME(MATRIX(6,NIX))(1:1),str
655        FORMAT(10X,A7,' reduced matrix element ',I3,' <->',I3,
     X      ' for k=',i2,a1,'  is',F10.4)
      	IC = ICK
	IA = abs(IA)
	I=MATRIX(6,NIX)
          R = RR * CCPT(I) / CC
      	JA=JEX(I,IC,IA)
      	KG=JEX(I+2,IC,IA)
      	JB=JEX(I,IC,IB)
      	PH = (-1)**NINT( (JA-JB + ABS(JA-JB))*0.5 )
      		RC = PH * SQRT(2.*JA+1.) * CLEB6(JA,KG,KK+Z8,Z8,JB,KG)
      	IF(ABS(RC)>1E-5.and.
     X     COPY(I,IC,IB)==0.and.COPY(I,IC,IA)==0) then
       		   BETA = str / (RR * RC)
         	   WRITE(KO,1705) IB,IA,KK,BETA,IC
          	IF(PT(I).AND.BETA.NE.0.0) THEN
                   BETA = str / (R * RC)
             	   WRITE(KO,1715) PTNAME(I),BETA,R
            	ENDIF
		endif
	NIX = NIX+1
            IF(J.GT.0) GO TO 650
660        CONTINUE
	  
	    do 670 JF=MF,NF
               CALL CHECK(JF,MLOC,24)
            HPOT(JF) = H
            PTYPE(1,JF) = PTYPE(1,NFD)
            PTYPE(2,JF) = TYPE
            PTYPE(3,JF) = ICK
	 	if(JF>MF) PTYPE(3,JF) = -1
            PTYPE(4,JF) = NFD
            PTYPE(5,JF) = ITYP
            PTYPE(6,JF) = 0
            IF(itt) PTYPE(6,JF) = 1
            IF(nosub) PTYPE(6,JF) = PTYPE(6,JF)+2
            PTYPE(7,JF) = 0
            PTYPE(8,JF) = 0
	      LDEP(JF) = 0
	      if(JL.eq.'J') LDEP(JF) = 1
	      if(JL.eq.'L') LDEP(JF) = 2
	  LSHAPE(JF) = LSHAPEI
	  VARYL(1,JF) = XLVARY
	  VARYL(2,JF) = ALVARY
670	    continue
C					The monopole is replaced completely!
               PTYPE(4,NFD) = -NF
        GO TO 190

      ELSE IF(TYPE==20 .OR. TYPE==21) THEN
         SHAPE = MAX(1,MIN(12,SHAPE))
         NF = NF + SHAPE
      ELSE
         IF(TYPEI.GE.0) NF = NF + 1
         WRIT25 = SHAPE.GE.10 .AND. SHAPE.LT.20
         IF(WRIT25) SHAPE = SHAPE-10
      ENDIF
         CALL CHECK(NF,MLOC,24)
      IF(TYPE.LE.9.or.TYPE.ge.18) NFD = NF
C
      IF(TYPE.EQ.0) THEN
         PT(2) = ABS(P1) .GT. 1E-5
         PT(1) = ABS(P2) .GT. 1E-5
         CCPT(1) = P2**TH
         CCPT(2) = P1**TH
      ELSE IF(ABS(P0).GT.1E-5) THEN
         PT(2) = .TRUE.
         PT(1) = .FALSE.
         CCPT(2) = P0**TH
         CCPT(1) = 0.0
      ENDIF
         CC = CCPT(1) + CCPT(2)
         A = CC**3
	LIN =KPOP(KP)
	if(LIN>=0) then
	if(TYPE==0) then
	 write(900+KP,889) P1,P2,LIN,P3,P4,LIN,TYPE,LIN,SHAPE
889	 format('  a2=',f8.3,' a1=',f8.3,/,
     &    '  potl(3:4,',i3,') =',2f12.5,12x,' potltype(',i3,')=',i2,
     &    '  potlshape(',I3,')=',I2)
	else
	 write(900+KP,890) LIN,P1,P2,P3,LIN,TYPE,LIN,SHAPE
890	 format('  potl(1:3,',i3,') =',3f12.5,' potltype(',i3,')=',i2,
     &         '  potlshape(',I3,')=',I2)
	endif
	KPOP(KP) = KPOP(KP)+1  ! for next time
	endif
C
17    FORMAT(/1X,I4,I5,'=',A11,1X,A1,I2,'=',A11,I3,6F10.4,2F10.3)
171   FORMAT(14X,'At ',I4,' is k =',I3,:,' (Potl. beta =',F8.5,')')
1715  FORMAT(' ',30X,'      Beta of ',A5,' =',F8.5,
     &       ' (using radius of',F6.3,' fm.).')
172   FORMAT(/41X,'   A#1       A#2        r0c       ac',48X,'h @',I2,
     &/1X,I4,I5,'=',A11,1X,I3,'=',A11,I3,2F10.3,4F10.4,2F10.3,F10.5)
1721  FORMAT(
     &/1X,I4,I5,'=',A11,1X,I3,'=',A11,I3,2F10.3,4F10.4,2F10.3,F10.5)
173   FORMAT(/1X,I4,I5,'=',A11,1X,I3,'=',A11,' at',I3)
174   FORMAT(3X,116('-'))
175   FORMAT(/41X,'   Def1      Def2      Def3   .....',
     &  /1X,I4,I5,'=',A11,1X,A1,I2,'=',A11,I3,6F10.4,2F10.3,'= RADIUS'/)
176   FORMAT(/41X,'  Mn(E1)    Mn(E2)    Mn(E3)  .....',
     &  /1X,I4,I5,'=',A11,1X,A1,I2,'=',A11,I3,6F10.4,2F10.3,'= RADIUS'/)
177   FORMAT(/41X,'  All-order deformations by diagonalisation',
     &  /1X,I4,I5,'=',A11,1X,A1/)
C
      IF(TYPE.EQ.0) THEN
      SHPNAM = 'CHARGE (WS)'
      if(SHAPE>=7) SHPNAM = SHP(SHAPE+1)
      if(final) then
      WRITE(KO,172) ICK,KP,TYPE,WHO(1),0,SHPNAM,NF,(DEF(K),K=1,7),A,H
      else
      WRITE(KO,1721) KP,TYPE,WHO(1),0,SHPNAM,NF,(DEF(K),K=1,7),A
      endif
C                          make type=0 coulomb form factor
      if(SHAPE<7) then
      Z1Z2 = COULCN
      CALL COULPO(RL(2),P3*CC,P4,H4,Z1Z2,N-1)
         RL(1) = RL(2)
C                   these coulomb forms to be multiplied by z1*z2 later
	FORMF(1:N,NF) = RL(1:N)
	else
!	  P5 = 1.0; P6 = 0.0 ! real & imag scaling for external form
	endif
      ELSE IF(TYPE==20.OR.TYPE==21) THEN
         DO 18 K=1,SHAPE
            M = MF + K -1
18       WRITE(KO,173) KP,TYPE,WHO(TYPE+1),K,NNKIND(K),M
      ELSE IF(TYPE.LE.8.or.TYPE.ge.18) THEN
      IF(SHAPE.LE.12) SHPNAM = SHP(SHAPE+1)
      IF(SHAPE.GE.20.and.SHAPE.lt.30) THEN
          SHPNAM = 'READ L FILE'
          OPEN(SHAPE,RECL=72,FORM='FORMATTED',
     X           ACCESS='DIRECT',STATUS='OLD')
	STOP19 = .true.
        ENDIF
      IF(SHAPE.eq.40) THEN
          SHPNAM = 'Parity +1-2'
	  DEF(1:2) = nint(DEF(1:2)); DEF(3:7) = 0
	  FORMDEF(1:7,NF) = DEF(1:7)
        ENDIF
      IF(SHAPE.eq.41) THEN
          SHPNAM = 'L-dep 0-5,a'
	  DEF(1:7) = nint(DEF(1:7))
	  FORMDEF(1:7,NF) = DEF(1:7)
        ENDIF
      IF(SHAPE.eq.42) THEN
          SHPNAM = 'J-dep 0-5,a'
	  DEF(1:7) = nint(DEF(1:7))
	  FORMDEF(1:7,NF) = DEF(1:7)
        ENDIF
      WRITE(KO,17) KP,TYPE,WHO(TYPE+1),cit,SHAPE,SHPNAM,NF,
     X         (DEF(K),K=1,7),A
      ENDIF
C
      DO 19 M=MF,NF
      HPOT(M) = H
      PTYPE(1,M) = KP
      PTYPE(2,M) = TYPE
      PTYPE(3,M) = 0
         IF(TYPE.EQ.20 .OR. TYPE.EQ.21) PTYPE(3,M) = M-MF+1
         IF(TYPE>=5 .and. TYPE<=7) PTYPE(3,M) = 2   ! tensor parts of optical potentials
      PTYPE(4,M) = 0
      PTYPE(5,M) = 0
      IF(SHAPE.GE.20) PTYPE(5,M) = SHAPE
      PTYPE(6,NF) = 0
      PTYPE(7,NF) = 0
      PTYPE(8,NF) = 0
        IF(itt) PTYPE(6,NF) = 1
        IF(nosub) PTYPE(6,NF) = PTYPE(6,NF)+2
	      LDEP(NF) = 0
	      if(JL.eq.'J') LDEP(NF) = 1
	      if(JL.eq.'L') LDEP(NF) = 2
	  LSHAPE(NF) = LSHAPEI
	  VARYL(1,NF) = XLVARY
	  VARYL(2,NF) = ALVARY
19    CONTINUE
195   IF(SHAPE.GE.7.AND.SHAPE.LE.9.and.(TYPE.le.8.or.TYPE.ge.10))
     x            GO TO 70
         IF(TYPE.EQ.0) THEN
            RR = P3 * CC
            GO TO 100
         ELSE IF(TYPE.LE.8.or.TYPE.ge.17) THEN
            RR = (P2*ABS(P1) + P5*ABS(P4))*CC / (ABS(P1)+ABS(P4)+1E-10)
           IF(P3.EQ.0.0) P3 = 1.00
           IF(P6.EQ.0.0) P6 = 1.00
           ENDIF
      DO 65 I=1,N
         R = (I-1)*H4
         IF(I.EQ.1) R = H4*.5
       IF(SHAPE.LE.6) THEN
         RH1 = (R - P2*CC)/P3
         RH2 = (R - P5*CC)/P6
        ENDIF
C                          select type of tensor, derivative & rl/im
      GO TO (199,20,25,30,30,35,35,35,20,65,50,50,50,50,
C            0   1  2  3  4  5  6  7  8  9  10 11 12 13,

     X      55,55,55,65,65,65,40,47,55,55,55,55,55,55,55,55,20),TYPE+1
C           14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30.. = type
      GO TO 65
199    IF(SHAPE<=9) GO TO 70  ! external coul form
C                                                    -volume q=0, rl&im
20    V = - CMPLX(P1 * CURVE(SHAPEI,0,R,RH1,P3),
     &            P4 * CURVE(SHAPEI,0,R,RH2,P6))
        if(I>MATCH-7) V=0.
      GO TO 60
C                                                    -surface q=0 rl&im
25    V =   CMPLX(P1 * CURVE(SHAPEI,1,R,RH1,P3),
     &            P4 * CURVE(SHAPEI,1,R,RH2,P6))
      GO TO 60
C                                                    -vso     q=1 rl&im
30    V =   CMPLX(P1 * CURVE(SHAPEI,1,R,RH1,P3)/P3,
     &            P4 * CURVE(SHAPEI,1,R,RH2,P6)/P6) * CONLS/(4*R)
      GO TO 60
C                                                    +tensor  q=2 rl&im
35    V =   CMPLX(P1 * CURVE(SHAPEI,2,R,RH1,P3),
     &            P4 * CURVE(SHAPEI,2,R,RH2,P6))
      GO TO 60
40    CONTINUE
      CALL  SSCC(R,RL(6),RL(5),RL(7),RL(8),RL(11),RL(12),RL(1),RL(2),
     &             RL(3),RL(4),RL(9),RL(10))
C                                                     super soft-core
      DO 45 K=1,SHAPE
45    FORMF(I,MF+K-1) = RL(K)
      GO TO 65
47    CONTINUE
      CALL  NNPOT(R,RL,SHAPE)
C                                                     users NN potential
      DO 49 K=1,SHAPE
49    FORMF(I,MF+K-1) = RL(K)
      GO TO 65
C
50    IF(SHAPE<=9) GO TO 70
C                                deform previous potential :
         CALL DEFORM(DEF,FORMF,NFD,MF,NF, MAXM,MLOC,PTYPE,N,H,RR,
     &               SHAPE.GE.11,SHAPE.le.12,COULCN,STREN,LAMBDA)
         GO TO 100
55    IF(SHAPE<=9) GO TO 70
C                                second derivative of previous potential :
         CALL DDEFORM(DEF,FORMF,NFD,MF,NF, MAXM,MLOC,PTYPE,N,H,
     &               COULCN,STREN,LAMBDA)
         GO TO 100
C
60    FORMF(I,NF) = FORMF(I,NF) + V
65    CONTINUE
      GO TO 100
C
70    IF(SHNG) REWIND INFORM
      DO 80 M=MIN(MF,NF),NF
         K = PTYPE(3,M)
         IF(K.EQ.0) DF = 1.0
         IF(K.GT.0) DF = DEF(K)
	 READ(INFORM,'(a)') COMMENT
         READ(INFORM,*) NPOINTS,RSTEP,RFIRST
	 datavals= '   real'
	 if(SHAPE==8) datavals= '   imag'
	 if(SHAPE==9) datavals= 'complex'
         WRITE(KO,71) COMMENT,M,NPOINTS,datavals,RSTEP,RFIRST,DF
71      FORMAT('  Input:',a80/
     x   '  Form at',I3,': Reading ',I4,1x,a7,' data points at '
     X   ,F8.4,' intervals, starting from r =',F8.4,', with ',
     X   'overall coefficient',F8.4)
         IF(TYPE>0.and.TYPE.LE.9.or.TYPE.ge.18) WRITE(KO,711) P1,P2
         IF(TYPE==0) WRITE(KO,711) P5,P6
711	FORMAT('    and multiplying real & imag parts by ',2f10.5,
     X	  ' respectively')
	if(NPOINTS>MMX*N) then
	  write(KO,*) ' Only room to read in ',MMX*N,' potential',
     x     ' points, not ',NPOINTS
	  stop
	  endif
      IF(SHAPE.LE.8) READ(INFORM,*) (RL(I),I=1,NPOINTS)
!      IF(SHAPE.EQ.9) READ(INFORM,*) (CO(I),I=1,NPOINTS) ! (r,i) complex format (old)
      IF(SHAPE.EQ.9) then				 ! r,i   complex format
	 READ(INFORM,*) (RL(I),IM(I),I=1,NPOINTS)
	 do I=1,NPOINTS
	  CO(I) = cmplx(RL(I),IM(I))
	  enddo
	endif
232   FORMAT(6E12.4)
      RSTEPI = 1./RSTEP
      DO 80 I=1,N
      R = (I-1)*H
      R1 = 0.0
      R2 = 0.0
         IF(I.le.MATCH) then
*      IF(SHAPE.EQ.7) R1 = RL(I)
*      IF(SHAPE.EQ.8) R2 = RL(I)
*      IF(SHAPE.EQ.9) R1 = real(CO(I))
*C     if(shape.eq.9) r2 = aimag(co(i))
*      IF(SHAPE.EQ.9) R2 =      CO(I)  * (0E0,-1E0)
      IF(SHAPE.EQ.7) R1 = FF4(R-RFIRST,RSTEPI,RL,NPOINTS)
      IF(SHAPE.EQ.8) R2 = FF4(R-RFIRST,RSTEPI,RL,NPOINTS)
      IF(SHAPE.EQ.9) R1 = real(FFC4((R-RFIRST)*RSTEPI,CO,NPOINTS))
      IF(SHAPE.EQ.9) R2 = FFC4((R-RFIRST)*RSTEPI,CO,NPOINTS)*(0.,-1.)
      endif
	if(TYPE==0) then
	  V = CMPLX(R1*P5,R2*P5)
        else if(TYPE.LE.9.or.TYPE.ge.18) then
	  V = CMPLX(R1*P1,R2*P2)
	else
          V = CMPLX(R1,R2) * DF
	endif
80    FORMF(I,M) = FORMF(I,M) + V
      RR = P3 * CC
C			Have potential form. Print it out if required:
100   DO 105 I=1,N
105      RL(I) = (I-1)*H
      IF(TRENEG.GE.1) THEN
          DO 115 M=MIN(MF,NF),NF
         IF(ABS(PTYPE(5,M))+PTYPE(3,M).EQ.0.AND.TRENEG.LT.3)GOTO115
            WRITE(KO,110) M
110         FORMAT(/' Potential Form at',I4,' is')
            WRITE(KO,120) (RL(I),FORMF(I,M),I=1,MATCH,MR)
115       CONTINUE
120         FORMAT(5(1X,F5.1,':',2F9.4,1X))
           ENDIF
          IF(WRIT25) THEN
C            WRITE OUT FORM AT NF TO FILE 25
            IF(IR0.EQ.0) OPEN(25,RECL=72,ACCESS='DIRECT',
     X         FORM='FORMATTED',STATUS='UNKNOWN',FILE='potent25')
            NR = (N-1)/3 + 1
            DO 150 IR=1,NR
               I0 = (IR-1)*3
150            WRITE(25,REC=IR+IR0,FMT=232)
     X           (FORMF(I+I0,NF),I=1,MIN(3,N-I0))
            WRITE(KO,*) 'Form at ',NF,' written to file 25 in records ',
     X              IR0+1,' to ',IR0+NR-1
            IR0 = IR0 + NR
            WRIT25 = .FALSE.
            ENDIF
C
190   if(LDEPP) then
	write(KO,191) JL,LSHAPEI,SHP(LSHAPEI+1),XLVARY,ALVARY
191	format(6x,'  Potential component multiplied by ',
     X         a1,'-dependent factor of shape=',
     x	       i1,'=',a11,':   X =',f8.4,'  A =',f8.4/)
	endif
      IF(KPI.GT.0) GO TO 10
200   NIX = NIX - 1
      IF(IR0.GT.0) WRITE(KO,*) 'FILE 25 needed ',IR0-1,' records'

      if(TRENEG>0) then
	written(34) = .true.
        DO 81 M=NF1,NF
       if(PTYPE(2,M).ne.17) 
     x   write(34,79) M,(PTYPE(I,M),I=1,3)
79     format('# Potential at',i4,': KP=',i3,' of type',i3,
     x			' and multipole =',i3)
         DO I=1,N
         write(34,'(f8.3,1p,2e15.6)') (I-1)*H,FORMF(I,M)
	 enddo
81       write(34,*) '&'
	call rewop(89)
	written(89) = .true.
        WRITE(89,'(''#'',f8.5)') HCM
        DO M=1,NF
        WRITE(89,'(''#'',4i5,f8.5)') M,MATCH,MR,1
        WRITE(89,141) (PTYPE(I,M),I=1,3)
141	format('# Optical potential KP=',i3,': type,multipole =',2i3)
	DO 83 I=1,MATCH,MR
83      WRITE(89,144) RL(I),FORMF(I,M)
144   FORMAT(1X,F8.3,1p,2g13.5)
	WRITE(89,*) '&'
	enddo
      endif
	do KP=1,KPMAX
	 if(KPOP(KP)>0) then
!	   write(900+KP,891)
891	   format('&')
	   close(900+KP)
	  endif
	enddo
      RETURN
      END
      SUBROUTINE DEFORM(DEF,FORMF,NFD,MF,NF, MAXM,MLOC,PTYPE,N,H,RR,
     &                  QUAD,VOLCON,COULCN,STREN,LAMBDA)
	use io
      IMPLICIT REAL*8(A-H,O-Z)
C                               (note that implicit is different here|)
C
C        deform the form factor at NFD into MF-to-NF inclusive.
C        by deformation lengths DEF(k),k=1,7
C        and putting into MF-NF the projections onto multipoles
C
C     If 'QUAD' then do numerical quadrature over angles,
C     otherwise use simple derivative form-factors
C
      COMPLEX*16 FORMF(MAXM,MLOC)
      REAL*8 H,AI(4),R,COULCN,STREN(MLOC)
      REAL*8 PL(8),DEF(0:8),SP(9),W(9),C(8),A4(4)
      INTEGER PTYPE(8,MLOC),LAMBDA(MLOC)
      LOGICAL QUAD,VOLCON
      DATA (SP(I),I=1,5) /-.96816024,-.8360311,-.61337143,-.32425342,0./
     &     ,(W(I),I=1,5) /.08127439,.18064816,.26061070,.31234708,
     &                    .33023936 /
      PI = 4 * ATAN(1.0)
      RSP = 1.0/SQRT(PI)
      SQFPI = SQRT(4.*PI)
      DO 2 I=1,4
      SP(5+I) =-SP(5-I)
2     W(5+I)  = W(5-I)
      IF(QUAD) THEN
      DELTA = 0.0
      DO 15 K=0,7
      IF(VOLCON.and.K.GE.1.AND.RR.GT.0.) 
     X        	DELTA = DELTA + DEF(K)**2 / (4.*PI * RR)
15    C(K+1) = SQRT((2.*K+1.)/4.)
         PL(1) = 1.
      DO 50 NU=1,9
         U = SP(NU)
         PL(2) = U
         DO 20 K=2,7
20       PL(K+1) = ((2*K-1)*U*PL(K-1+1) - (K-1)*PL(K-2+1))/DBLE(K)
      SH  = - DELTA
      DO 25 K=1,7
25       SH = SH + C(K+1)*RSP * PL(K+1) * DEF(K)
         DO 40 I=1,N
            R = (I-1)*H - SH
            CALL SPLINT(R/H,N,II,JM,AI)
            DO 30 J=1,JM
30          A4(J) = AI(J)
        DO 40 M=MF,NF
            K = PTYPE(3,M)
            CNS = W(NU) * PL(K+1) * C(K+1)
        DO 40 J=1,JM
40      FORMF(I,M)=FORMF(I,M) + CNS * A4(J) * FORMF(II+J-1,PTYPE(4,M))
50      CONTINUE
      ELSE
C                use ordinary derivative of form factor :
      DO 70 M=MF,NF
         J = PTYPE(4,M)
         K = PTYPE(3,M)
         TYPEP = PTYPE(5,M)
      IF(K.GT.0) THEN
C                        non-zero multipole
      IF(TYPEP.GT.0) THEN
C                        nuclear
C         SH = - 0.5 * DEF(K) / (H * SQFPI)
C         DO 60 I=2,N-1
C60       FORMF(I,M) = (FORMF(I+1,J) - FORMF(I-1,J)) * SH
         SH = - DEF(K) / (12.0 * H * SQFPI)
         DO 60 I=3,N-2
60       FORMF(I,M) = (FORMF(I-2,J) - FORMF(I+2,J) +
     X             8.*(FORMF(I+1,J) - FORMF(I-1,J)) ) * SH
      ELSE
C                        Coulomb
C
C     DEF(k) for Coulomb deformations
C           = M(Ek) = reduced matrix element in units of e.fm**k
C                     by defn. of Alder & Winter (NOT by Brink&Satchler)
C           = +/- sqrt[ (2.I+1) * B(Ek, I -> I') ]
C
C   so is proportional to square root of B(Ek) for a Transition I -> I'
C   (The sqrt(2I+1) factor makes M(Ek) the same for up & down couplings,
C     more precisely : M(Ek,I->I') = (-1)**(I+k-I') * M(Ek,I'->I),
C     and the phase factor is +1 for electric transitions).
C
C         M(Ek)  =  3 * Z * Beta * RR**k / (4*pi)
C                     * (-1)**k * sqrt(2I'+1) * <I'K k0 | IK>
C
C         M(E2,I->I)  =  3 * Beta * RR**2 / (4*pi)
C                        * Z * (3*K*K-I(I+1)) *
C                        * sqrt((2I+1)/((2I-1)I(I+1)(2I+3))) )
C
C                   for states of spin 'I' with projection 'K'
C
         SH = DEF(K) * SQFPI * COULCN / (2*K+1)
         NU = NINT(RR/H) + 1
         U = SH / RR**(2*K+1)
         DO 62 I=1,NU
            R = (I-1)*H
62       FORMF(I,M) = U * R**K
         DO 63 I=NU+1,N
            R = (I-1)*H
63       FORMF(I,M) = SH / R**(K+1)
         STREN(M) = SH
         LAMBDA(M)= K
      ENDIF
      ELSE
         DO 68 I=1,N
68       FORMF(I,M) = FORMF(I,J)
      ENDIF
70    CONTINUE
      ENDIF
	
      RETURN
      END
      SUBROUTINE DDEFORM(DEF,FORMF,NFD,MF,NF, MAXM,MLOC,PTYPE,N,H,
     &                  COULCN,STREN,LAMBDA)
	use io
      IMPLICIT REAL*8(A-H,O-Z)
C                               (note that implicit is different here|)
C
C        double differentiate the form factor at NFD into MF-to-NF inclusive.
C        by deformation lengths DEF(k),k=1,7
C
C        use simple derivative form-factors
C
      COMPLEX*16 FORMF(MAXM,MLOC)
      REAL*8 H,R,COULCN,STREN(MLOC)
      REAL*8 DEF(0:8)
      INTEGER PTYPE(8,MLOC),LAMBDA(MLOC)
      PI = 4. * ATAN(1.0)
      SQFPI = SQRT(4.*PI)
C                use 2nd derivative of form factor :
      DO 70 M=MF,NF
         J = PTYPE(4,M)
         K = PTYPE(3,M)
         TYPEP = PTYPE(5,M)
      IF(K.GT.0) THEN
C                        non-zero multipole
      IF(TYPEP.GT.0) THEN
C                        nuclear
         SH = - DEF(K) / (H*H * SQFPI)
         DO 60 I=2,N-1
60       FORMF(I,M) = (FORMF(I-1,J) - 2.*FORMF(I,J) + FORMF(I+1,J))*SH
      ELSE
         WRITE(6,*) 'SECOND-ORDER COULOMB DEFORMATIONS NOT IMPLEMENTED!'
	 stop
      ENDIF
      ELSE
C                        zero multipole
         DO 68 I=1,N
68       FORMF(I,M) = FORMF(I,J)
      ENDIF
70    CONTINUE
	
      RETURN
      END
      FUNCTION CURVE(SHAPE,DERIV,R,RH,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (EXPMAX = 75.0)
      INTEGER SHAPE,DERIV
      CURVE = 0.0
      IF(SHAPE.GT.6) RETURN
         RHL = RH
         IF(RHL.GT.EXPMAX) RHL = EXPMAX
         IF(RHL.LT.-EXPMAX)RHL =-EXPMAX
      IF(SHAPE.NE.2) E = EXP(-RHL)
C
C      RH is linear function of R such that d(RH)/dR = 1/A
C
      IF(DERIV.GT.0) GO TO 200
C                                               0th derivative
      GO TO (9,10,11,12,13,14,15,16),SHAPE+2
9     CURVE = 1
C				    fourier-bessel
      if(abs(RH)>1e-10) CURVE = SIN(RH)/RH
      RETURN
10    CURVE = E/(1.+E)
C                                   woods-saxon
      RETURN
11    CURVE = (E/(1.+E))**2
C                                   woods-saxon square
      RETURN
C  12    CURVE = EXP(-RH**2)
12    CURVE = EXP(-MIN(RH**2,EXPMAX))
C                                   gaussian
      RETURN
13    CURVE = E/R
C                                   yukawa
      RETURN
14    CURVE = E
C                                   exponential
      RETURN
15    CURVE = - VCRSC(R)
C                                   Reid Soft-Core, central shape, T=0
      RETURN
16    CURVE = - VCRSC0(R)
C                                   Reid Soft-Core, central shape, T=1
      RETURN
C
200    IF(DERIV.GT.1) GO TO 300
C                                              1st derivative
C                              (normalised to -1 where e=1 is possible)
      GO TO (19,20,21,22,23,24,25),SHAPE+2
19    CURVE = (COS(RH)-SIN(RH)/RH)/RH  
C				    fourier-bessel
      RETURN
20    CURVE = - E/(1.+E)**2 * 4
C                                   woods-saxon
      RETURN
21    CURVE = - E**2/(1.+E)**3 * 8
C                                   woods-saxon square
      RETURN
22    CURVE = - EXP(-RH**2) * RH * 2
C                                   gaussian
      RETURN
23    CURVE = - E/R**2 * (R/A + 1.)
C                                   yukawa
      RETURN
24    CURVE = - E
C                                   exponential
      RETURN
25    CURVE = 0.5 * VLSRSC(R) * (4*R/2.000)
C                                   Reid Soft-Core, spin-orbit shape
      RETURN
C                                               2nd derivative
C                                               (normalised to +1
C                                                where e=1 is possible)
300    GO TO (29,30,31,32,33,34,35,36),SHAPE+2
29    CURVE = 1
 	stop 'd2(FB) not implemented'
	RETURN
30    CURVE = E*(E-1.)/(1.+E)**3 * 8
C                                   woods-saxon
      RETURN
31    CURVE = E**2/(1.+E)**3 * 4
C                                   woods-saxon square
      RETURN
32    CURVE = EXP(-RH**2) * RH * 2
C                                   gaussian
      RETURN
33    CURVE = E/R * (1/A**2 + 2/R**2)
C                                   yukawa
      RETURN
34    CURVE = E
C                                   exponential
      RETURN
35    CURVE = 12.0 * VTRSC(R)
C                                   Reid Soft-Core, tensor shape
      RETURN
36    CURVE = 12.0 * EXP(-RH**2)
C                                   12 * Volume Gaussian e.g. for GPT force
      RETURN
      END
      SUBROUTINE SSCC(R,V3P0,V1P1,V3P1,V3P2,V1D2,V3D2,V1S0,V3S1,V3S3D,
     1V3D1,V3P3F,V3F2)
      IMPLICIT REAL*8(A-H,O-Z)
C                               (note that implicit is different here|)
C****
C****
      REAL*8 PP1(2,2),PP2(2,2),PP3(2,2),PP4(2,2),PP5(2,2),PP6(2,2),
     1 PP7(2,2),MU
      REAL*8 PP8(2),PP9(2),PP10(2),PP11(2),PP12(2),PP13(2),PP14(2),
     1PP15(2),PP16(2),PP17(2)
      REAL*8 VC0(2,2),VCC(2,2),VL2(2,2)
      REAL*8 VLS(2),VQ(2),VT0(2),VTT(2)
      DATA MU/.7/
      DATA VC0/31.389,-10.463,-10.463,3.488/
      DATA VT0/-10.463,3.488/
      DATA PP1/75.653,215.32,375.,144.83/
      DATA PP2/3.,0.85807,0.47552,0.88787/
      DATA PP3/-286.26,-883.6,-1001.6,-241.34/
      DATA PP4/2.0254,3.5042,3.6071,3.3788/
      DATA PP5/15.633,17.626,14.,6.65/
      DATA PP6/2.01,2.6463,2.5,1.965/
      DATA PP7/0.72581,-0.35261,-0.35,-0.959/
      DATA PP8/-40.466,520./
      DATA PP9/5.768,5.661/
      DATA PP10/-40.408,-54.85/
      DATA PP11/4.0676,4.0141/
      DATA PP12/-58.951,36./
      DATA PP13/1.3171,1.0805/
      DATA PP14/395.18,-110./
      DATA PP15/4.3098,3.95290/
      DATA PP16/-3.9904,0./
      DATA PP17/2.4583,0./
      DDEXP(T) = EXP(MIN(100D0,MAX(-100D0,T)))
      YC(T)=DDEXP(-T)/X
      YT(T)=(1.+(3.+3./T)/T)*DDEXP(-T)/X
      YL2(T)=(1.+2./T)*DDEXP(-T)/X**3
      YLS(T)=-(1.+T)*DDEXP(-T)/X**3
      IF(R.NE.0.) GO TO 3
C****
C**** comportement a l origine
C****
C     v3p0=V1p1=V3p1=V3p2=V1d2=V3d2=V1s0=V3S1=V3s3d=V3d1=V3p3f=V3f2=0.
      V3P0=0
      V1P1 = 0
       V3P1 = 0
       V3P2 = 0
       V1D2 = 0
       V3D2 = 0
       V3S3D= 0
       V3D1 = 0
       V3P3F= 0
       V3F2 = 0
       V1S0 = 0
       V3S1 = 0
      RETURN
C****
    3 CONTINUE
      X=MU*R
      DO 5 IT=1,2
      VTT(IT)=PP12(IT)*DDEXP(-R**4/PP13(IT)**2)+
     1(PP14(IT)*YT(PP15(IT)*X)+VT0(IT)*YT(X))*(1.-DDEXP(-R**6))
      VLS(IT)=(PP8(IT)*YLS(PP9(IT)*X)+PP10(IT)*YLS(PP11(IT)*X))
     1*(1.-DDEXP(-R**6))
      VQ(1)=PP16(1)*YL2(PP17(1)*X)*(1.-DDEXP(-R**6))
      VQ(2)=0.
      DO 4 IS=1,2
      VCC(IS,IT)=PP1(IS,IT)*DDEXP(-R**4/PP2(IS,IT)**2)+
     1(PP3(IS,IT)*YC(PP4(IS,IT)*X)+VC0(IS,IT)*YC(X))*(1.-DDEXP(-R**4))
      VL2(IS,IT)=(PP5(IS,IT)*YL2(PP6(IS,IT)*X)+PP7(IS,IT)*YL2(X))
     1*(1.-DDEXP(-R**6))
    4 CONTINUE
    5 CONTINUE
      V1S0=VCC(1,2)
      V3S1=VCC(2,1)
      V3S3D=2.82847*VTT(1)
      V3D1=(VCC(2,1)-2.*VTT(1)-3.*VLS(1)+6.*VL2(2,1)+21.*VQ(1))
      V3P0=(VCC(2,2)+2.*VL2(2,2)-2.*VLS(2)-4*VTT(2)+10*VQ(2))
      V1P1=(VCC(1,1)+2.*VL2(1,1))
      V3P1=(VCC(2,2)+2*VL2(2,2)-VLS(2)+2*VTT(2)-5*VQ(2))
      V3P2=(VCC(2,2)+2*VL2(2,2)+VLS(2)-.4*VTT(2)+VQ(2))
      V1D2=(VCC(1,2)+6*VL2(1,2))
      V3D2=(VCC(2,1)+6*VL2(2,1)-VLS(1)+2*VTT(1)-21.*VQ(1))
      V3P3F=2.9393876*VTT(2)
      V3F2=(VCC(2,2)-1.6*VTT(2)-4.*VLS(2)+12.*VL2(2,2)-4.*VQ(2))
      RETURN
      END
      SUBROUTINE NNPOT(R,V,N)
      REAL*8 R,V(12),VPOT(2,2)
      CHARACTER*3 PNAME
      COMMON/CHOICE/IDPAR
      COMMON/EMANHP/PNAME
c     DO 10 I=1,N
c10    V(I) = 0.0
c     1     2     3       4    5    6    7    8    9      10   11       12
c     1S0,  3S1,  3S-3D,  3D1, 1P1, 3P0, 3P1, 3P2, 3P-3F, 3F2, 1D2, and 3D2
c REID93:               (IDPAR=0 for REID68)
       IDPAR=1
      PNAME='1S0'
       CALL REIDSC(R,'NN',VPOT)
       V(1) = VPOT(1,1)
       if(N.eq.1) return
      PNAME='3C1'
       CALL REIDSC(R,'NN',VPOT)
       V(2) = VPOT(1,1)
       V(3) = VPOT(1,2)
       V(4) = VPOT(2,2)
       if(N.le.4) return
      PNAME='1P1'
       CALL REIDSC(R,'NN',VPOT)
       V(5) = VPOT(1,1)
       if(N.eq.5) return
      PNAME='3P0'
       CALL REIDSC(R,'NN',VPOT)
       V(6) = VPOT(1,1)
       if(N.eq.6) return
      PNAME='3P1'
       CALL REIDSC(R,'NN',VPOT)
       V(7) = VPOT(1,1)
       if(N.eq.7) return
      PNAME='3C2'
       CALL REIDSC(R,'NN',VPOT)
       V(8) = VPOT(1,1)
       V(9) = VPOT(1,2)
       V(10) = VPOT(2,2)
       if(N.le.10) return
      PNAME='1D2'
       CALL REIDSC(R,'NN',VPOT)
       V(11) = VPOT(1,1)
       if(N.eq.11) return
      PNAME='3D2'
       CALL REIDSC(R,'NN',VPOT)
       V(12) = VPOT(1,1)
       if(N.eq.12) return
c     1     2     3       4    5    6    7    8    9      10   11       12
c     1S0,  3S1,  3S-3D,  3D1, 1P1, 3P0, 3P1, 3P2, 3P-3F, 3F2, 1D2, and 3D2
      RETURN
      END

C     ******************************************************************
C
C VCRSC.......CENTRAL INTERACTION REID SOFT CORE, T=0
C VCRSC0......CENTRAL INTERACTION REID SOFT CORE, T=1
C VTRSC.......TENSOR PART
C VLSRSC......LS PART
C
C     ******************************************************************
      FUNCTION VCRSC (R)
      IMPLICIT REAL*8(A-H,O-Z)
      X=0.7*R
      E=EXP(-X)
      VCRSC=E*(-10.463+E*(105.468+E*E*(-3187.8+9924.3*E*E)))/X
      RETURN
      END
C     ******************************************************************
      FUNCTION VCRSC0(R)
      IMPLICIT REAL*8(A-H,O-Z)
      X=0.7*R
      E=EXP(-X)
      VCRSC0=E*(-10.463-E**3*(1650.6-E**3*6484.2))/X
      RETURN
      END
C     ******************************************************************
      FUNCTION VTRSC (R)
      IMPLICIT REAL*8(A-H,O-Z)
      X=0.7*R
      Y=1.0/X
      E=EXP(-X)
      VTRSC=E*Y*(-10.463*(1.+3.*Y*(1.+Y)-3.*Y*(4.+Y)*E**3)+
     1E**3*(351.77-1673.5*E*E))
      RETURN
      END
C     ******************************************************************
      FUNCTION VLSRSC (R)
      IMPLICIT REAL*8(A-H,O-Z)
      X=0.7*R
      E=EXP(-X)
      VLSRSC=E**4*(708.91-2713.1*E*E)/X
      RETURN
      END
C     ******************************************************************
*From martr@sci.kun.nl Wed Jul  6 17:12:04 1994
*Received: from wn1.sci.kun.nl by marie.ph.surrey.ac.uk; Wed, 6 Jul 94 17:11:55 BST
*Received: from wn2.sci.kun.nl by wn1.sci.kun.nl  via wn2.sci.kun.nl [131.174.80.2] with SMTP 
*       id SAA26945 (8.6.9/2.2 for <I.Thompson@ph.surrey.ac.uk>); Wed, 6 Jul 1994 18:14:02 +0200
*Message-Id: <199407061614.SAA26945@wn1.sci.kun.nl>
*To: Ian Thompson <I.Thompson@ph.surrey.ac.uk>
*Subject: Re: Rread93.f fails to compile 
*Date: Wed, 06 Jul 1994 18:14:01 +0200
*From: Mart Rentmeester <martr@sci.kun.nl>
*Status: R
*
*
*Dear dr. Thompson,
*
*This is a new version of Rreid93.f. I hope I have now removed all errors.
*
*regards,
*
*Mart Rentmeester
*
*
*
      SUBROUTINE REIDSC(R,TYPE,VPOT)
************************************************************************
**    Reid soft core phenomenological potentials
**    Version: July 1994 
**    E-mail: thefalg@sci.kun.nl
**    Reference: Stoks et al. Phys.Rev. C49 (1994) June
**             : Roderick V. Reid, Annals of Physics 50 (1968) 411
**
**    IDPAR=0: original 1968 version
**          1: updated version, including one-pion-exchange with
**             neutral-pion and charged-pion masses; coupling F^2=0.075.
**             Tensor potential is regularized to equal zero at R=0
**
**    INPUT:  R       in fermi
**    -----   TYPE    'PP', 'NP', 'PN', or 'NN' (character*2)
**
**    OUTPUT: VPOT    2x2 potential matrix in MeV on LSJ-basis
**    ------
**
**    COMMON/CHOICE/IDPAR  has to be filled beforehand.
**    COMMON/EMANHP/PHNAME has to be filled beforehand.
**           PHNAME is character*3 and contains the name of the phase
**           shift in the spectral notation.
**           - singlets:           1S0  1P1  1D2  1F3  1G4 ...
**           - triplets uncoupled: 3P0  3P1  3D2  3F3  3G4 ...
**           - triplets coupled:        3C1  3C2  3C3  3C4 ...
**             where 3C1 denotes 3S1 - epsilon1 - 3D1 channel
**                   3C2 denotes 3P2 - epsilon2 - 3F2 channel ...
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER PHNAME*3,PHNAM0*3, TYPE*2,TYP0*2
      PARAMETER(IPTMAX=2000)
      INTEGER SPIN
      LOGICAL FIRST, CALC
      REAL*8 VPOT(2,2)
      COMMON/EMANHP/PHNAME
      COMMON/CHOICE/IDPAR
      COMMON/POT9SV/POTMAT(IPTMAX,5),IPT
      DATA ISO0/-1/, R0/-1D0/, PHNAM0/'***'/, TYP0/'XX'/, IDPAR0/-1/
      DATA NOWRIT/1/
      SAVE FIRST,CALC, NCHAN,SPIN,L,J,ISO

      IF(PHNAME.NE.PHNAM0) THEN
        IF(NOWRIT.EQ.1) PHNAM0=PHNAME
        NCHAN=1
        IF(PHNAME(2:2).EQ.'C') NCHAN=2
        IF(PHNAME(1:1).EQ.'1') SPIN=0
        IF(PHNAME(1:1).EQ.'3') SPIN=1
        READ(PHNAME,'(2X,I1)') J
        L=J
        IF(PHNAME.EQ.'3P0') L=1
        IF(NCHAN.EQ.2) L=J-1
        ISO=MOD(SPIN+L+1,2)
      ENDIF
 
      IF(NOWRIT.EQ.1) THEN
***     Potential is not stored in POTMAT
        IF(IDPAR.EQ.0) CALL REID0(R,L,SPIN,J,NCHAN,ISO,TYPE,VPOT)
        IF(IDPAR.EQ.1) CALL REID1(R,L,SPIN,J,NCHAN,ISO,TYPE,VPOT)
        GOTO 13
      ENDIF
 
      IF(IDPAR.NE.IDPAR0 .OR. TYPE.NE.TYP0 .OR. PHNAME.NE.PHNAM0
     .   .OR. ISO.NE.ISO0) FIRST=.TRUE.
      IF(FIRST) THEN
        CALC=.TRUE.
        IPT=1
      ENDIF
 
      IF(.NOT.FIRST .AND. R.LT.R0) THEN
        CALC=.FALSE.
        IPT=1
      ENDIF
      FIRST=.FALSE.
      IDPAR0=IDPAR
      TYP0=TYPE
      PHNAM0=PHNAME
      ISO0=ISO
      R0 = R
      IF(.NOT.CALC) THEN
   11   RH=POTMAT(IPT,1)
        IF(RH.LT.R) THEN
          IPT=IPT+1
          IF(IPT.GT.IPTMAX) GOTO 12
          GOTO 11
        ELSEIF(DABS(RH-R).LT.1D-6) THEN
          RH = POTMAT(IPT,1)
          VPOT(1,1) = POTMAT(IPT,2)
          VPOT(1,2) = POTMAT(IPT,3)
          VPOT(2,1) = POTMAT(IPT,4)
          VPOT(2,2) = POTMAT(IPT,5)
          GOTO 13
        ENDIF
   12   WRITE(6,*)'*** REIDSC: stored grid file differs with input'
        STOP
      ELSE
        IF(IDPAR.EQ.0) CALL REID0(R,L,SPIN,J,NCHAN,ISO,TYPE,VPOT)
        IF(IDPAR.EQ.1) CALL REID1(R,L,SPIN,J,NCHAN,ISO,TYPE,VPOT)
        POTMAT(IPT,1) = R
        POTMAT(IPT,2) = VPOT(1,1)
        POTMAT(IPT,3) = VPOT(1,2)
        POTMAT(IPT,4) = VPOT(2,1)
        POTMAT(IPT,5) = VPOT(2,2)
        IPT=IPT+1
        IF(IPT.GT.IPTMAX) THEN
            WRITE(6,*)'*** REIDSC: POTMAT is too small'
            STOP
        ENDIF
      ENDIF
 
   13 RETURN
      END
************************************************************************
      SUBROUTINE REID0(R,L,SPIN,J,NCHAN,ISO,TYPE,VPOT)
C*       Reference : Roderick V. Reid, Annals of Physics 50 (1968) 411.
C*       Extension to J>2 waves : Day, Phys. Rev. C24 (1981) 1203.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VPOT(2,2),NEUTM
      INTEGER SPIN
      CHARACTER PHNAME*3, TYPE*2
      COMMON/EMANHP/PHNAME

      E1(X) = FDEXP(-X)/X
      E2(X) = FDEXP(-2D0*X)/X
      E3(X) = FDEXP(-3D0*X)/X
      E4(X) = FDEXP(-4D0*X)/X
      E6(X) = FDEXP(-6D0*X)/X
      E7(X) = FDEXP(-7D0*X)/X

**    Nucleon masses, PDG 1990
      PROTM = 938.27231D0
      NEUTM = 939.56563D0
 
      RM = 938.903D0
      H = 10.463D0
      AMU = 0.7D0
      IF(TYPE.EQ.'PP') SCALE=RM/PROTM
      IF(TYPE.EQ.'NN') SCALE=RM/NEUTM
      IF(TYPE.EQ.'NP'.OR. TYPE.EQ.'PN')
     .  SCALE = RM/(2D0*PROTM*NEUTM/(PROTM+NEUTM))
      X = AMU*R
 
      IF(J.GE.6 .AND. SPIN.EQ.0) THEN
**      The singlet higher J partial waves :
        IF(ISO.EQ.0) VPOT(1,1) =
     .               3D0*H*E1(X) - 634.39D0*E2(X) + 2163.4D0*E3(X)
        IF(ISO.EQ.1) VPOT(1,1) =
     .              -H*E1(X) - 1650.6D0*E4(X) + 6484.2D0*E7(X)
      ELSEIF(J.GE.6 .AND. SPIN.EQ.1) THEN
**      The triplet higher J partial waves :
        IF(ISO.EQ.0) THEN
          VC = -H*E1(X) + 105.468D0*E2(X) - 3187.8D0*E4(X) +
     .          9924.3D0*E6(X)
          VT = -H * ((1D0+3D0/X+3D0/X/X)*E1(X) - (12D0/X+3D0/X/X)*E4(X))
     .        + 351.77D0*E4(X) - 1673.5D0*E6(X)
        ELSEIF(ISO.EQ.1) THEN
          VC = H/3D0*E1(X) - 933.48D0*E4(X) + 4152.1D0*E6(X)
          VT = H*((1D0/3D0+1D0/X+1D0/X/X)*E1(X) - (4D0/X+1D0/X/X)*E4(X))
     .         - 34.925D0*E3(X)
        ENDIF
        VLS = 0D0
        IF(NCHAN.EQ.2) GOTO 100
        IF(L.EQ.(J-1)) VPOT(1,1) = VC - 2D0*(J-1)/(2*J+1)*VT
        IF(L.EQ. J)    VPOT(1,1) = VC + 2D0*VT
        IF(L.EQ.(J+1)) VPOT(1,1) = VC - 2D0*(J+2)/(2*J+1)*VT
      ELSEIF(ISO.EQ.0) THEN
**      Singlet isospin 0 :
        IF(PHNAME.EQ.'1P1') THEN
          VPOT(1,1) = 3D0*H*E1(X) - 634.39D0*E2(X) + 2163.4D0*E3(X)
        ELSEIF(PHNAME.EQ.'1F3' .OR. PHNAME.EQ.'1H5') THEN
          VPOT(1,1) = 3D0*H * (E1(X) - 16D0*E4(X))
**      Triplet isospin 0 uncoupled :
        ELSEIF(PHNAME.EQ.'3D2') THEN
          VPOT(1,1) = -220.12D0*E2(X) + 871D0*E3(X) - 3D0*H*
     .               ((1D0+2D0/X+2D0/X/X)*E1(X) - (8D0/X+2D0/X/X)*E4(X))
        ELSEIF(PHNAME.EQ.'3G4') THEN
          VPOT(1,1) = -3D0*H * ( (1D0+2D0/X+2D0/X/X)*E1(X) -
     .                   (16D0+8D0/X+2D0/X/X)*E4(X) ) + 3133.04D0*E6(X)
**      Triplet isospin 0 coupled :
        ELSEIF(PHNAME.EQ.'3C1') THEN
C         VC = -H*E1(X) + 105.468D0*E2(X) - 3187.8D0*E4(X) +
C    .          9924.3D0*E6(X)
C         VT = -H * ((1D0+3D0/X+3D0/X/X)*E1(X) - (12D0/X+3D0/X/X)*E4(X))
C    .        + 351.77D0*E4(X) - 1673.5D0*E6(X)
C         VLS = 708.91D0*E4(X) - 2713.1D0*E6(X)
**        We use the alternate version (to have a correct deuteron)
          VC = -H*E1(X) + 102.012D0*E2(X) - 2915.0D0*E4(X) +
     .          7800.0D0*E6(X)
          VT = -H * ((1D0+3D0/X+3D0/X/X)*E1(X) - (12D0/X+3D0/X/X)*E4(X))
     .        + 163.016D0*E4(X)
          VLS = 251.572D0*E4(X)
        ELSEIF(PHNAME.EQ.'3C3') THEN
          VC = -H*E1(X) - 103.4D0*E2(X) - 419.6D0*E4(X) + 9924.3D0*E6(X)
          VT = -H * ((1D0+3D0/X+3D0/X/X)*E1(X) - (12D0/X+3D0/X/X)*E4(X))
     .        + 351.77D0*E4(X) - 1673.5D0*E6(X)
          VLS = 650D0*E4(X) - 5506D0*E6(X)
        ELSEIF(PHNAME.EQ.'3C5') THEN
          VC = -H*E1(X) + 105.468D0*E2(X) - 3187.8D0*E4(X) +
     .          9924.3D0*E6(X)
          VT = -H * ((1D0+3D0/X+3D0/X/X)*E1(X) - (12D0/X+3D0/X/X)*E4(X))
     .        + 351.77D0*E4(X) -1673.5D0*E6(X)
          VLS = 0D0
        ENDIF
      ELSEIF(ISO.EQ.1) THEN
**      Singlet isospin 1 :
        IF(PHNAME.EQ.'1S0') THEN
          VPOT(1,1) = -H*E1(X) - 1650.6D0*E4(X) + 6484.2D0*E7(X)
        ELSEIF(PHNAME.EQ.'1D2') THEN
          VPOT(1,1) = -H*E1(X) - 12.322D0*E2(X) - 1112.6D0*E4(X) +
     .                6484.2D0*E7(X)
        ELSEIF(PHNAME.EQ.'1G4') THEN
          VPOT(1,1) = -H*E1(X) - 39.025D0*E2(X) + 6484.2D0*E7(X)
**      Triplet isospin 1 uncoupled :
        ELSEIF(PHNAME.EQ.'3P1') THEN
          VPOT(1,1) = -135.25D0*E2(X) + 472.81D0*E3(X) +
     .             H*((1D0+2D0/X+2D0/X/X)*E1(X) - (8D0/X+2D0/X/X)*E4(X))
        ELSEIF(PHNAME.EQ.'3F3') THEN
          VPOT(1,1) = -729.25D0*E4(X) + 219.8D0*E6(X) +
     .             H*((1D0+2D0/X+2D0/X/X)*E1(X) - (8D0/X+2D0/X/X)*E4(X))
        ELSEIF(PHNAME.EQ.'3H5') THEN
          VPOT(1,1) = H * ((1D0+2D0/X+2D0/X/X)*E1(X)-
     .                      (8D0/X+2D0/X/X)*E4(X))
**      Triplet isospin 1 coupled :
        ELSEIF(PHNAME.EQ.'3P0') THEN
          VPOT(1,1) = 27.133D0*E2(X) - 790.74D0*E4(X) + 20662D0*E7(X) -
     .            H*((1D0+4D0/X+4D0/X/X)*E1(X) - (16D0/X+4D0/X/X)*E4(X))
        ELSEIF(PHNAME.EQ.'3C2' .OR. PHNAME.EQ.'3C4') THEN
          VC = H/3D0*E1(X) - 933.48D0*E4(X) + 4152.1D0*E6(X)
          VT = H*((1D0/3D0+1D0/X+1D0/X/X)*E1(X) - (4D0/X+1D0/X/X)*E4(X))
     .         - 34.925D0*E3(X)
          IF(PHNAME.EQ.'3C2') VLS = -2074.1D0*E6(X)
          IF(PHNAME.EQ.'3C4') VLS = -1037.05D0*E6(X)
        ENDIF
      ENDIF
 
 100  IF(NCHAN.EQ.2) THEN
        VPOT(1,1) = VC - 2D0*(J-1)/(2*J+1)*VT + (J-1)*VLS
        VPOT(2,2) = VC - 2D0*(J+2)/(2*J+1)*VT - (J+2)*VLS
        VPOT(1,2) = 6D0*SQRT(real(J)*real(J+1))/(2*J+1)*VT
        VPOT(2,2) = SCALE*VPOT(2,2)
        VPOT(1,2) = SCALE*VPOT(1,2)
        VPOT(2,1) = VPOT(1,2)
      ENDIF
      VPOT(1,1) = SCALE*VPOT(1,1)
 
      RETURN
      END
************************************************************************
      SUBROUTINE REID1(R,L,SPIN,J,NCHAN,ISO,TYPE,VPOT)
C*      This subroutine calculates the Nijmegen version of the REID pot.
C*      It uses regularized (!) Yukawa functions.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VPOT(2,2), A(5,5),B(5,5), PARSPP(5,5),PARSNP(5,5)
      INTEGER SPIN
      CHARACTER PHNAME*3, TYPE*2
      COMMON/EMANHP/PHNAME
      DATA F0PI/0.075D0/, FCPI/0.075D0/, HBC/197.327053D0/, ICAL/0/
      DATA PIOM,PIOMC,PIOMS/134.9739D0,139.5675D0,139.5675D0/
      DATA PARSPP/
     1  .1756084D0,-.1414234D2, .1518489D3,-.6868230D3, .1104157D4
     2,-.4224976D2, .2072246D3,-.3354364D3,-.1989250D1,-.6178469D2
     3, .2912845D2, .1511690D3, .8151964D1, .5832103D2,-.2074743D2
     4,-.5840566D0,-.1029310D2, .2263391D2, .2316915D2,-.1959172D1
     5,-.2608488D1, .1090858D2,-.4374212D0,-.2148862D2,-.6584788D0/
      DATA PARSNP/
     1 -.2234989D2, .2551761D3,-.1063549D4, .1609196D4,-.3505968D1
     2,-.4248612D1,-.5352001D1, .1827642D3,-.3927086D3, .5812273D2
     3,-.2904577D1, .3802497D2, .3395927D0, .8318097D0, .1923895D1
     4, .0913746D0,-.1274773D2, .1458600D3,-.6432461D3, .1022217D4
     5,-.0461640D0, .7950192D1,-.1925573D1, .5066234D2, .83598955D1/
      SAVE A,B
 
      Y(N,M,X) = FDEXP(-N*X)/X
     .         - FDEXP(-M*X)/X*(1D0+(M**2-N**2)*X/(2D0*M))
      YP(N,M,X)= FDEXP(-N*X)/X
     .         - FDEXP(-M*X)/X*(1D0+(M**2-N**2)*M*X/(2D0*N**2))
      W(N,M,X) = FDEXP(-N*X)/X*(1D0/(N*X)+1D0/(N*X)**2)
     .         - FDEXP(-M*X)/X*(1D0/(M*X)+1D0/(M*X)**2)*M**2/N**2
     .         - FDEXP(-M*X)/X*(M**2-N**2)/(2D0*N**2)
      Z(N,M,X) = FDEXP(-N*X)/X*(1D0+3D0/(N*X)+3D0/(N*X)**2)
     .         - FDEXP(-M*X)/X*(1D0+3D0/(M*X)+3D0/(M*X)**2)*M**2/N**2
     .         - FDEXP(-M*X)*(1D0+1D0/(M*X))*M*(M**2-N**2)/(2D0*N**2)
 
      IF(ICAL.EQ.0) THEN
        DO 1 I1=1,5
          DO 1 I2=1,5
            A(I1,I2)=PARSPP(I2,I1)
  1         B(I1,I2)=PARSNP(I2,I1)
        ICAL=1
      ENDIF
 
      X0 = PIOM /HBC * R
      XC = PIOMC/HBC * R
      VSPIS= F0PI*(PIOM/PIOMS)**2*PIOM/3D0*YP(1,8,X0)
      VSPI = F0PI*(PIOM/PIOMS)**2*PIOM/3D0*Y(1,8,X0)
      VTPI = F0PI*(PIOM/PIOMS)**2*PIOM/3D0*Z(1,8,X0)
      IF(TYPE.EQ.'NP'.OR. TYPE.EQ.'PN') THEN
        VSPIS= (4D0*ISO-2D0)*FCPI*(PIOMC/PIOMS)**2*PIOMC/3D0*YP(1,8,XC)
     .       - VSPIS
        VSPI = (4D0*ISO-2D0)*FCPI*(PIOMC/PIOMS)**2*PIOMC/3D0*Y(1,8,XC)
     .       - VSPI
        VTPI = (4D0*ISO-2D0)*FCPI*(PIOMC/PIOMS)**2*PIOMC/3D0*Z(1,8,XC)
     .       - VTPI
      ENDIF
      PIOMM=(PIOM+2D0*PIOMC)/3D0
      X = PIOMM/HBC * R
 
      IF(PHNAME.EQ.'1S0') THEN
        IF(TYPE.EQ.'PP' .OR. TYPE.EQ.'NN') THEN
          VPOT(1,1) = PIOMM * (A(1,1)*Y(2,8,X) + A(1,2)*Y(3,8,X)
     .        + A(1,3)*Y(4,8,X) + A(1,4)*Y(5,8,X) + A(1,5)*Y(6,8,X) )
        ELSEIF(TYPE.EQ.'NP' .OR. TYPE.EQ.'PN') THEN
          VPOT(1,1) = PIOMM * (B(1,1)*Y(3,8,X) + B(1,2)*Y(4,8,X)
     .               + B(1,3)*Y(5,8,X) + B(1,4)*Y(6,8,X) )
        ENDIF
      ELSEIF(PHNAME.EQ.'1D2') THEN
        VPOT(1,1) = PIOMM * (A(2,1)*Y(4,8,X) + A(2,2)*Y(5,8,X)
     .             + A(2,3)*Y(6,8,X) )
      ELSEIF(PHNAME.EQ.'1G4') THEN
        VPOT(1,1) = PIOMM * A(2,4)*Y(3,8,X)
      ELSEIF(PHNAME.EQ.'3P0') THEN
        VPOT(1,1) = PIOMM * (A(3,1)*Y(3,8,X) + A(3,2)*Y(5,8,X)
     .             + A(2,5)*Z(3,8,X)/3D0 )
      ELSEIF(PHNAME.EQ.'3P1') THEN
        VPOT(1,1) = PIOMM * (A(3,3)*Y(3,8,X) + A(3,4)*Y(5,8,X)
     .             + A(3,5)*Z(3,8,X)/3D0 )
      ELSEIF(PHNAME.EQ.'3F3') THEN
        VPOT(1,1) = PIOMM * A(4,5)*Y(3,8,X)
      ELSEIF(PHNAME.EQ.'3C2' .OR. PHNAME.EQ.'3C4') THEN
        VC = PIOMM * (A(4,1)*Y(3,8,X) + A(4,2)*Y(4,8,X)
     .             + A(4,3)*Y(5,8,X) + A(4,4)*Y(6,8,X) )
        VT = PIOMM * (A(5,1)*Z(4,8,X) + A(5,2)*Z(6,8,X) )/3D0
        IF(PHNAME.EQ.'3C2')
     .       VLS = PIOMM * (A(5,3)*W(3,8,X)+A(5,4)*W(5,8,X))
        IF(PHNAME.EQ.'3C4') VLS = PIOMM * A(5,5)*W(3,8,X)
        VPOT(1,1) = VC + (J-1)*VLS - 2D0*(J-1)/(2*J+1)*VT
        VPOT(2,2) = VC - (J+2)*VLS - 2D0*(J+2)/(2*J+1)*VT
        VPOT(1,2) = 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VT
        VPOT(2,1) = 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VT
      ELSEIF(PHNAME.EQ.'1P1') THEN
        VPOT(1,1) = PIOMM * (B(2,1)*Y(3,8,X) + B(2,2)*Y(4,8,X)
     .             + B(2,3)*Y(5,8,X) + B(2,4)*Y(6,8,X) )
      ELSEIF(PHNAME.EQ.'1F3') THEN
        VPOT(1,1) = PIOMM * (B(1,5)*Y(3,8,X) + B(2,5)*Y(5,8,X) )
      ELSEIF(PHNAME.EQ.'3D2') THEN
        VPOT(1,1) = PIOMM * (B(3,1)*Y(3,8,X) + B(3,2)*Y(5,8,X)
     .             + B(3,3)*Z(3,8,X)/3D0 )
      ELSEIF(PHNAME.EQ.'3G4') THEN
        VPOT(1,1) = PIOMM * B(3,4)*Y(3,8,X)
      ELSEIF(PHNAME.EQ.'3C1' .OR. PHNAME.EQ.'3C3') THEN
        VC = PIOMM * (B(4,1)*Y(2,8,X) + B(4,2)*Y(3,8,X)
     .      +B(4,3)*Y(4,8,X)+ B(4,4)*Y(5,8,X) + B(4,5)*Y(6,8,X) )
        VT = PIOMM * (B(3,5)*Z(4,8,X) + B(5,5)*Z(6,8,X) )/3D0
        IF(PHNAME.EQ.'3C1')
     .       VLS = PIOMM *(B(5,1)*W(3,8,X)+B(5,2)*W(5,8,X))
        IF(PHNAME.EQ.'3C3')
     .       VLS = PIOMM *(B(5,3)*W(3,8,X)+B(5,4)*W(5,8,X))
        VPOT(1,1) = VC + (J-1)*VLS - 2D0*(J-1)/(2*J+1)*VT
        VPOT(2,2) = VC - (J+2)*VLS - 2D0*(J+2)/(2*J+1)*VT
        VPOT(1,2) = 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VT
        VPOT(2,1) = 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VT
      ELSEIF(J.GE.5 .AND. SPIN.EQ.0) THEN
        IF(ISO.EQ.1) THEN
          VPOT(1,1) = PIOMM * (A(1,1)*Y(2,8,X) + A(1,2)*Y(3,8,X)
     .        + A(1,3)*Y(4,8,X) + A(1,4)*Y(5,8,X) + A(1,5)*Y(6,8,X) )
        ELSEIF(ISO.EQ.0) THEN
          VPOT(1,1) = PIOMM * (B(2,1)*Y(3,8,X) + B(2,2)*Y(4,8,X)
     .               + B(2,3)*Y(5,8,X) + B(2,4)*Y(6,8,X) )
        ENDIF
      ELSEIF(J.GE.5 .AND. SPIN.EQ.1) THEN
        IF(ISO.EQ.1) THEN
          VC = PIOMM * (A(4,1)*Y(3,8,X) + A(4,2)*Y(4,8,X)
     .               + A(4,3)*Y(5,8,X) + A(4,4)*Y(6,8,X) )
          VT = PIOMM * (A(5,1)*Z(4,8,X) + A(5,2)*Z(6,8,X) )/3D0
        ELSEIF(ISO.EQ.0) THEN
          VC = PIOMM * (B(4,1)*Y(2,8,X) + B(4,2)*Y(3,8,X)
     .        +B(4,3)*Y(4,8,X)+ B(4,4)*Y(5,8,X) + B(4,5)*Y(6,8,X) )
          VT = PIOMM * (B(3,5)*Z(4,8,X) + B(5,5)*Z(6,8,X) )/3D0
        ENDIF
        IF(NCHAN.EQ.1) THEN
          IF(L.EQ.(J-1)) VPOT(1,1) = VC - 2D0*(J-1)/(2*J+1)*VT
          IF(L.EQ. J)    VPOT(1,1) = VC + 2D0*VT
          IF(L.EQ.(J+1)) VPOT(1,1) = VC - 2D0*(J+2)/(2*J+1)*VT
        ELSEIF(NCHAN.EQ.2) THEN
          VPOT(1,1) = VC - 2D0*(J-1)/(2*J+1)*VT
          VPOT(2,2) = VC - 2D0*(J+2)/(2*J+1)*VT
          VPOT(1,2) = 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VT
          VPOT(2,1) = 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VT
        ENDIF
      ENDIF
 
      IF(NCHAN.EQ.1) THEN
        IF(SPIN.EQ.0) THEN
          IF(L.EQ.0) VPOT(1,1) = VPOT(1,1) - 3D0*VSPIS
          IF(L.NE.0) VPOT(1,1) = VPOT(1,1) - 3D0*VSPI
        ELSEIF(L.EQ.J) THEN
          VPOT(1,1) = VPOT(1,1) + VSPI + 2D0*VTPI
        ELSEIF(PHNAME.EQ.'3P0') THEN
          VPOT(1,1) = VPOT(1,1) + VSPI - 2D0*(J+2)/(2*J+1)*VTPI
        ENDIF
      ELSE
        IF(L.EQ.0) VPOT(1,1) = VPOT(1,1)+VSPIS-2D0*(J-1)/(2*J+1)*VTPI
        IF(L.NE.0) VPOT(1,1) = VPOT(1,1)+VSPI -2D0*(J-1)/(2*J+1)*VTPI
        VPOT(2,2) = VPOT(2,2) + VSPI - 2D0*(J+2)/(2*J+1)*VTPI
        VPOT(1,2) = VPOT(1,2) + 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VTPI
        VPOT(2,1) = VPOT(2,1) + 6D0*DSQRT(J*(J+1)*1D0)/(2*J+1)*VTPI
      ENDIF
 
      RETURN
      END
  



************************************************************************
      FUNCTION FDEXP(X)
      IMPLICIT REAL*8(A-Z)
	real*8, intent(IN):: X
      IF(X.LE.-100D0) THEN
        FDEXP=0D0
      ELSE
        FDEXP=EXP(X)
      ENDIF
      RETURN
      END

!***********************************************************************
!     
!    Copyright (c) 2012, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Ian Thompson, I-Thompson@llnl.gov
!     
!    LLNL-CODE-XXXXX All rights reserved.
!
!    Copyright 2012, I.J. Thompson
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
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
!     
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!     
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!        
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!        
!***********************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C          Code to determine the parameters for the                   C
C             global Koning-Delaroche OMP                             C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  Written: 14 October 2004 by Jutta Escher                           C
C  Modified: 18 Dec 2008 by Ian Thompson                              C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C                                                                     C 
C  Reference:                                                         C 
C                                                                     C 
C     A.J. Koning and J.P. Delaroche                                  C 
C     Local and global optical models from 1 keV to 200 MeV           C 
C     Nuclear Physics A 713 (2003) 310                                C 
C                                                                     C 
C                                                                     C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE KDParam(NA,NZ,NTYPE,ENERG,
     #	VV,RV,AV,WV,RW,AW,WD,RD,AD,VSO,RSO,ASO,WSO,WRSO,WASO,RC,D3)
C
	implicit real (A-H,O-Z)
      PARAMETER(INPUTFILE=5)
      PARAMETER(DTOLERANCE=1.E-12)
C
      WRITE(105,*) ' Welcome to KDParamCode'
C
C  INPUT 
C
      NN=NA-NZ
      WRITE(105,*) ' A=',NA,' Z=',NZ,' N=',NN

      ALPHA=1.D0*(NN-NZ)/NA
      IF(NTYPE.EQ.1) THEN
        WRITE(105,*) ' incoming neutron'
      ELSEIF(NTYPE.EQ.2)  THEN
        WRITE(105,*) ' incoming proton'
      ELSE    
        WRITE(105,*) ' mistake'
	stop 'NTYPE'
      ENDIF
C
C  Calculation of relevant paramenters
C
      IF(NTYPE.EQ.1) THEN
        V1=59.30-21.0*ALPHA-0.024*NA
        V2=0.007228-0.00000148*NA
        V3=0.00001994-0.00000002*NA
        V4=0.0000000071
        W1=12.195+0.0167*NA
        W2=73.55+0.0795*NA
        D1=16.0-16.0*ALPHA
        D2=0.0180+0.003802/(1+EXP((NA-156.0)/8.0))
        D3=11.5
        VSO1=5.922+0.0030*NA
        VSO2=0.0040
        WSO1=-3.1
        WSO2=160.0
        EF=-11.2814+0.02646*NA
        RC=0.0
        VCAV=0.0
      ENDIF   ! NTYPE=1 (NEUTRONS)
C
      IF(NTYPE.EQ.2) THEN
        V1=59.30+21.0*ALPHA-0.024*NA
        V2=0.007067+0.00000423*NA
        V3=0.00001729-0.00000001136*NA
        V4=0.0000000071
        W1=14.667+0.009629*NA
        W2=73.55+0.0795*NA
        D1=16.0+16.0*ALPHA
        D2=0.0180+0.003802/(1+EXP((NA-156.0)/8.0))
        D3=11.5
        VSO1=5.922+0.0030*NA
        VSO2=0.0040
        WSO1=-3.1
        WSO2=160.0
        EF=-8.4075+0.01378*NA
        RC=1.198+0.697/(NA**(2.0/3.0))+12.994/(NA**(5.0/3.0))
        VCAV=1.73*NZ/(RC*(NA**0.333333333))
      ENDIF   ! NTYPE=2 (PROTONS)
C
!     WRITE(105,*) ' Relevant parameters for this projectile+target:'
!     WRITE(105,*) '  ALPHA=',ALPHA
!     WRITE(105,*) '  V1=   ',V1
!     WRITE(105,*) '  V2=   ',V2
!     WRITE(105,*) '  V3=   ',V3
!     WRITE(105,*) '  V4=   ',V4
!     WRITE(105,*) '  W1=   ',W1
!     WRITE(105,*) '  W2=   ',W2
      WRITE(105,*) '  D1=   ',D1
!     WRITE(105,*) '  D2=   ',D2
!     WRITE(105,*) '  D3=   ',D3
!     WRITE(105,*) '  VSO1= ',VSO1
!     WRITE(105,*) '  VSO2= ',VSO2
!     WRITE(105,*) '  WSO1= ',WSO1
      WRITE(105,*) '  WSO2= ',WSO2
      WRITE(105,*) '  EF=   ',EF
!     WRITE(105,*) '  VCAV= ',VCAV
      WRITE(105,*) '  RC=   ',RC
      WRITE(105,*) ' '
C
C  Calculation of OMP paramenters for given energy
C
        WRITE(105,*) ' E=',ENERG
C
      DELVC=VCAV*V1*(V2-2.0*V3*(ENERG-EF)+3.0*V4*(ENERG-EF)**2)
      VV=V1*(1.D0-V2*(ENERG-EF)+V3*(ENERG-EF)**2-V4*(ENERG-EF)**3)+DELVC
      WV=W1*(ENERG-EF)**2/((ENERG-EF)**2+W2**2)
      WD=D1*(ENERG-EF)**2*EXP(-D2*(ENERG-EF))/((ENERG-EF)**2+D3**2)
      write(105,*) 'WD:',D1,(ENERG-EF)**2,EXP(-D2*(ENERG-EF)),
     x                     ((ENERG-EF)**2+D3**2)
      VSO=VSO1*EXP(-VSO2*(ENERG-EF))
      WSO=WSO1*(ENERG-EF)**2/((ENERG-EF)**2+WSO2**2)
C
      RV=1.3039-0.4054/(NA**0.333333333)
      AV=0.6778-0.0001487*NA
      RD=1.3424-0.01585*(NA**0.333333333)
      IF(NTYPE.EQ.1) AD=0.5446-0.0001656*NA      
      IF(NTYPE.EQ.2) AD=0.5187-0.0005206*NA
      RSO=1.1854-0.647/(NA**0.333333333)      
      ASO=0.59
C
C Writing out the OMP parameters
C
      WRITE(105,*) ' '
      WRITE(105,*) ' The OMP parameters for this case and energy are:'
      WRITE(105,*) '  ENERG=',ENERG
      WRITE(105,*) '  '
      WRITE(105,*) '  DELVC=',DELVC
      WRITE(105,*) '  VV=   ',VV
      WRITE(105,*) '  WV=   ',WV
      WRITE(105,*) '  WD=   ',WD
      WRITE(105,*) '  VSO=  ',VSO
      WRITE(105,*) '  WSO=  ',WSO
      WRITE(105,*) '  '
      WRITE(105,*) '  RC=   ',RC
      WRITE(105,*) '  RV=   ',RV
      WRITE(105,*) '  AV=   ',AV
      WRITE(105,*) '  RD=   ',RD
      WRITE(105,*) '  AD=   ',AD
      WRITE(105,*) '  RSO=  ',RSO
      WRITE(105,*) '  ASO=  ',ASO
      WRITE(105,*) '  '
C
	RW = RV
	AW = AV
	WRSO = RSO
	WASO = ASO
      RETURN
      END

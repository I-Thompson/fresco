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
       implicit real*8(a-h,o-z)
       parameter (MAXN=100,MAXCH=10)
       real*8 JTOTAL,HCM,ENLAB,JVAL(MAXCH)
       integer PARITY,LVAL(MAXCH),ITC(1,1),PART(MAXCH),EXCIT(MAXCH),
     x		WDISK
       complex*16 psi(MAXN,MAXCH),smat(MAXCH)


       WRITE(17,657) N,HCM,ENLAB,JTOTAL,PARITY
657   FORMAT(I4,2F8.4,F8.1,I3)
      DO 690 IT=1,ITCM
          IF(ABS(MOD(WDISK,2)).EQ.1.AND.IT.NE.ITC(PEL,EXL)) GOTO 690
          DO 670 C=1,NCH
            IF(ITC(PART(C),EXCIT(C)).NE.IT) GOTO 670
            IF(ABS(MOD(WDISK,2)).EQ.1.AND.C.ne.JIN) GOTO 670
          WRITE(17,660) IT,LVAL(C),JVAL(C),JTOTAL,LVAL(JIN),JVAL(JIN),
     X                   SMAT(C)
  660         FORMAT(2I4,2F6.1,I4,F6.1,2F15.10)
              WRITE(17,680) (PSI(I,C),I=1,M)
  670      CONTINUE
  680      FORMAT(1P,6E12.4)
  690  CONTINUE
         WRITE(17,660) -1
       end

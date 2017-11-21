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

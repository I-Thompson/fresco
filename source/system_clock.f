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

C  NON-F90 SUBROUTINES NEEDED ON SOME MACHINES
      FUNCTION SECOND()
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 TARRAY(2),ETIME
      SECOND = 0.0
C IBM-----------------
C     CALL CPTIME(I)
C     SECOND = I/100.0

C SUN/ALPHA---------
!      SECOND = ETIME(TARRAY)

C F90 real-time clock (NOT cpu time!)
 	call system_clock(ic,icr,icm)
	if(icm*icr.ne.0) SECOND = real(ic)/real(icr)

      RETURN
      END

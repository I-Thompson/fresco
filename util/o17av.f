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

      real vr(12)
      complex veff(6),v
      character*80 head

      do 5 i=1,2
      read(5,1) head
1     format(A80)
5     write(6,1) head
      
C10    read(5,*,end=99) r,veff
10    read(5,*,end=99) r,VR
      v = 0.0
      do 15 i=1,6
      veff(i) = cmplx(vr(2*i-1),vr(2*i))
15    v = v + veff(i)
      v = v / 6.
      write(6,20) r,v
20    format(F9.3,2F11.4)
      go to 10
99    stop
      end

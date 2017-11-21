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

C read fam amplitude files from Fresco,
C and plot individual m-dependent cross sections
      real ja,jb,jap,jbp,ma,mja,mb,mjb
      complex f1(200)
      do ist=1,200
      read(5,*,end=99)ja,jb,jap,jbp,nangl
      ni1=int((2*ja+1)*(2*jb+1)*(2*jap+1)*(2*jbp+1))
c     print *,ni1
	faci = 10.0/((2.*ja+1.)*(2.*jb+1.))


      do ith=1,nangl
      read(5,*)angl
      read(5,90)(f1(i),i=1,ni1)
90    format( 6e12.5)
	xs = 0.0
      do i=1,ni1
        xs = xs + abs(f1(i))**2*faci
      enddo
      write(6,60)  angl,xs
60	format(f8.3,g14.5)
      enddo
      enddo
99	stop
      END
 

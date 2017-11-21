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

      real pot(1000)
      iset = 0
1      read(5,*,end=30)
       read(5,*,end=30) np,h,rfirst
!      write(0,*) np,h,rfirst
10	read(5,*) (pot(i),i=1,np)
	iset = iset+1
	fac = 1.0
	do 20 i=1,np
20	write(6,*) (i-1)*h+rfirst,pot(i)*fac
	write(6,*) '&'
	go to 1
30	stop
	end


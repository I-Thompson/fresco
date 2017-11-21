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

	complex*16 f(3)

	OPEN(1,ACCESS='DIRECT',FORM='FORMATTED',RECL=72,
     X			STATUS='NEW')
	IR=0
1      IR=IR+1
       read(5,232,end=300) (F(I),I=1,3)
       write(1,REC=IR,FMT=232) (F(I),I=1,3)
	go to 1
232    FORMAT(6E12.4)
300	stop
	end

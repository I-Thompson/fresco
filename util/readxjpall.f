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

	real jtotal,jtotlast
	real,allocatable:: SIGJ(:),FUSL(:),xsj(:),xstotal(:)
	character par
	
	read(5,*) ITCM,IF1,IF2
	allocate (SIGJ(0:ITCM),FUSL(IF1:IF2))
 	itplot = min(20,ITCM)+1
	allocate (xsj(itplot),xstotal(itplot))
	write(6,5) ITCM,if1,if2,itplot
5	format('# reading ',3i4,' states; print J-dependence to #',i4)

	xsj(:) = 0.0
	jtotlast = -1
	xstotal(:) = 0.

10	read(5,1445,end=20) JTOTAL,par,LVAL,JIN
     X       , (SIGJ(IT),IT=0,ITCM),(FUSL(I),I=IF1,IFT)
 1445 FORMAT(F7.1,A1,I6,I2,10G12.4,:,/(16X,10G12.4))
!	write(7,1445) JTOTAL,par,LVAL,JIN
	
	if(abs(jtotal-jtotlast)>0.1) then ! write out previous
	 if(jtotlast>-.1) then  
 	   write(6,15) jtotlast,xsj(:)
15	   format(f7.1,1x,21g13.4)
	   xstotal(:) = xstotal(:) + xsj(:)
	 endif
	 jtotlast = jtotal
	 xsj(:) = 0.0
	 endif

	xsj(:) = xsj(:) + sigj(0:itplot-1)
	go to 10

20	write(6,25) xstotal(:)
25	format('# Total cross sections =',/'#',7x,21g13.5)
	end

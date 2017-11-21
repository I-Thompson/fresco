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

       real JTOTAL,JIN,SIGJ(0:1000),FUSL(2),LASTJ
	character psign

	LASTJ = -1
	IF1=1
	IFT=1
	ITCM=26
	TOTJ = 0.
	XSLJ = 0.
	FUSJ = 0.

10     read(38,1445,end=99) JTOTAL,PSIGN,LVAL,JIN
     X       , (SIGJ(IT),IT=0,ITCM),(FUSL(I),I=IF1,IFT)
 1445 FORMAT(F7.1,A1,I4,I2,10G12.4,/(14X,10G12.4))
  	write(6,*) JTOTAL,PSIGN,SIGJ(0),FUSL(IF1)

	if(abs(JTOTAL-LASTJ)>.1) then
	 
	  if(LASTJ>0.) then
	  write(56,1446) LASTJ,FUSJ,TOTJ,XSLJ
 1446   FORMAT(f8.1,11G12.4)

	  TOTJ = 0.
	  FUSJ = 0.
	  XSLJ = 0.
	  endif
	endif
	 TOTJ = TOTJ + sigJ(0)
	 FUSJ = FUSJ + FUSL(IF1)

	 do i=1,ITCM
	  XSLJ = XSLJ + SIGJ(i)
	  enddo
	 LASTJ = JTOTAL

	  go to 10
99 	stop
	end

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
!    along with FRESCO. If not, see <http://www.gnu.org/licenses/>.
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

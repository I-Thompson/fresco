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
c
	subroutine lab2cm(TH,A1,A2,A3,A4,ECMI,ECMF,XCM_LAB)
	implicit real*8 (a-h,o-z)
		
	x=(A1*A3*ECMI/((A2*A4)*(ECMF)))**.5
	
	thetalab = TH*pi/180.d0
	coslab=cos(thetalab)
        c2=coslab**2
        s2=1d0-c2
*
c As 0 < TH < 180 deg (at most), we have always sinlab > 0, so it is 
c ok to take only the positive square root of s2.
*
        sinlab=sqrt(s2)
*
        thetacm=thetalab+asin(x*sinlab)
*
        coscm=cos(thetacm)
        thetacm=thetacm*180.d0/pi !switching to degrees
*
        a=1.+x*coscm
        b=abs(a)
	XCM_LAB=b/(1.+x*x+2.*x*coscm)**1.5
!        sigCM=sigLAB*XCM_LAB

	TH = thetacm   ! return cm degrees
	return
	end

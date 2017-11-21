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


!  1: plane wave
!  2: resonance, weighting exp(-i.delta)
!  3: resonance, weighting exp(-i.delta)*sin(delta)

	real bnorms(3),nm(3)
	complex g(3),bin(3),hm,hp,sm
	rmax= 100.0
	h = 0.1
	dk = 0.002
	xk1 = 1.0; 	xk2 = 2.0
!	xk1 = 0.5;  xk2 = 3.0
!	xk1 = 0.5;  xk2 = 5.0
	widk = 0.1
	resk = 1.5
	nk = (xk2-xk1)/dk
	pi = 4.0*atan(1.0)
	
	nr=rmax/h
	bnorms(:) = 0
	do ir=1,nr
	r = ir*h
	
	bin(:)=0.0
	nm(:)= 0.0
	do ik=0,nk
	 xk = xk1 + ik*dk
	 
	 delta = atan(widk/(resk - xk))
	 sm = exp(cmplx(0.,2*delta))
	 if(ir==1) write(20,4) xk,delta*180/pi
	 
	 s = sin(xk*r)
	 c = cos(xk*r)
	 hp = cmplx(c,s)
	 hm = cmplx(c,-s)
	 
	 g(1) = 1.0
	 g(2) = exp(cmplx(0.,-delta))
	 g(3) = g(2) * sin(delta)
	 
	 nm(:) = nm(:) + abs(g(:))**2
	 
	 bin(1) = bin(1) + g(1) * s
	 bin(2) = bin(2) + g(2) * (0.,.5)*(hm - sm * hp)
	 bin(3) = bin(3) + g(3) * (0.,.5)*(hm - sm * hp)
	 enddo
	 bin(:) = bin(:) * dk
	 nm(:) = nm(:) * dk
	 if(ir==1) write(6,*) nm
	 
	 nm(:) = sqrt(2/(pi * nm(:)))
	 bin(:) = bin(:) * nm(:)
	 bnorms(:) = bnorms(:) + abs(bin(:))**2 * h

	 write(1,4) r,real(bin(:))
	 write(2,4) r,abs(bin(:))**2
	 write(4,4) r,abs(bin(:))
	 
	 ab = nm(1) * (cos(xk1*r) - cos(xk2*r))/r
	 write(12,4) r,real(bin(1)),ab
	 write(3,4) r,real(bin(1)),bnorms(1)
4	format(f8.3,6f12.5)
	enddo
	write(6,*) 'Bin norms:',bnorms(:)
	end

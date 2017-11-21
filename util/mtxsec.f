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
C and plot individual m-dependent cross sections for residual m.
      real jix(2),jex(2),mres
      complex f1(400)
      if = 50		! file number for M=0
	jmax=10    	! max |M| to look for
	njmax = 2*jmax+1
      do ist=1,1000
      read(5,*,end=99)jix,jex,nangl,nearfa,enlab

      NSA = 2*JIX(1) + 1  ! projectile
      NJA = 2*JIX(2) + 1  ! target
      NSB = 2*JEX(1) + 1  ! ejectile
      NJB = 2*JEX(2) + 1  ! residual

	ni1 = NSA * NJA * NSB * NJB
	if(ni1>400) stop 'ni1'

	faci = 10.0/((2.*jix(1)+1.)*(2.*jix(2)+1.))  ! factor for cross section (mb)
      do ith=1,nangl
      read(5,*)angl				     ! scattering angle (cm degrees)
      read(5,90)(f1(i),i=1,ni1)
90    format( 6e12.5)

      do ijres=1,njmax
	mres = -jmax + ijres-1
 	if(mres>=0.0) then
	 ifo = if + int(mres)
	 if(ith==1) write(ifo,91) ist,jex(2),mres
91	format('#   Level pair',i4,': residual spin =',f5.1,' m =',f5.1)

c
c  The f1(i) array can in fortran be declared equivalently as f1(NJB,NSB,NJA,NSA),
c      and accessed as f1(IJB,IMB,IJA,IMA)
c

       xs = 0.0
         IAM = 0
      DO 222 IMA = 1,NSA
        MA = -JIX(1) + IMA-1
      DO 222 IJA = 1,NJA
        MJA = -JIX(2) + IJA-1
      DO 222 IMB = 1,NSB
        MB = -JEX(1) + IMB-1
      DO 222 IJB = 1,NJB
        MJB = -JEX(2) + IJB-1
         IAM = IAM + 1
	if(abs(mjb-mres)<.01) xs = xs + abs(f1(IAM))**2*faci
258	format(2g12.4)
222     continue  ! m's
	write(ifo,258) angl,xs
	endif
      enddo  ! mres
      enddo  ! ith

        do ijres=1,njmax
	mres = -jmax + ijres-1
	 ifo = if + int(mres)
	if(mres>=0.0) write(ifo,*) 'END'    ! end-of-dataset marker, for xmgrace
	enddo

      enddo  ! ist
99	stop
      END
 

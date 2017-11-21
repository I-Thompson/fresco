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
c		This short program reads the WDISK=1 or 2 output
c		for FRESCO to file 17,
c		and formats for plotting by xvgr/xmgr
c
c		Run by Unix command:
c			scatwf < fort.17
c					 (or whatever output file 17 used)
c		Produces file:   plot.psi
c
c		Use   xvgr plot.psi       to see real parts (legends correct)
c		 or
c		      xvgr -nxy plot.psi  to see real & imag parts.
c				 (in which case the legend numbers are wrong!)
c
       implicit real*8(a-h,o-z)
       parameter (MAXN=4000)
       real*8 JTOTAL,JVAL,JVIN,psir(MAXN),psii(MAXN)
       integer PARITY
      	CHARACTER*1 PSIGN(3)
      	DATA PSIGN / '-','?','+' /

	iwf = 0
	open(1,file='plot.psi')
1      READ(5,2,end=20) N,HCM,ENLAB,JTOTAL,PARITY
2   	FORMAT(I4,2F8.4,F8.1,I3)
       WRITE(1,3)
3	format('@legend ON')
       WRITE(1,4) ENLAB,JTOTAL,PARITY
       WRITE(6,4) ENLAB,JTOTAL,PARITY
	if(N.gt.MAXN) stop 'MAXN'
4   	FORMAT('#  E =', F8.4,'; J,pi=',F5.1,I3)

10          READ(5,*) IT,LVAL,JVAL,JTOTAL,LVIN,JVIN,SMATR,SMATI
	  if(IT.eq.-1) go to 1
          WRITE(1,12) IT,LVAL,JVAL,JTOTAL,LVIN,JVIN,SMATR,SMATI
12	format('#', 2I4,2F6.1,I4,F6.1,2F15.10)
          WRITE(6,12) IT,LVAL,JVAL,JTOTAL,LVIN,JVIN,SMATR,SMATI
	write(1,14) iwf,ENLAB,JTOTAL,PSIGN(PARITY+2),
     X			IT,LVAL,JVAL,LVIN,JVIN
14  	FORMAT('@legend string',I3,' "E=',F6.2,f5.1,a1,
     x	  ' #',I2,': lj',I3,f5.1,' <',I3,f5.1,'"')
           READ(5,*) (PSIR(I),PSII(I),I=1,N)
      	do 18 i=1,N-1
      	r = (i-1)*hcm
      	WRITE(1,15) r,PSIR(I),PSII(I)
15    	format(f8.4,1p,2e12.3)
18    	continue
      	WRITE(1,*) '&'
	iwf = iwf + 1
         go to 10 

20	stop
      end

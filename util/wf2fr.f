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

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nrmax=10000)
      DIMENSION D(nrmax),RR(nrmax),WF(nrmax)
C
	write(0,*) ' This program creates wf/potl*wf file for FRESCO'
	i = 1
1	read(5,*,err=10,end=10) rr(i),d(i)
	write(1,*) rr(i),d(i)
	if(rr(i).ne.0.) d(i) = d(i)/rr(i)
	write(2,*) rr(i),d(i)
	i = i+1
	go to 1
10	n = i-1
        dr = rr(2) - rr(1)
	write(0,*) 'number of points =',n,', spacing =',dr
	if(n.gt.nrmax) stop 'nrmax1'

	rfr = 25
	hfr = 0.04
	nfr = nint(rfr/hfr)+2
	alpha = 0.306
	if(nfr.gt.nrmax) stop 'nrmax2'
	
        write(6,*)  ' Fresco bound state created by wf2fr'
        write(6,*) nfr,hfr,0
	do 20 i=2,nfr
	r = (i-1)*hfr
	wf(i) = fival(r,rr,d,n,alpha)
20	write(3,*) r,wf(i)
	wf(1) = wf(2)
	write(6,232) (wf(i),i=1,nfr)
	write(6,232) (0.,i=1,nfr)
232	format(1p,6e12.4)
	stop
	end
      function fival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
************************************************************************
*     real 4-point lagrange interpolation routine.
*     interpolates thr function value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) contains
*     the function value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      real*8 fdis(ndm),y1,y2,y3,y4
      dimension xv(ndm)
      if(r.gt.xv(ndm)) go to 9
      do 5 k=1,ndm-2
 5    if(r.lt.xv(k)) go to 6
      k=ndm-2
 6    nst=max(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      fival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      return
 9    fival=fdis(ndm) * exp(alpha*(xv(ndm)-r))
      return
      end

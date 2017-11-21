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

      parameter(mset=500, mlen=1000)
      real x(mlen,mset),y(mlen,mset)
      real xpl(mlen),refs(mset)
      integer n(mset)
      character*30 infile,outfile
      character*80 card
      character first
      
       write(6,*) 'Type names of input and output files (in 2 lines)'
       read(5,1) infile
1	format(a)
       OPEN(2,ACCESS='SEQUENTIAL',FORM='FORMATTED',STATUS='OLD',
     X       FILE=infile)
       read(5,1) outfile
       OPEN(3,ACCESS='SEQUENTIAL',FORM='FORMATTED',STATUS='UNKNOWN',
     X       FILE=outfile)
	write(6,*) 'Read file ',infile,' and write ',outfile
       nset = 1
       nel   = 1
       n(nset) = 0
9     continue
	refs(nset) = nset
10    continue
      read(2,'(a)',err=20,end=30) card
      if(card(20:22).eq.'Lab') then
!	write(21,'(a10)') card(32:42)
	read(card,11,err=20) refs(nset)
11	format(31x,f10.4)
	go to 10
	endif
12     first = card(1:1)
      if(first.eq.'@' .or. first.eq.'#') go to 10
      read(card,*,err=20) x(nel,nset),y(nel,nset)
      n(nset) = nel
      nel = nel + 1
      if(nel.gt.mlen) then
         write(6,*) 'ONLY ',mlen,' values allowed'
         stop
         endif
      go to 10
20    nset = nset+1
      if(nset.gt.mset) then
         write(6,*) 'ONLY ',mset,' sets allowed'
         stop
         endif
      if(n(nset).lt.4) then
         write(6,*) 'NEED at least 4 points/set for splines'
         stop
         endif
      nel = 1
      go to 9
30    if(nel.eq.1) nset = nset -1
      write(6,*) 'Read in ',nset,' sets, each of ',n,' values'
!	write(6,*) 'Ref values:', refs(1:nset)
      write(6,*) 'How many x values do you want to plot?'
      read(5,*) ninc
      if(ninc.eq.0) stop
       write(6,*) 'Enter the ',ninc,' x values to plot:'
       read(5,*) (xpl(i),i=1,ninc)
!       write(6,*) 'Plot at x=', (xpl(i),i=1,ninc)

       write(6,*) 'Scale factor?'
       read(5,*) scale
       write(6,*) 'scaling by ',scale
       rewind 3
	do 100 incl=1,ninc
	xplot = xpl(incl)
	write(3,41) xplot
41	format('# Plotted for x =',f10.4)

	do 80 is=1,nset
	ref = refs(is)

          yplot = fival(xplot,x,y(1,is),n(is))

80       write(3,*) ref,yplot*scale
	write(3,*) '&'
100	continue
        stop
       end
**********************************************************************
      function fival(r,xv,fdis,ndm)
************************************************************************
*     real 4-point lagrange interpolation routine.
*     interpolates thr function value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) contains
*     the function value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      real fdis(ndm),y1,y2,y3,y4
      dimension xv(ndm)
      if(r.gt.xv(ndm)) go to 9
      do 5 k=1,ndm-2
 5    if(r.lt.xv(k)) go to 6
      k=ndm-2
 6    nst=max(k-2,1)
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
 9    fival=fdis(ndm)
      return
      end

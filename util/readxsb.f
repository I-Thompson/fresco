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

*  Version of readxst, with output to fort.8 too

      parameter(meng=100, mangle=1000)
      integer PEL,EXL,PTY(2,2,meng),nex(2),bin(meng)
      real JEX(2,2,meng),ENEX(2,2,meng),width(meng),asum(meng)
      real theta(mangle),xs(mangle,meng),xss(mangle,meng),ens(meng)
      real esum(mangle)
      
      character*80 card
      character first
      character*30 filename,xst,xsec

	write(6,*) 'Enter filename root, and projectile binding energy:'
	read(5,*) filename,BE
	xst = filename(1:lnblnk(filename))//'.xst'
	xsec = filename(1:lnblnk(filename))//'.xsec'
	write(6,1) filename,xst,xsec
1	format(' Filename =',a30,', so read files ',a30,' and ',a30)
	open(3,file=xst,form='formatted')
	open(4,file=xsec,form='formatted')

      read(3,*) 
      read(3,'(3i3,1p,e12.4)') PEL,EXL,NCHAN,ENLAB
      write(6,'(3i3,1p,e12.4)') PEL,EXL,NCHAN,ENLAB
      do IC=1,NCHAN
      read(3,'(2i4)') IC,NEX(IC)
      write(6,'(2i4)') IC,NEX(IC)
        do IA=1,NEX(IC)
        if(IC.eq.PEL.and.IA.eq.EXL) then
          read(3,'(2(f5.1,i3,f8.4),1p,e12.4)')
     X    (JEX(j,IC,IA),PTY(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGT
        else
          read(3,'(2(f5.1,i3,f8.4),1p,e12.4)')
     X    (JEX(j,IC,IA),PTY(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGR
        endif
        enddo
      enddo
      
       nset = 1
       nel   = 1
       n = 0
10    continue
      read(4,'(a)',err=20,end=30) card
      first = card(1:1)
      if(first.eq.'@' .or. first.eq.'#') go to 10
      read(card,*,err=20) theta(nel),xs(nel,nset)
      n = max(n,nel)
      nel = nel + 1
      if(nel.gt.mangle) then
         write(0,*) 'ONLY ',mangle,' angles allowed'
         stop
         endif
      go to 10
20    nset = nset+1
      if(nset.gt.meng) then
         write(0,*) 'ONLY ',meng,' energies allowed'
         stop
         endif
      nel = 1
      go to 10
30    if(nel.eq.1) nset = nset -1
	dtheta = theta(2) - theta(1)
      write(0,*) 'Read in ',nset,' channels, each of ',n,' angles'
      write(0,*) ' dtheta =',dtheta
      ninc = nset-1

C		find energies of bins
	ien = 0
	do 40 IA=2,NEX(1)
           ien = ien+1
	   ens(ien) = enex(1,1,IA)-BE
	if(IA.lt.NEX(1)) then
	 if(enex(1,1,IA).gt.enex(1,1,IA+1)) go to 45
	endif
40	continue
45	continue
	 nen = ien
	 write(6,50)  nen,(ens(ien),ien=1,nen)
50	format(' ',i4,' bins, at',(10f8.4))

	do 55 ien=2,nen
	width(ien) = ens(ien)-ens(ien-1)
	if(width(ien).gt.width(ien-1).and.ien.gt.2) then
	  width(ien) = 2.0*(width(ien)-0.5*width(ien-1))
	endif
55	continue
	width(1) = width(2)
	 write(6,57)  nen,(width(ien),ien=1,nen)
57	format(' ',i4,' widths :',(10f8.4))

	bin(1) = 0
	do 70 ia=2,nex(1)
	   do 60 ien=1,nen
	   if(abs(ens(ien)+BE-enex(1,1,IA)).lt.1e-5) go to 65
60	   continue
	   write(6,*) ' Bin energy for fresco state ',ia,' not found!'
65	  bin(ia) = ien
70	continue
	write(6,80) (bin(ia),ia=1,nex(1))
80	format(' Fresco bins are summed into energy states:',
     x		(/40i3))

          do 85 iang=1,n
	  esum(iang) = 0.0
	  do 85 ien=1,nen
85	  xss(iang,ien) = 0.0
	  do 86 ien=1,nen
86	  asum(ien) = 0.0

	do 100 ia=2,nex(1)
	  if(bin(ia).eq.0) go to 100
          do 90 iang=1,n
90	  xss(iang,bin(ia)) = xss(iang,bin(ia)) + 
     X		xs(iang,ia)/width(bin(ia))
100	continue

c	write(6,*) 'XSS:',(xss(n,ien),ien=1,nen)
c	write(6,*) 'XS:',(xs(n,ia),ia=1,nex(1))
 	pi = 4.0*atan(1.0)
	rad = 180.0/pi
 	do 150 iang=1,n
	  th = theta(iang)/rad
	  wt = 2.*pi*sin(th) / rad
 	do 150 ien=1,nen
	 esum(iang) = esum(iang) + xss(iang,ien) * width(ien)
	 asum(ien ) = asum(ien ) + xss(iang,ien) * dtheta * wt
150	continue
!        write(6,160) (asum(ien),ien=1,nen)
!160	format(' Cross section summed over angles:',(/10f9.1))
	open(7,form='formatted')
	do 170 ien=1,nen
170	write(7,175) ens(ien),asum(ien)
175	format(f8.4,f12.3)
        write(6,*) ' Cross section summed over angles in fort.7'
        write(6,*) ' Cross sections for breakup in fort.8'

	write(8,*) nen,n
	do 190 iang=1,n
	write(8,180) theta(iang)
180	format(' theta=',f12.5)
	do 185 ien=1,nen
185	write(8,186) ens(ien),xss(iang,ien),0.,0.
186	format(f10.5,1p,e12.5,2e13.5)
190	continue
	
	stop
	end

      function lnblnk(s)
      character*30 s
      l=len(s)
      do 1 i=l,1,-1
      if(s(i:i).ne.' ') go to 5
1     continue
      i=0
5     lnblnk=i
      end

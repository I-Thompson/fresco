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
C readxst4

      parameter(meng=200, mangle=1000, mpset=10)
      integer PEL,EXL,PTY(2,2,meng),nex(2),bin(meng,mpset),ipa(meng)
      real JEX(2,2,meng),ENEX(2,2,meng),width(meng,mpset),asum(meng)
      real theta(mangle),xs(mangle,meng),xss(mangle,meng)
      real esum(mangle),sigt(meng,2),ens(meng,mpset)
	logical ialast(meng)     
      character*80 card
      character first,psign(3)
	data psign / '-','?','+' /


      read(5,*) 
      read(5,'(3i3,e12.4,f8.4)') PEL,EXL,NCHAN,ENLAB,BE
      write(1,'(3i3,1p,e12.4,f8.4)') PEL,EXL,NCHAN,ENLAB,BE
	if(NCHAN.gt.2) stop 'NCHAN>2'
      do IC=1,NCHAN
      read(5,'(2i4)') ICI,NEX(IC)
      write(1,'(2i4)') IC,NEX(IC)
        do IA=1,NEX(IC)
          read(5,'(2(f5.1,i3,f8.4),1p,e12.4)') (JEX(j,IC,IA),
     X       PTY(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGT(IA,IC)
	 ialast(IA) = .false.
        enddo
      enddo
      

C		find energies of bins for each partial wave set:
	IA0=2
	IA1=nex(1)
	npset=0
	do 82 ipset=1,10
*	write(1,48) ipset,JEX(1,1,IA0),psign(PTY(1,1,IA0)+2),IA0,IA1
	ien = 0
	do 40 IA=IA0,NEX(1)
           ien = ien+1
	   ens(ien,ipset) = enex(1,1,IA)-BE
	   write(2,*)  ien,ipset,IA, ens(ien,ipset) , enex(1,1,IA),BE
	   ipa(IA) = ipset
	if(IA.lt.NEX(1)) then
	 if(enex(1,1,IA).gt.enex(1,1,IA+1)) go to 45
	endif
40	continue
*	if(ien.eq.0) write(1,*) ' No more found from ',IA0,NEX(1)
	if(ien.eq.0) go to 83
	 IA = NEX(1)
45	continue
	 IA1 = IA
	 ialast(IA1) = .true.
	write(1,48) ipset,JEX(1,1,IA0),psign(PTY(1,1,IA0)+2),IA0,IA1
48	format(/' Partial wave set #',i2,': ',f5.1,A1,' from',i3,' to',i3)
	 nen = ien
	 write(1,50)  nen,(ens(ien,ipset),ien=1,nen)
50	format(' ',i4,' bins, at',(10f8.4))

	do 55 ien=2,nen
	width(ien,ipset) = ens(ien,ipset)-ens(ien-1,ipset)
	if(width(ien,ipset).gt.width(ien-1,ipset).and.ien.gt.2) then
	  width(ien,ipset) = 2.0*(width(ien,ipset)-0.5*width(ien-1,ipset))
	endif
55	continue
	width(1,ipset) = width(2,ipset)
	 write(1,57)  nen,(width(ien,ipset),ien=1,nen)
57	format(' ',i4,' widths :',(10f8.4))
	npset=ipset
60	continue
  	bin(1,ipset) = 0
	do 70 ia=IA0,IA1
	   do 62 ien=1,nen
	   if(abs(ens(ien,ipset)+BE-enex(1,1,IA)).lt.1e-5) go to 65
62	   continue
	   write(1,*) ' Bin energy for fresco state ',ia,' not found!'
65	  bin(ia,ipset) = ien
70	continue
	write(1,80) (bin(ia,ipset),ia=IA0,IA1)
80	format(' Fresco bins are summed into energy states:',
     x		(/40i3))
	 IA0 = IA
82	continue
83	continue

C			Just print out density SIGT cross sections
        write(1,*) ' Energy distribution of cross sections in stdout'
	write(6,165)
165	format('@legend ON'/'@legend x1 0.2'/'@legend y1 0.8')
	IC = 1
	do 210 ia=2,nex(1)
	  ipset = ipa(ia)
	  ien = bin(ia,ipset)
	if(ia.eq.2.or.ialast(ia-1)) then
	write(6,170) ipset-1,ipset,JEX(1,1,ia),psign(PTY(1,1,ia)+2)
170	format('@legend string ',i2,' " Bin set #',i2,': ',f5.1,A1,'"')
	endif
	
	  if(ien.eq.0) go to 210
          SIGT(IA,IC) = SIGT(IA,IC)/width(ien,ipset)
	  write(6,175) ens(ien,ipset),SIGT(IA,IC)
175	format(f8.4,f12.3)
	  if(ialast(ia)) write(6,*) '&'
210	continue

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

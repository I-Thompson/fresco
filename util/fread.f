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
	program fread
      	implicit real*8(a-h,o-z)
	parameter (mp=1000)
	character*80 headng
	real*8 elab(5),jbord(6),jmax,jtmin,jtmax,mass(2),charge(2)
	real*8 enex(2),jex(6),def(0:8),mek,j,dalitz(8)
	real*8 tnt(4,mp),coef(mp),qscale(11)
	character name*8(2),pwf,ex*2(3),cset,rela,iso,ch1
	integer jump(5,3),cp,labe(4),icop(2),band(2),cpot,pp,pade,
     x	   typei,shapei,type,shape,matrix(3),qq
	logical dryin

	ki = 5
	ko = 3
      	read(ki,1005) headng
      	write(ko,1005) headng
 1005 	format(a80)
      	read(ki,1010) hcm,rmatch,rintp,hnl,rnl,centre,hnn,rnn,rmi
      	write(ko,1010) hcm,rmatch,rintp,hnl,rnl,centre,hnn,rnn,rmi
 1010 	format(9f8.3)
      	if (rmatch.lt.0) then
         read(ki,1011) rasym,accrcy,switch,ajswtch
         write(ko,1011) rasym,accrcy,switch,ajswtch
 1011 	format(f8.2,f8.5,f8.2,f8.2)
	endif
	   m = nint(abs(rmatch)/hcm) + 1
	      rmatch = (m-1)*hcm
	     n = m + 1
	      md = m - 6
	   mr= nint(rintp/hcm)
	      if(mr.le.0) mr = 4
	      rintp = mr * hcm
	    nln = (n-1)/mr + 1
	      if(hnl.eq.0.) hnl = hcm
	   mlt = nint( hnl/hcm )
	   nlt = nint( hcm/hnl )
			   hnl = hcm
	      if(nlt.gt.1) hnl = hcm/nlt
	      if(mlt.gt.1) hnl = hcm*mlt
	      nlt = max(nlt,1)
	      mlt = max(mlt,1)
	   hsp = hcm * mlt
	   nlo = nint(rnl/hsp) + 1
	      if(nlo.lt.4) nlo = 4
	   if(hnn.lt.eps) hnn = rintp
	   nnn = (rnn-rmi) / hnn

	maxn = n
	maxnln = nln
	maxnlo = nlo
	maxnnn = nnn

      	read(ki,1030) jtmin,jtmax,absend,dryin,cset,rela,nearfa,
     x         (jump(i,1),jbord(i),i=2,5)
      	write(ko,1030) jtmin,jtmax,absend,dryin,cset,rela,nearfa
 1030 	format(2f4.0,f8.4,l2,1x,a1,1x,a1,i2,5(i4,f4.0))

	lmax1 = jtmax+20

      read(ki,1055) kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc
 1055 format(2i1,f6.3,f8.3,f6.3,i2,3f8.3)
      write(ko,1055) kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc

      read(ki,1070)
     x     ips,it0,it1,itl,pade,iso,nnu,maxl,minl,mtmin,epc,erange,dk
 1070 format(f6.4,i2,i4,2i2,a1,i3,2i4,i2,f6.2,2f8.4)
      write(ko,1070)
     x     ips,it0,it1,itl,pade,iso,nnu,maxl,minl,mtmin,epc,erange,dk
	maxit = 0
	if(pade.gt.0) maxit = abs(it1)
	maxnnu = nnu

      read(ki,1100) chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm
 1100 format(40i2)
      write(ko,1100) chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm

      ic = 1
	it = 1
	mxx = 0
   30 read(ki,1130) name(1),mass(1),charge(1),nex,
     x   pwf,name(2),mass(2),charge(2),qval
      if(mass(1)+mass(2).lt.1e-5) go to 62
      write(ko,1130) name(1),mass(1),charge(1),nex,
     x   pwf,name(2),mass(2),charge(2),qval
 1130 format(a8,2f8.4,i4,a1,1x,a8,2f8.4,f8.4)
	mxx = max(mxx,abs(nex))

      do 60 ia=1,abs(nex)
      read(ki,1150) jex(1),icop(1),band(1),enex(1),
     x         (jex(1+i),i=2,4,2),  cpot,
     x         jex(2),icop(2),band(2),enex(2),
     x         (jex(2+i),i=2,4,2), (ex(i),i=1,3)
      write(ko,1150) jex(1),icop(1),band(1),enex(1),
     x         (jex(1+i),i=2,4,2),  cpot,
     x         jex(2),icop(2),band(2),enex(2),
     x         (jex(2+i),i=2,4,2), (ex(i),i=1,3)
 1150 format(f4.1,2i2,f8.4,2f4.1,i4,2x,f4.1,2i2,f8.4,2f4.1,3a2)
	it = it+1
   60	continue
	go to 30
   62	write(ko,*)
	mxp = ic
	mxpex = it
	
	nf = 0
	ityp = -1
	nix = 1
70    read(ki,72) kpi,typei,it,shapei,(def(k),k=1,7)
72    format(i3,i2,a1,i2,7f8.4)
      if(kpi.eq.0) then
	write(ko,*)
	go to 80
	endif
      write(ko,72) kpi,typei,it,shapei,(def(k),k=1,7)
      type = abs(typei)
      shape = abs(shapei)
	kp = abs(kpi)
      if(type.ge.10.and.type.lt.15) then
               if(shape.eq.0.or.shape.gt.13) shape = 10
c                        remove next line when coul quadrature written.
               if(ityp.eq.0) shape = min(shape,10)
            do 13 k=0,7
            if(k.ne.0 .and. abs(def(k)).lt.1e-9) go to 13
            if(k.eq.0 .and. shape.ne.12.and.shape.ne.13) go to 13
               nf= nf + 1
13        continue
       if(type.lt.12.or.type.ge.15) go to 1350
1301    read(ki,1302) j,(matrix(i),i=2,3),mek
1302     format(4x,3i4,f8.4)
        write(ko,1302) j,(matrix(i),i=2,3),mek
         if(j.eq.0) go to 1350
            matrix(1) = abs(j)
1340        nix = nix + 1
            if(j.gt.0) go to 1301
1350        continue
	else  if(type.eq.8 .or. type.eq.9) then
         shape = max(1,min(12,shape))
         nf = nf + shape
      	else
         if(typei.ge.0) nf = nf + 1	        
	endif

      if(type.le.9.or.type.ge.15) ityp = type
	go to 70
80	mpair = nix-1
	

	msp = 0	
85     read(ki,852) kn1,kn2,ic1,ic2,ini,kind,ch1,nn,l,lmax,sn,iak,j,ib,
     &         kbpot,krpot,e,isc,ipc,nfl,nam,ampl
852   format(2i3,4i2,1x,a1,3i2,f4.1,i2,f4.1,i2,2i3,f8.4,4i3,f8.4)
      if(kn1*ic1*ic2*ini.eq.0) then
	write(ko,*)
	go to 105
	endif
      write(ko,852) kn1,kn2,ic1,ic2,ini,kind,ch1,nn,l,lmax,sn,iak,j,ib,
     &         kbpot,krpot,e,isc,ipc,nfl,nam,ampl
      kn2 = max(kn2,kn1)
	msp = max(msp,kn2)
      if(kind.eq.4.and.nn.ge.2) then
	read(ki,*) dalitz
	write(ko,*) dalitz
	endif
      if(kind.gt.5) then
      	nk = nn
        mpair = max(mpair,nk)
      	if(nk.gt.0) then
	  if(nk.le.mp) then
            read(ki,8626) ((tnt(i,jj),i=1,4),coef(jj),jj=1,nk)
8626  		format(3(4i3,f8.4))
            write(ko,8626) ((tnt(i,jj),i=1,4),coef(jj),jj=1,nk)
	    else
		write(6,*) ' mp in fread should be at least ',nk
	    endif
	endif
	endif
	go to 85
	
  105	cp = 1
  110 read(ki,1220) icto,icfrom,kind,qq,irem,
     x         kpcore,betar,betai,jmax,rmax
      write(ko,1220) icto,icfrom,kind,qq,irem,
     x         kpcore,betar,betai,jmax,rmax
 1220 format(3i4,3i2,2f8.2,2f4.0)
  115 if(icto.eq.0.or.icfrom.eq.0) go to 120
	ip3 = kpcore
	 if(kind.le.2) then
	   do 125 ik=1,qq
      	   read(ki,121)  ib,ia,k,no,kp,a
      	   write(ko,121)  ib,ia,k,no,kp,a
121      	format(4x,5i4,f8.4)
        	if(ib.eq.0) go to 126
	   nf = nf+1
	 
125	   continue
126	   continue
	 else if(kind.le.8) then
8735  	   read(ki,8739) in,ib,ia,kn,a
   	   write(ko,8739) in,ib,ia,kn,a
8739  	   format(4x,4i4,f8.4)
      	    if(in.eq.0) go to 876
	   if(in.gt.0) go to 8735
876	   continue
	 endif
	 if((kind.eq.3.or.kind.eq.4).and.ip3.ge.4) then
         read(ki,232) (qscale(i+1),i=max(0,-qq),abs(qq))
         write(ko,232) (qscale(i+1),i=max(0,-qq),abs(qq))
232   		format(6e12.4)
	 endif
      cp = cp + 1
      if(kind.ne.7.or.abs(irem).ne.2) go to 110
	kind=8
	go to 115
  120  maxcpl = cp-1

	mloc = nf

      read(ki,1255) (elab(i),labe(i),i=1,3),elab(4)
      write(ko,1255) (elab(i),labe(i),i=1,3),elab(4)
 1255 format(3(f8.4,i8),f8.4)

      write(0,1480) mxx,mxp,mxpex,maxit,maxcpl,maxnnu
 1480 format(/' parameters:    mxx    mxp  mxpex  maxit ',
     x'maxcpl maxnnu ' / '   required:',15i7/)
      write(0,1490)
     x maxn,maxnln,maxnlo,msp,lmax1,mpair
 1490 format(/' parameters:   maxn maxnln maxnlo ',
     x'   msp  lmax1  mpair' / '   required:',15i7/)
      write(0,1495) mloc
 1495 format(/' parameters:   mloc (from potentials only)' 
     X	/ '   required:',15i7/)

	stop
	end

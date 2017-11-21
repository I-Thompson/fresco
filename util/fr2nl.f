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
C  Convert old-style Fresco column input to namelist style.
C
	program fsread
      	implicit real*8(a-h,o-z)
	parameter (mp=1000)
	character*80 headng,TMP,MASFIL,line
	real*8 elab(4),jbord(6),jmax,jtmin,jtmax
	real*8 def(0:8),j,triton(8),p(0:8)
	real*8 coef(mp),ips,ebeta(2)
 	real*8 masst,massp,zp,zt,jp,jt,kkp,kkt
	complex qscale(0:11)
	character*8 namep,namet
	character*2 ex(2)
	character cpwf,cset,rela,iso,ch1,citt
	integer jump(6),cp,nlab(3),cpot,pp,pade,pel,exl,pset,
     x	   typei,shapei,type,shape,qq,copyp,bandp,copyt,bandt,ptyp,ptyt
	equivalence (bandp,ptyp),(bandt,ptyt)
        integer chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm,melfil,pcon,tnt(4,mp),plane,
     x        infam,outfam,buttle,vsearch,echan,enodes,eigens
	logical dry,pwf,fexch,ignore,itt,fatal,nosol,psiren,keep,
     X	        rlmasses,pralpha,bloch
c NAMELISTS:
	namelist/fresco/hcm,rmatch,rintp,hnl,rnl,centre,hnn,
     X   rnn,rmin,rsp, rasym,accrcy,switch,ajswtch, sinjmax,plane,
     X   jtmin,jtmax,absend,dry,rela,nearfa,jump,jbord,pset,jset,
     x   kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc,
     x   ips,it0,iter,fatal,iblock,pade,psiren,iso,nnu,maxl,minl,mtmin,
     X     epc,erange,dk, nosol,nrbases,nrbmin,pralpha,pcon,rmatr,ebeta,
     X   chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm,melfil,
     X	 inh,TMP,MASFIL,unitmass,finec,pel,exl,lab,lin,lex,elab,nlab
	namelist/partition/namep,massp,zp,nex,pwf,namet,masst,zt,qval
	namelist/states/ jp,  copyp,ptyp,bandp,ep,kkp,tp, cpot,
     X	jt,copyt,ptyt,bandt,et,kkt,tt, fexch,ignore,keep,infam,outfam
	namelist/pot/ kp,type,itt,shape,p
!	namelist/pot/ kp,type,itt,shape,def,mnep,mnet,ap,at,rc,ac,p,
!     X		 p0,p1,p2,p3,p4,p5,p6,p7,v,r0,rv,vr0,a,av,
!     X		 w,wr0,rw,aw,wa,r0w,vso,rso,rso0,aso,vsoi,rsoi,asoi,
!     X           wd,wdr,wda,wdr0,awd,defp,deft,vd,vdr,vda
	namelist/step/ib,ia,k,str
	namelist/overlap/ kn1,kn2,ic1,ic2,in,kind,ch1,nn,l,lmax,sn,
     &         ia,j,ib,kbpot,krpot,be,isc,ipc,nfl,nam,ampl,keep,
     & 	       dm,nk,er,e
	namelist/dalitz/ triton
	namelist/twont/ tnt,coef
	namelist/coupling/icto,icfrom,kind,ip1,ip2,ip3,p1,p2,jmax,rmax,
     X			  kfrag,kcore
	namelist/inel/ ib,ia,k,no,kp,a
	namelist/cfp/ in,ib,ia,kn,a,keep
	namelist/scale/ qscale

	ki = 5
	ko1 = 1		! fully-formatted old-syle input
	ko2 = 2		! new-style with system namelist output
	ko3 = 6 	! new-style with customised namelist output
	kos = 10 	! new-style with customised namelist output, scratch
	open(ko1,recl=80,form='formatted',delim='apostrophe')
	open(ko2,recl=200,form='formatted',delim='apostrophe')
!	open(ko3,recl=80,form='formatted',delim='apostrophe')
	open(kos,recl=80,form='formatted',delim='apostrophe',   
     X			status='scratch')
	open(11,recl=80,form='formatted',delim='apostrophe')
	eps = 1d-5
	keep = .false.
	rlmasses = .false.
      	read(ki,1005) headng
      	write(ko1,1005) headng
	write(ko3,1001) headng
 1001   format(a80,/,'NAMELIST')
 1005 	format(a80)
      	read(ki,1010) hcm,rmatch,rintp,hnl,rnl,centre,hnn,rnn,rmi,rsp
      	write(ko1,1010) hcm,rmatch,rintp,hnl,rnl,centre,hnn,rnn,rmi,rsp
 1010 	format(10f8.3)
      	if (rmatch.lt.0) then
         read(ki,1011) rasym,accrcy,switch,ajswtch
         write(ko1,1011) rasym,accrcy,switch,ajswtch
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
	rmin = rmi

	maxn = n
	maxnln = nln
	maxnlo = nlo
	maxnnn = nnn

      	read(ki,1030) jtmin,jtmax,absend,dry,cset,rela,nearfa,
     x         (jump(i),jbord(i),i=1,6)
      	write(ko1,1030) jtmin,jtmax,absend,dry,cset,rela,nearfa
 1030 	format(2f4.0,f8.4,l2,1x,a1,1x,a1,i2,7(i4,f4.0))
         jset = 0
         pset = 0
         IF(CSET.EQ.'T') THEN
           JSET = 1
         ELSE IF (LGE(CSET,'0') .AND. LLE(CSET,'9')) THEN
            JSET = ICHAR(CSET) - ICHAR('0')
         ELSE IF(CSET.EQ.'P') THEN
           PSET = 1
         ELSE IF(CSET.EQ.'M'.or.CSET.eq.'N') THEN
           PSET = -1
         ENDIF
	lmax1 = jtmax+20

      read(ki,1055) kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc
 1055 format(2i1,f6.3,f8.3,f6.3,i2,3f8.3)
      write(ko1,1055) kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc

      read(ki,1070) ips,it0,iter,iblock,
     x   pade,iso,nnu,maxl,minl,mtmin,epc,erange,dk,inh,plane,
     x   smallchan,smallcoup,numnode
 1070 format(f6.4,i2,i4,2i2,a1,i3,2i4,i2,f6.2,2f8.4,2i2,2f8.6,i4)
      write(ko1,1070) ips,it0,iter,iblock,
     x   pade,iso,nnu,maxl,minl,mtmin,epc,erange,dk,inh,plane,
     x   smallchan,smallcoup,numnode
        fatal = iter.gt.0
        iter = abs(iter)
        nosol = iter < it0
	psiren = pade.gt.0
	maxit = 0
	if(pade.gt.0) maxit = abs(iter)
	pade=abs(pade)
	maxnnu = nnu
        vsearch=0; echan=0; enodes=0; bloch=.false.; phase=0.

       if(iblock<0.and.iblock>-9) then
        read(ki,1076) eigens,nrbases,nrbmin,buttle,pralpha,pcon,
     x                  meigs,rmatr,ebeta,weak
 1076    format(i1,i3,i4,i1,l1,I2,i2,f6.2,3f8.4)
         if(eigens>0) then
          read(ki,10761) vsearch,echan,enodes,bloch,phase
10761     format(3i4,l2,f6.0)
          endif
        endif

      read(ki,1100) chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm,melfil
 1100 format(40i2)
      write(ko1,1100) chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm,melfil
!	 inh = 0
	 TMP='/tmp'
	 MASFIL = '/opt/PHS1IT/lib/m88lrb'
!				NEW-STYLE CONSTANTS:
!	 unitmass = 1d0
	 unitmass=1.007335d0    ! unless any masses not integer
 	 finec = 137.03599d0	! 1/alpha (fine-structure constant)
!				OLD-STYLE CONSTANTS:
!	 finec=137.5648d0

	 
      ic = 1
	it = 1
	mxx = 0
   30	 read(ki,1130) namep,massp,zp,nex,cpwf,namet,masst,zt,qval
	 pwf = cpwf.eq.'T'.or.cpwf.eq.' '
      if(namep=='        '.or.namet=='        ') go to 62
      write(ko1,1130) namep,massp,zp,nex,cpwf,namet,masst,zt,qval
      write(ko2,nml=partition) 
      write(kos,1131) namep,massp,nint(zp),nex,pwf,
     X                namet,masst,nint(zt),qval
 1130 format(a8,2f8.4,i4,a1,1x,a8,2f8.4,f8.4)
 1131 format(' &PARTITION namep=''',a8,''' massp=',f8.4,' zp=',i3,
     X ' nex=',i2,' pwf=',L1,/
     X       '            namet=''',a8,''' masst=',f8.4,' zt=',i3,
     X       ' qval=',f8.4,'/')
	mxx = max(mxx,abs(nex))
	rlmasses = rlmasses .or.abs(massp-nint(massp))>eps
     X		    	.or.abs(masst-nint(masst))>eps

      do 60 ia=1,abs(nex)
       read(ki,1150)  jp,  copyp,bandp ,ep,  kkp,tp, cpot,
     X		jt,  copyt,bandt ,et,  kkt,tt,  ex, infam,outfam
	fexch = ex(1)(1:1).eq.'T' .or. ex(1)(2:2).eq.'T'
	ignore = ex(2)(1:1).eq.'T' .or. ex(2)(2:2).eq.'T'
       write(ko1,1150) jp,  copyp,bandp ,ep,  kkp,tp, cpot,
     X		jt,  copyt,bandt ,et,  kkt,tt,  ex, infam,outfam
       write(ko2,nml=states) 
       if(copyp==0) write(kos,11511) jp,  bandp ,ep,cpot
       if(copyp/=0) write(kos,11512) copyp,cpot
       if(abs(kkp)+abs(tp)+abs(kkt)+abs(tt)>0
     X    .or.fexch.or.ignore)
     X    write(kos,1152) kkp,tp,kkt,tt,fexch,ignore,infam,outfam
       if(copyt==0) write(kos,1153) jt,  bandt ,et
       if(copyt/=0) write(kos,1154) copyt
 1150 format(f4.1,2i2,f8.4,2f4.1,i4,2x,f4.1,2i2,f8.4,2f4.1,2a2,2i4)
!1151 format(' &STATES jp=',f4.1,' copyp=',i2,' ptyp=',i2,' ep=',f8.4,
!    X		'  cpot=',i3)
11511 format(' &STATES jp=',f4.1,' ptyp=',i2,' ep=',f8.4,'  cpot=',i3)
11512 format(' &STATES copyp=',i2,19x,'  cpot=',i3)
 1152 format('         kkp=',f4.1,' tp=',f4.1,' kkt=',f4.1,' tt=',f4.1,
     X	 '  fexch=',l1,' ignore=',l1,'  infam=',i4,' outfam=',i4)
 1153 format('         jt=',f4.1,' ptyt=',i2,' et=',f8.4,'/')
 1154 format('         copyt=',i2,'/')
	it = it+1
   60	continue
	go to 30
   62	write(ko1,*) 
   	write(ko2,'('' &partition /   ! End of defining partitions''/)')
   	write(kos,'('' &partition /   ! END OF DEFINING PARTITIONS''/)')
	mxp = ic
	mxpex = it
	
	nf = 0
	ityp = -1
	nix = 1
70    read(ki,72) kpi,typei,citt,shapei,(p(k),k=1,7)
72    format(i3,i2,a1,i2,7f8.4)
       if(kpi.eq.0) then
	write(ko1,*)
   	write(ko2,'('' &pot /   ! End of defining potentials''/)')
   	write(kos,'('' &pot /   ! END OF DEFINING POTENTIALS''/)')
	go to 80
	endif
      write(ko1,72) kpi,typei,citt,shapei,(p(k),k=1,7)
	  itt = citt .ne. ' '
      write(ko2,pot)
      if(itt.or.shapei/=0) then
	write(kos,721) kpi,typei,itt,shapei
	else
        write(kos,721) kpi,typei
	endif
      if(maxval(abs(p(1:7)))>999.) then
        write(kos,724) (p(k),k=1,7)
        else if(sum(abs(p(4:7)))>eps) then
        write(kos,722) (p(k),k=1,7)
	else
        write(kos,723) (p(k),k=1,3)
	endif
721    format(' &pot kp=',i2,' type=',i2,:,' itt=',L1,' shape=',i2)
722    format('      p(1:7)=',f9.3,6f8.4,' /')
723    format('      p(1:3)=',f9.3,2f8.4,' /')
724    format('      p(1:7)=',3g12.4,/ '             ',4g12.4,' /')
      type = abs(typei)
      shape = abs(shapei)
	kp = abs(kpi)
      if(type.ge.10.and.type.lt.16) then
               if(shape.eq.0.or.shape.gt.13) shape = 10
c                        remove next line when coul quadrature written.
               if(ityp.eq.0) shape = min(shape,10)
            do 13 k=0,7
            if(k.ne.0 .and. abs(p(k)).lt.1e-9) go to 13
            if(k.eq.0 .and. shape.ne.12.and.shape.ne.13) go to 13
               nf= nf + 1
13        continue
       if(type.lt.12.or.type.ge.18) go to 1350
1301    read(ki,*) ib,ia,k,str
1302     format(4x,3i4,f8.3)
        write(ko1,1302) ib,ia,k,str
        write(ko2,nml=step) 
        write(kos,1303) ib,ia,k,str
1303     format('   &step ib=',i2,' ia=',i2,' k=',i1,' str=',g11.4,' /')
         if(ib.eq.0) go to 1350
1340        nix = nix + 1
            if(ib.gt.0) go to 1301
1350        continue
	else  if(type.eq.20 .or. type.eq.21) then
         shape = max(1,min(12,shape))
         nf = nf + shape
      	else
         if(typei.ge.0) nf = nf + 1	        
	endif

      if(type.le.9.or.type.ge.18) ityp = type
	if(kpi>0) go to 70
80	mpair = nix-1
	

	msp = 0	
85     read(ki,852) kn1,kn2,ic1,ic2,ini,kind,ch1,nn,l,lmax,sn,ia,j,ib,
     &         kbpot,krpot,be,isc,ipc,nfl,nam,ampl
852   format(2i3,4i2,1x,a1,3i2,f4.1,i2,f4.1,i2,2i3,f8.4,4i3,f8.4)
      if(kn1*ic1*ic2*ini.eq.0) then
	write(ko1,*)
   	write(ko2,'('' &overlap /   ! End of defining overlaps''/)')
   	write(kos,'('' &overlap /   ! END OF DEFINING OVERLAPS''/)')
	go to 109
	endif
	er=0; dm=0; e=0; nk=0
      write(ko1,852) kn1,kn2,ic1,ic2,ini,kind,ch1,nn,l,lmax,sn,ia,j,ib,
     &         kbpot,krpot,be,isc,ipc,nfl,nam,ampl
      write(ko2,overlap)
      write(kos,8521) kn1,kn2,ic1,ic2,ini
	if(ch1/=' '.or.lmax>0.or.ia>0.or.ib>0) then
      	  write(kos,8523) kind,ch1,nn,l,lmax,sn,ia,j,ib
	else
      	  write(kos,8524) kind,nn,l,sn,j
	endif
	if(kbpot/=0.or. nfl/=0.or.nam/=0.or.abs(ampl)>eps) then
          write(kos,8525) kbpot,krpot,be,isc,ipc,nfl,nam,ampl
	 else
          write(kos,8526) kbpot,be,isc,ipc
	 endif
8521   format(' &OVERLAP kn1=',i3,' kn2=',i3,' ic1=',i1,' ic2=',i1,
     X	 ' in=',i2)
8523   format('    kind=',i1,' ch1=''',a1,''' nn=',i2,' l=',i1,
     X   ' lmax=',i1,' sn=',f4.1,' ia=',i2,' j=',f4.1,' ib=',i2)
8524   format('          kind=',i1,' nn=',i2,' l=',i1,' sn=',f4.1,
     X		' j=',f4.1)
8525   format('    kbpot=',i2,' krpot=',i2,' be=',f8.4,' isc=',i2,
     X   ' ipc=',i1,' nfl=',i3,' nam=',i3,' ampl=',f8.4,' /')
8526   format('    kbpot=',i2,' be=',f8.4,' isc=',i1,' ipc=',i1,' /')
      kn2 = max(kn2,kn1)
	msp = max(msp,kn2)
      if(kind.eq.4.and.nn.ge.2) then
	read(ki,*) triton
	write(ko1,*) triton
	write(ko2,dalitz)
	write(kos,dalitz)
	endif
      if(kind.gt.5) then
      	nk = nn
        mpair = max(mpair,nk)
      	if(nk.gt.0) then
	  if(nk.le.mp) then
            read(ki,8626) ((tnt(i,jj),i=1,4),coef(jj),jj=1,nk)
8626  		format(3(4i3,f8.4))
            write(ko1,8626) ((tnt(i,jj),i=1,4),coef(jj),jj=1,nk)
            write(ko2,twont)
	     write(kos,'('' &twont'')')
             write(kos,8627) (jj,(tnt(i,jj),i=1,4),jj,coef(jj),jj=1,nk)
 8627  	     format('   tnt(1:4,',i3,')=',4i4,' coef(',i3,')=',f8.4)
  	     write(kos,'('' /'')')
	    else
		write(0,*) ' mp in fsread should be at least ',nk
	    endif
	endif
	endif
	go to 85
	
  109	cp = 1
  110 read(ki,1220) icto,icfrom,kind,qq,irem,
     x         kpcore,betar,betai,jmax,rmax
	if(kind>8) kind=kind-8
      write(ko1,1220) icto,icfrom,kind,qq,irem,
     x         kpcore,betar,betai,jmax,rmax
 1220 format(3i4,3i2,2f8.2,2f4.0)
	ip1 = qq
	ip2 = irem
	ip3 = kpcore
	p1 = betar
	p2 = betai
	kfrag =0
	kcore =0
      write(ko2,coupling)
      if(abs(p1)+abs(p2)+abs(jmax)+abs(rmax)>eps) then
          write(kos,112) icto,icfrom,kind,ip1,ip2,ip3,p1,p2,jmax,rmax
	else
          write(kos,113) icto,icfrom,kind,ip1,ip2,ip3
	endif
112   format(' &COUPLING icto=',i2,' icfrom=',i2,' kind=',i1,' ip1=',i2,
     X	' ip2=',i2,' ip3=',i2,/,'   p1=',f10.4,' p2=',f8.4,' jmax=',f6.1
     X  ,' rmax=',f6.2,'/')
113   format(' &COUPLING icto=',i2,' icfrom=',i2,' kind=',i1,' ip1=',i2,
     X	' ip2=',i2,' ip3=',i2,' /')

  115 if(icto.eq.0.or.icfrom.eq.0) go to 120
	 if(kind>1.and.kind.le.8) then
8735  	   read(ki,8739) in,ib,ia,kn,a
   	   write(ko1,8739) in,ib,ia,kn,a
	   write(ko2,nml=cfp)
   	   write(kos,8740) in,ib,ia,kn,a
8739  	   format(4x,4i4,f8.4)
8740  	   format('   &cfp  in=',i2,' ib=',i3,' ia=',i3,' kn=',i3,
     X		  '  a=',g12.4,' /')
      	    if(in.eq.0) go to 876
	   if(in.gt.0) go to 8735
876	   continue
	 endif
	 if((kind.eq.3.or.kind.eq.4).and.ip3.ge.10) then
         read(ki,232) (qscale(i),i=max(0,-qq),abs(qq))
         write(ko1,232) (qscale(i),i=max(0,-qq),abs(qq))
232   		format(6e12.4)
         write(ko2,scale)
	 write(kos,880)
 880	 format('   &Scale ')
	 do i=max(0,-qq),abs(qq)
	 write(kos,881) i,qscale(i)
	 enddo
 881	 format('          qscale(',i1,') = (',g12.4,',',g12.4,') ')  
	 write(kos,882) 
 882	 format('   / ')
	 endif
      cp = cp + 1
!      if(kind.ne.7.or.abs(irem).ne.2) go to 110
!	kind=8
!	go to 115
	go to 110
  120  maxcpl = cp-1
	pel = abs(icfrom)
        exl = kind
        lab = ip1
        lin = ip2
        lex = ip3

	mloc = nf
         write(kos,'('' ! *******  END OF FRESCO INPUTS *******  '')')

      read(ki,1255) (elab(i),nlab(i),i=1,3),elab(4)
      write(ko1,1255) (elab(i),nlab(i),i=1,3),elab(4)
 1255 format(3(f8.4,i8),f8.4)
	write(ko2,
     X'(// ''**** now move the remaining lines to the beginning ***'')')

        if(rlmasses) unitmass=1d0

	write(ko2,1500) headng
 1500   format(/a80,/,'NAMELIST')
 	write(11,99) hcm,rmatch,rintp,hnl,rnl,centre,hnn,
     X   rnn,rmin,rsp, rasym,accrcy,switch,ajswtch, sinjmax,
     X   jtmin,jtmax,absend,dry,rela,nearfa,jump,jbord,pset,jset,
     x   kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc,
     x   ips,it0,iter,fatal,iblock,pade,psiren,iso,nnu,maxl,minl,mtmin,
     X        epc,erange,dk,plane, nosol,
     X   chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x        lampl,veff,kfus,wdisk,bpm,melfil,
     X	 inh,TMP,MASFIL,unitmass,finec,pel,exl,lab,lin,lex,elab,nlab
99	format(6g12.4)
 
 	write(ko2,nml=fresco)

	write(ko3,100)hcm,rmatch,rintp
100	format(' &FRESCO','  hcm=',f6.3,' rmatch=',f8.3,' rintp=',f6.2)
	if(nlo>0) write(ko3,101)hnl,rnl,centre
101     format('     hnl=',f6.3,' rnl=',f6.2,' centre=',f6.2)
	if(nnn>0) write(ko3,102)hnn,rnn,rmin
102     format('     hnn=',f6.3,' rnn=',f6.2,' rmin=',f6.2)
	if(rsp>0) write(ko3,1025)rsp
1025    format('     rsp=',f7.3)
	if(abs(rasym)>eps) 
     x	 write(ko3,103) rasym,accrcy,switch,ajswtch,sinjmax
103     format('     rasym=',f9.2,' accrcy=',f8.4,' switch=',f8.2,
     X		' ajswtch=',f8.2,' sinjmax=',f8.0)
	write(ko3,104) jtmin,jtmax,absend,dry,rela,nearfa
104     format('     jtmin=',f5.1,' jtmax=',f8.1,' absend=',f8.4,
     X         ' dry=',L1,' rela=''',a1,''' nearfa=',I2)
	if(jump(1)/=0) 
     X		write(ko3,105) (jump(i),i=1,6),(jbord(i),i=1,6)
105	format('      jump(1:6)=',6i5,' jbord=',6f7.1)
	if(pset/=0.or.jset/=0) write(ko3,107) pset,jset
107     format('     pset=',i3,'  jset=',i3)
    	write(ko3,200) thmin,thmax,thinc,koords,kqmax,pp,cutl,cutr,cutc
200	format('     thmin=',f6.2,' thmax=',f6.2,' thinc=',f6.2,
     X    ' koords=',i1,' kqmax=',i1,' pp=',i1,/,
     X         '     cutl=',f6.2,' cutr=',f6.2,' cutc=',f6.2)
	write(ko3,201) ips,it0,iter,fatal,iblock,pade,nosol,psiren
201	format('     ips=',f7.4,'  it0=',i1,' iter=',i3,' fatal=',l1,
     X	   ' iblock=',i2,' pade=',i1,' nosol=',l1,' psiren=',l1)
	write(ko3,202) iso,nnu,maxl,minl,mtmin,epc
202	format('     iso=''',a1,''' nnu=',i3,' maxl=',i4,' minl=',i3,
     X	   ' mtmin=',i2,' epc=',f7.4)
	if(abs(erange)+abs(dk)>eps) write(ko3,203) erange,dk
203     format('     erange=',f8.4,'  dk=',f8.4)
	if(plane>0) write(ko3,2031) plane
2031    format('     plane=',i3)
        if(iblock<0) write(ko3,2035) nrbases,nrbmin,buttle,pralpha,pcon,
     x				rmatr,ebeta,weak,meigs
2035	format('     nrbases=',i2,' nrbmin=',i2,' buttle=',i2,
     x      ' pralpha=',l1,' pcon=',i1,'  rmatr=',f8.2,'  ebeta=',2f8.4,
     x      ' weak =',1p,e12.2,' meigs =',i4)
	if(eigens>0) write(ko3,2036) eigens, vsearch,echan,enodes,
     x                               bloch,phase
2036	format('     eigens=',i2,' vsearch=',i3,' echan =',i3,
     x      ' enodes=',i2,' bloch=',l1,'  phase=',f8.2)
	write(ko3,204) chans,listcc,treneg,cdetr,smats,xstabl,nlpl,
     x        waves,lampl,veff,kfus,wdisk,bpm,melfil
204	format('     chans=',i2,' listcc=',i2,' treneg=',i2,
     X	       ' cdetr=',i2,' smats=',i2,' xstabl=',i2,' nlpl=',i2,/,
     x         '     waves=',i2,' lampl=',i2,' veff=',i2,' kfus=',i2,
     X         ' wdisk=',i2,' bpm=',i2,' melfil=',i2)
	write(ko3,205) unitmass,finec,inh,pel,exl,lab,lin,lex
205	format('     unitmass=',f9.6,' finec=',f9.5,:,/,
     X         '     inh=',i1,' pel=',i1,' exl=',i2,' lab=',i1,
     X         ' lin=',i1,' lex=',i2)
!	TMP=''''//trim(TMP)//''''
!	write(ko3,206) TMP
!206	format('     TMP=',a80)
!	MASFIL=''''//trim(MASFIL)//''''
!	write(ko3,207) MASFIL
!207	format('     MASFIL=',a80)
	if(sum(abs(elab(2:4)))+sum(abs(nlab))>eps) then
	write(ko3,208) elab,nlab
	else
	write(ko3,209) elab(1)
	endif
208	format('     elab(1)=',f11.4,3f8.4,' nlab=',3i3,' /'/)
209	format('     elab(1)=',f11.4,' /'/)


	rewind(kos)
	do 20 i=1,999999
	  read(kos,1005,end=21) line
	  write(ko3,1005) line
20	continue
21	continue

	write(ko3,1470) 
1470 	format(/' Please check that your input numbers have sufficient'
     X /' range and accuracy in the formats chosen for the namelists,',
     X /' and that ''unitmass'' is correct for your masses.'/)

	write(ko3,1471) 0.0481960    
1471 	format(/' To exactly reproduce previous FRESCO runs using',
     X   ' 2*amu/hbar^2 = ',f9.7,/,
     X   '  you will need to enter to this version the numbers:')
	unitmass=1.007335d0; FINEC=137.5648d0
	write(ko3,205) unitmass,finec
	write(ko3,1471) 0.0478326 
	unitmass=0.999740d0; FINEC=137.0455d0
	write(ko3,205) unitmass,finec
	write(ko3,*) 

      write(ko1,1480) mxx,mxp,mxpex,maxit,maxcpl,maxnnu
 1480 format(/' parameters:    mxx    mxp  mxpex  maxit ',
     x'maxcpl maxnnu ' / '   required:',15i7/)
      write(ko1,1490)
     x maxn,maxnln,maxnlo,msp,lmax1,mpair
 1490 format(/' parameters:   maxn maxnln maxnlo ',
     x'   msp  lmax1  mpair' / '   required:',15i7/)
      write(ko1,1495) mloc
 1495 format(/' parameters:   mloc (from potentials only)' 
     X	/ '   required:',15i7/)

	write(0,*) ' fort.1 has fully-formatted old-syle input'
	write(0,*) ' fort.2 has new-style system namelist output'
	write(0,*) ' fort.6 has new-style customised namelist output'
	end


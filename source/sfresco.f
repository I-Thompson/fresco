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

	program sfresco
	use parameters
	use factorials
	use drier
	use io
	use searchpar
	use searchdata
	use fresco1, only:lmax,jbord,jump,thmin,thmax,thinc,elab,rterms,
     x    peli=>pel,exli=>exl,labi=>lab,lini=>lin,lexi=>lex
	implicit real*8(a-h,o-z)
	integer kp,pline,col,nafrac,channel,q,par,nparameters,ndatasets,
     x          dataset,term,type,maxrank,points,pel,exl
	real*8 jtot,energy,afrac,potential,datanorm,width
	real*8 value,step,valmin,valmax,srch_error(mvars),esave(maxen),
     x         stepi,error,plist(1)
	character*50 input_file,output_file,data_file,search_file,
     x       dat_file(mds),plot_file,tag
	character*3 cmlab
	character*8 stype(0:4),errortype
	character*10 name
	character*100 cmd
        character psign(3)
	logical xvals,abserr,lab,undef,noerror,ranks(0:4),nopot,
     x		nodat,op,loge
	external fcn,futil
        data psign / '-','?','+' /

	namelist /variable/ name,kind,step,valmin,valmax, nul,
     x                      kp,pline,col,potential, dataset,datanorm,
     x     	            nafrac,afrac, energy,jtot,par,channel,width,
     x			    term,nopot
	namelist /data/ type,data_file,points,xmin,delta,lab,energy,
     X          idir,iscale,abserr,ic,ia,k,q,angle,jtot,par,channel,
     x          pel,exl,labe,lin,lex,value,error,ib
	
        data stype/ ' ','in Fm**2','in barns', 'in mb', 'in micbn'/
! DATA:
!    type      = 0       angular distribution for fixed energy
!              = 1       excitation & angular cross sections 
!              = 2       excitation cross section for fixed angle
!              = 3       excitation total cross section
!                         (ic=0: ia=0 is fusion;
!                                ia=1 is reaction xs,
!                                ia=-1 is angle-integrated elastic
!                                ia=-2 is total cross section (elastic+reaction)
!              = 4       excitation phase shift for fixed partial wave
!              = 5 	 search factor for bound state
!              = 6       value and error for a search parameter
!              = 7       value and error for pole energy in Brune basis (not yet implemented)
!              = 8       value and error for total formal width of pole in Brune basis (not yet implemented)
!              = 9       value and error for ANC of bound state

!    idir      = 0       cross-section data are given in absolute units.
!              = 1       cross-section data are ratio to rutherford.
!              = 2       cross sections are given in absolute units but will
!                             converted to ratio to rutherford.
!              =-1       cross sections are given as spectroscopic factors but will
!                             converted to absolute

!    iscale    =-1       dimensionless (eg ratio to rutherford if idir>0)
!              = 0       absolute cross-section units are fermi-squared/sr.
!              = 1       absolute scale is barn/sr.    (MeV-barn for idir=-1)
!              = 2       absolute scale is mb/sr. (MeV-mb for idir=-1)
!              = 3       absolute scale is micro-b/sr. (MeV-microbarn for idir=-1)
	
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
!				Change stdout recl on some machines
!	call stdout(6,140)
	nul = -124578   ! 'undefined'
	fine = 10000.    ! this will be the penalty for FAILs, eg in EIGCC
	interactive = .false.
	number_calls = 0

	written = .false.; written(3) = .true.
1	write(6,1001) 
 1001   format('SFRESCO - FRES 3.4: Search Coupled Reaction Channels'
     x          /)
!
!####### Read in search input
	write(6,*) 'Please give search file name'
	read(5,*,end=999,err=1) search_file
	write(6,1002) search_file
 1002   format(' SFresco search file = ',a50/)
        open(303,file=search_file,status='old')

	read(303,*) input_file,output_file,nparameters,ndatasets
        if(nparameters>mvars) then
          write(6,*) 'ONLY ROOM FOR ',mvars,' SEARCH VARIABLES'
          stop
          endif
        if(ndatasets>mds) then
          write(6,*) 'ONLY ROOM FOR ',mds,' DATASETS'
          stop
          endif
	nvars=nparameters ; datasets=ndatasets
	write(6,1003) input_file,output_file,nparameters,ndatasets
 1003   format(' Fresco input file = ',a50/,
     x         '   and output file = ',a50/,
     X         ' Search on',i3,' variables for ',i3,' datasets,')
	call MNINIT(5,6,33)
	open(33,file='minuit-saved.dat')

!
!####### Read in specification of search parameters
       	if(nparameters>0) write(6,'(/''  Define SEARCH VARIABLES:''/)')
        i0 = ichar('0')
	undef = .false.
	rterms=.false.
	stepi = 0.01
	srch_datanorm(:) = 0
        do ip=1,nparameters
        name='Var'//char(mod(ip/10,10)+i0)//char(mod(ip,10)+i0)
        kind=0;valmin=0;valmax=0; kp=0;pline=0;col=0;
        potential=nul; step=stepi; width=0
        nafrac=0;afrac=nul; energy=nul;jtot=0;par=0;channel=1;term=1
        dataset=1;datanorm=1.0;    nopot=.false.
        read(303,nml=variable)

        srch_kind(ip) = kind; srch_name(ip) = name
        srch_minvalue(ip)=valmin; srch_maxvalue(ip)=valmax
	if(abs(step-stepi)<1e-10.and.step>abs(width) ! make default steps small for small widths
     x             .and.abs(width)>1e-20) step=step*abs(width)
        srch_step(ip)=step
	srch_error(ip)=0  ! initially

         if(kind==1) then !   potential parameter
 	  write(6,1010) ip,name,kp,pline,col
1010	  format('   Variable',i3,'=',a10,' is potential KP=',i2,
     X       ', line=',i2,' col=',i2)
          srch_value(ip) = potential; 
          srch_kp(ip) = kp; srch_pline(ip)=pline; srch_col(ip)=col

         else if(kind==2) then ! spectroscopic amplitude
          write(6,1012) ip,name,nafrac
1012	  format('   Variable',i3,'=',a10,' is Afrac #',i3)
          srch_value(ip) = afrac; srch_nafrac(ip) = nafrac; 

         else if(kind==3) then ! R-matrix energy
 	  write(6,1013) ip,name,term,jtot,psign(par+2),nopot
1013	  format('   Variable',i3,'=',a10,' is energy of R-matrix term'
     X       ,i3,' at J/pi =',f5.1,a1,' [NoPot=',L1,']')
          srch_value(ip) = energy; srch_rterm(ip) = term
	  srch_jtot(ip) = jtot; srch_par(ip) = par
	  srch_nopot(ip) = nopot
	  rterms=.true.

         else if(kind==4) then ! R-matrix width
 	  write(6,1014) ip,name,term,channel
1014	  format('   Variable',i3,'=',a10,' is width of R-matrix term',
     X       i3,' in channel',i3)
          srch_value(ip) = width; srch_rterm(ip) = term
	  srch_r_ch(ip) = channel

         else if(kind==5) then ! dataset normalisation
          write(6,1018) ip,name,dataset
1018	  format('   Variable',i3,'=',a10,
     X       ' is normalisation for dataset ',i3)
          srch_value(ip) = datanorm; srch_datanorm(ip) = dataset 
        endif

         if(abs(srch_value(ip)-nul)>.001) then
          write(6,1019) srch_value(ip),step,valmin,valmax
1019	  format('     value ',f8.4,10x,' step ',f7.4,
     X           ', min,max ',2f8.4/)
	 else
          write(6,1020)  step,valmin,valmax
1020	  format('     value from Fresco input',
     X           ', step ',f7.4,', min,max ',2f8.4/)
	  undef = .true.
	 endif
       enddo
	  peli=0; exli=0; labi=0; lini=0; lexi=0  ! we have not read in fresco input yet!
   
!
!####### Read in experimental data sets
       	ndof = 0
	ranks(:) = .false.; maxrank=-1
	num_energies=0
        do id=1,ndatasets
	  type=0; angle=0; jtot=-1;par=0;channel=1
          data_file="="; xmin=0;delta=-1;idir=0;iscale=-1;lab=.false.
          ic=1;ia=1;k=0;q=0; abserr=.false.;  points=-1; ib=0
          energy=elab(1)
	  pel=peli; exl=exli; labe=labi; lin=lini; lex=lexi
	  write(6,*) ' Read definition of data set ',id
          read(303,nml=data)
          if(idir.eq.1) iscale=-1
          errortype='relative'; if(abserr) errortype='absolute'
          cmlab='CM '; if(lab) cmlab='LAB'
	  dat_file(id) = data_file
          write(6,1025) id,data_file,type,stype(iscale+1),
     x      		errortype,cmlab,ic,ia,idir
1025      format('  Read DATA SET ',i2,' from file ',a50/
     x     '    of type ',i2,'  ',a8,' with ',a8,' errors,',
     x     1x,a3,' for state/excit',2i3,',  idir=',i2/)
	  ranks(k) = .true.
	  maxrank = max(maxrank,k)
	  if(type==0.and.energy<0) write(6,1029) k,q
	  if(type==0.and.energy>0) write(6,1030) k,q,energy
	  if(type==1) write(6,1031) k,q
	  if(type==2) write(6,1032) k,q,angle,cmlab
          if(type<=3.and.ib>0) write(6,1033) 'gamma decay to state',ib
          if(type==3.and.ic>=1) write(6,1033)
          if(type==3.and.ic==0) then
                if(ia==0) write(6,1033) 'fusion'
                if(ia==1) write(6,1033) 'reaction'
                if(ia==-1) write(6,1033) 'elastic'
                if(ia==-2) write(6,1033) 'total'
                if(ia>=2) write(6,1033) 'outgoing',ia
                endif
	  if(type==4) write(6,1034) jtot,psign(par+2),channel
	  if(type==5) write(6,1035) 
	  if(type==6) write(6,1036) par,value,error,abserr
	  if(type==9) write(6,1039) 
1029	  format('  Angular cross section T',2i1)
1030	  format('  Angular cross section T',2i1,
     x				' for energy',f8.3,' MeV')
1031	  format('  Excitation and Angular cross sections T',2i1)
1032	  format('  Excitation cross section T',2i1,
     x				' for angle',f8.3,' deg ',a3)
1033      format('  Excitation cross section',:,' for ',a,:,' in',i3)
1034	  format('  Excitation phase shift in channel',f5.1,a1,' #',i2)
1035	  format('  Known search parameters for bound states kn=x')
1036	  format('  Bound state search parameter',i3,' value and error',
     x        ' are',2f10.5,' (abserr=',l1,')')
1037      format('  R-matrix Brune energy of term ',i3,
     x        ' value and error are',/2f10.5,' (abserr=',l1,')', i4)
1038      format('  R-matrix Brune total formal width of','term ',i3,
     x        ' value and error are',/,2f10.5,' (abserr=',l1,')', i4)
1039      format('  Known ANC of bound state components kn=x')

	  if(type<0.or.type>9) write(6,*) 'Unrecognised data type ',type
	  data_type(id) = type

	if(type==6) then
	      datavals(1,id)=value
	      if(abserr) error = error*value
 	      dataerr(1,id)=error
 		datalen(id) = 1
 	      ip = 2
        else if(type==7 .or. type==8) then
              datavals(1,id)=value
              if(.not.abserr) error = error*value
              dataerr(1,id)=error
              data_term(id) = term
              rm_Brune(term) = .true.
              par = 0
              do ip=1,nparameters
               if(srch_kind(ip)==3 .and. srch_rterm(ip)==term) par=ip
              enddo
              data_par(id) = par

                datalen(id) = 1
              ip = 2
          if(type==7) write(6,1037) term,value,error,abserr,par
          if(type==8) write(6,1038) term,value,error,abserr,par
	else
              factor=1.0
	      if(iscale<0) factor=1. 
	      if(iscale==0) factor=10.
	      if(iscale>0) factor=1000.0/10.0**(3*(iscale-1)) 
	  inf=306
	  if(data_file=="=") inf=303
	  if(data_file=="<") inf=5
          if(inf==306) then
             write(6,*) ' To open file for dataset ',id,': ',data_file
             open(inf,file=data_file,status='old')
             endif
          xvals = delta<=0.; x = xmin
          if(points<0) points=99999
	  write(6,*) ' Read data set ',id
          do 10 ip=1,points
	    if(type==1) then
              read(inf,end=111,fmt=*,err=11) x,a,val,err
            else if(xvals.or.type==5.or.type==9) then
              read(inf,end=111,fmt=*,err=11) x,val,err
            else
              read(inf,end=111,fmt=*,err=11) val,err
              x = x + delta            
            endif
              if(x<0) go to 11
	      val =factor*val
	      if(abserr) then
	         err =factor*err
	        else
	         err = val * err
	        endif  ! Now all errors are absolute (and mb, except for r/ruth)
	      datavals(ip,id)=val
 	      dataerr(ip,id)=abs(err)

            if(type==0) then        ! angular distribution for fixed energy
              datangles(ip,id) = x
	      data_energies(ip,id) = energy
             else if(type==1) then  ! excitation and angular cross sections
 	      data_energies(ip,id) = x
              datangles(ip,id) = a
             else if(type==2) then  ! excitation cross section for fixed angle
 	      data_energies(ip,id) = x
	      datangles(ip,id) = angle
             else if(type==3) then  ! excitation total cross section 
 	      data_energies(ip,id) = x
             else if(type==4) then  ! excitation phase shift
 	      data_energies(ip,id) = x
             else if(type==5) then  ! bound state search factor
 	      bs_no(ip,id) = nint(x)
             else if(type==9) then  ! bound state ANC
 	      bs_no(ip,id) = nint(x)
	     endif
  10      continue
          ip = points+1
	  go to 111
  11      backspace inf
	endif
 111      datalen(id) = ip-1
          if(inf==306) then
!             write(6,*) ' Close file for dataset ',id,': ',data_file
             close(inf)
             endif
  	  ndof = ndof + datalen(id)
          data_idir(id)=idir; data_idir1(id)=idir; data_lab(id)=lab; 
          data_ic(id) = ic; data_ia(id) = ia;  data_ib(id) = ib
          data_rank_k(id) = k; data_rank_q(id) = q; 
	  data_ch(id) = channel; data_jtot(id)=jtot; data_par(id)=par
	  if(data_type(id)/=6) then
          if(energy<0) write(6,*) ' ',datalen(id),' data points:'
          if(energy>0) write(6,*) ' ',datalen(id),' data points',
     x              ' for lab energy ',real(energy)
     	  endif
	  if(datalen(id)==0) then
	    write(6,*) '   NO DATA POINTS !! Stop'
	    stop
	    endif
            if(pel.le.0) pel = 1
            if(exl.le.0) exl = 1
            if(labe.eq.0) labe = pel
            if(lin.eq.0) lin = 1
            if(lex.le.0) lex = 1
	    data_pel(id) = pel
	    data_exl(id) = exl
	    data_labe(id) = labe
	    data_lin(id) = lin
	    data_lex(id) = lex
            if(pel+exl+labe+lin+lex>5) then
	      write(6,1250) pel,exl,labe,lin,lex
 1250       format('     Incoming partition',I3,' in state #',I2,
     X             ',   Lab Energy for part.',I3,' Nucleus',I2,
     X             ' in Excitation pair',I2/)
              endif
	 if(type==0) then
          write(6,*) ' Angle '//cmlab//' Datum     Absolute error'
	  do ip=1,datalen(id)
	  write(6,12) datangles(ip,id),datavals(ip,id),dataerr(ip,id)
  12 	   format(1x,f8.3,2g12.4)
	  enddo
	 else if(type==1) then
          write(6,*) '   Energy Angle '//cmlab//
     x               ' Datum     Absolute error'
	  do ip=1,datalen(id)
	  write(6,13) data_energies(ip,id),datangles(ip,id),
     x			datavals(ip,id),dataerr(ip,id)
  13 	   format(1x,2f8.3,2g12.4)
	  enddo
	 else if(type==5.or.type==9) then
          write(6,*) '   Bound state   Target     Absolute error'
	  do ip=1,datalen(id)
	  write(6,131) bs_no(ip,id),datavals(ip,id),dataerr(ip,id)
  131 	   format(1x,i8,3x,2f12.4)
	  enddo
	 else if(type==6) then
          write(6,*) '   Search param  Target     Absolute error'
	  do ip=1,datalen(id)
	  write(6,132) data_par(id),srch_name(par),
     x                 datavals(ip,id),dataerr(ip,id)
  132 	   format(1x,i2,':',a8,2f12.4)
	  enddo
	 else
          write(6,*) '  Energy   Datum     Absolute error'
	  do ip=1,datalen(id)
	  write(6,12) data_energies(ip,id),datavals(ip,id),
     x			dataerr(ip,id)
	  enddo
	 endif
          close(1)
	  neng=1
	if(type<5) then
	  if(energy<0) neng=0
	  if(type>0) neng=datalen(id)
	  do ip=1,neng
	    if(type>0) energy=data_energies(ip,id)
	  ien=0
	  do ie=1,num_energies
	   if(abs(energy-energy_list(ie))<1e-5) then
 	     ien=ie; go to 14		! found existing energy
	   endif
	  enddo
	  ien = num_energies+1 ! list new energy
	  num_energies = ien
	  if(ien>maxen) then
	    write(6,*) 'Should increase PARAMATER maxen!'
	    stop
	    endif
	  energy_list(ien) = energy
   14     continue
   	  enddo  
	endif

        enddo     
        close(2)  
        write(6,*) 

!			Pre-read input to find SOME array-parameter limits
!
!####### Read in main fresco input
	ki = 306
	ko = 3
	koe = 307
!	call machine(mach)
	mach = 1
        open(ki,file=input_file,status='old')
	open(ko,form='formatted',delim='apostrophe')
	open(koe,file=trim(output_file)//'-init')
	write(6,*) ' FRESCO output to ',trim(output_file)//'*'
        write(6,*) 

	call freadf(ki,ko,koe,TMP,lmax,jbord,jump)	
	close(ki); close(koe)
        NANGL = (abs(THMAX) - THMIN)/abs(THINC) + 1 + 0.5
        allocate (theoryplot(max(mdl,NANGL),datasets))
	gettheoryplot = .false.
	  if(num_energies==0.and.elab(1)>0.) then
	    num_energies = 1
	    energy_list(1) = elab(1)
	    endif
	  if(num_energies>0) then
	     write(6,15) (energy_list(ien),ien=1,num_energies)
   15	  format(' Calculate scattering at energies'/(1x,10f8.3))
             else
	     write(6,*) ' No scattering energies'
	     endif
        write(6,*) 
	
	ki = 3
	ko = 308
	open(ko,file=output_file)
	koi= ko
	written(ko) = .true.
!	rewind ki
	noerror = .true.
	final = .false.
	final = .true.
	MAXCH = 0; MAXICH = 0  	 !	no arrays allocated
       	
C			DO IT!

	call fr

	final = .false.
	if(undef) then
	  do ip=1,nvars
	  do iof=ko,6,6-ko
	  write(iof,1048) ' ',ip,srch_name(ip),srch_value(ip)
	  enddo
	  enddo
	 endif
	 totchisq = sum(data_chisq(1:datasets))
	if(abs(totchisq/ndof)<1e5) then
	write(6,1040) totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	else
	write(6,10401) totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	endif
1040	format(/'  Total ChiSq/N =',f11.4,:' from ',10f10.3/
     x      (:/'#',26x,6f10.2))
10401	format(/'  Total ChiSq/N =',1p,e11.4,:' from ',10e10.3/
     x      (:/'#',26x,6e10.3))
1041	format(a1,' Variable ',i3,':',a10,'=',g12.4,
     x         ': ChiSq/N =',f10.3,:' from ',10f10.3,
     x                             (:/'#',26x,6f10.3))
10411	format(a1,' Variable ',i3,':',a10,'=',g12.4,
     x         ': ChiSq/N =',1p,e10.3,:' from ',10e10.3,
     x                                (:/'#',26x,6e10.3))
1042	format(/a1,' ChiSq/N =',f10.3,:' from ',6f10.3,
     x                              (:/'#',26x,6f10.3))
10421	format(/a1,' ChiSq/N =',1p,e10.3,:' from ',6e10.3,
     x                                 (:/'#',26x,6e10.3))

        do id=1,ndatasets   ! indicate data conversions at the end of the first run:
         idir = data_idir(id)
           if(idir==2) then  !  convert to ratio to Rutherford, the first time
              data_idir(id) = 1
              endif
           if(idir==-1.and.data_type(id)<=2) then  !  convert cross sections to absolute, the first time
              data_idir(id) = 0
              !write(191,*) ' set idir ',id,' to zero',IK
              endif
        enddo

	do ip=1,nvars
	if(srch_kind(ip)==2) then
	  write(6,1043) ip,srch_name(ip),srch_afrac_overlap(ip)
1043	  format(' Note: variable ',i3,'=',a10,' is overlap ',a50)
	endif
!			Give variable names to Minuit
	call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
	enddo
!
!####### INTERACTIVE SECTION

100	continue
	write(6,advance='no',fmt='(''sfresco> '')')

	read(5,'(a80)',end=999) cmd
	if(cmd(1:2)=='EX'.or.cmd(1:2)=='ex') then
	  stop

	else if(cmd(1:1)=='Q'.or.cmd(1:3)=='q') then
	  do ip=1,nvars
	  do iof=ko,6,6-ko
	   if(srch_error(ip)>1e-20.or.noerror) then
	    t = abs(srch_value(ip))
	    if(t<1e-3.or.t>1e3) then
	    write(iof,10481) ' ',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip)
	    else
	    write(iof,1048) ' ',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip)
	    endif
     	   else
	    write(iof,1049) ' ',ip,srch_name(ip),srch_value(ip)
	   endif
1048	  format(a1,'   Var ',i3,'=',a10,' value ',f12.6,:,
     x          ', step ',f8.4,', error ',f8.4)
10481	  format(a1,'   Var ',i3,'=',a10,' value ',1p,e12.4,:,
     x          ', step ',g8.1,', error' ,g9.2)
1049	  format(a1,'   Var ',i3,'=',a10,' value ',f12.6,' fixed')
	  enddo
	  enddo

	else if(cmd(1:3)=='SET'.or.cmd(1:3)=='set') then
	  read(cmd(4:100),*,err=980,end=981) ip,val
	  if(ip<1.or.ip>nvars) go to 982
	  do iof=ko,6,6-ko
	  write(iof,1050) ip,srch_name(ip),val
1050	  format(/' SET',i3,'=',a10,' to ',g12.4,:,' from',g12.4)
	  enddo
	  srch_value(ip) = val
	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
	 call fr
	 totchisq = sum(data_chisq(1:datasets))
         if(datasets>1) then
	  if(abs(totchisq/ndof)<1e6) then
 	   write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	   else
 	   write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	  endif
	 else
	  if(abs(totchisq/ndof)<1e6) then
 	    write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof
	    else
 	    write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof
	    endif
	 endif
	  if(abs(totchisq/ndof)<1e6) then
 	   write(ko,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	   else
 	   write(ko,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	   endif


	else if(cmd(1:3)=='REA'.or.cmd(1:3)=='rea') then
	  read(cmd(5:100),*,err=980,end=981) plot_file
	  write(6,*) ' Replace parameters by those in plot file ',
     x        '<'//trim(plot_file)//'>'
	  open(304,file=plot_file)
	  if(index(plot_file,'snap')>0) then
10377	    read(304,10382,end=10380) ic,chi
	    write(6,10379) ic,chi
10379	    format('  Reading snap at #',i5,' with Chisq =',1p,e12.3)
	    read(304,10381,end=10380) (srch_value(ip),ip=1,nvars)
	    go to 10377
10380	    do ip1=1,nvars
!      	    write(6,10485) ip1,name,srch_name(ip1),srch_value(ip1)
	      call MNPARM(ip1,srch_name(ip1),
     x          srch_value(ip1),srch_step(ip1),
     x		srch_minvalue(ip1),srch_maxvalue(ip1),ierflg)
	      if(ierflg>0) write(6,*) 'MNPARM for ',ip1,
     x        			      ': IERFLG =',ierflg
            enddo
10381	 	format(5e12.5)
10382	 	format(i6,1p,e12.4,1x,0p,6f10.5)
	  else
	  do ip=1,nvars
	  read(304,10482,err=102,end=102)  ip1,name,val
10482	  format(1x,7x,i3,1x,a10,7x,e12.4)
	  if(ip1<1.or.ip1>nvars) go to 982
!	  do iof=ko,6,6-ko
	   iof=ko
	  if(abs(srch_value(ip1)-val)>1d-12) 
     x 	    write(iof,10485) ip1,name,srch_name(ip1),val,srch_value(ip1)
10485	  format(' READ',i3,'=',2a10,' to ',g12.4,:,' from',g12.4)
!	  enddo
	  srch_value(ip1) = val
	  call MNPARM(ip1,srch_name(ip1),srch_value(ip1),srch_step(ip1),
     x		srch_minvalue(ip1),srch_maxvalue(ip1),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip1,': IERFLG =',ierflg
	  enddo
	  endif
	  go to 103
102	  write(6,*) ' Read plot file ended'
103	 call fr
	 totchisq = sum(data_chisq(1:datasets))
	 close(1)
         if(datasets>1) then
	  if(abs(totchisq/ndof)<1e6) then
 	   write(6,1042) ' ',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	   else
 	   write(6,10421) ' ',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	  endif
	 else
	  if(abs(totchisq/ndof)<1e6) then
 	    write(6,1042) ' ',totchisq/ndof
	    else
 	    write(6,10421) ' ',totchisq/ndof
	    endif
	 endif
	  if(abs(totchisq/ndof)<1e6) then
 	   write(ko,1042) ' ',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	   else
 	   write(ko,10421) ' ',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	   endif

!	else if(cmd(1:3)=='REL'.or.cmd(1:3)=='rel') then
!	  read(cmd(4:100),*,err=980,end=981) ip,step
!	  if(ip<1.or.ip>nvars) go to 982
!	  do iof=ko,6,6-ko
!	  write(iof,10501) ip,srch_name(ip),step
!10501	  format(/' REL',i3,'=',a10,' with step size =',f10.6)
!	  enddo
!	  srch_step(ip) = step
!	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
!     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
!	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg

	else if(cmd(1:3)=='FIX'.or.cmd(1:3)=='fix') then
	  read(cmd(4:100),*,err=980,end=981) ip
	  if(ip<1.or.ip>nvars) go to 982
	  do iof=ko,6,6-ko
	  write(iof,10502) ip,srch_name(ip)
10502	  format(/' FIX',i3,'=',a10)
	  enddo
	  srch_step(ip) = 0.0
	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg

	else if(cmd(1:4)=='STEP'.or.cmd(1:4)=='step') then
	  step = 0.01
	  read(cmd(5:100),*,err=980,end=981) ip,step
	  if(ip<1.or.ip>nvars) go to 982
	  do iof=ko,6,6-ko
	  write(iof,10503) ip,srch_name(ip),step
10503	  format(/' STEP ',i3,'=',a10,' to ',g12.4)
	  enddo
	  srch_step(ip) = step
	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg

       else if(cmd(1:5)=='ESCAN'.or.cmd(1:5)=='escan') then
	  read(cmd(6:100),*,err=980,end=981) emin,emax,de
	  ne = nint((emax-emin)/abs(de)) + 1
	  loge = de < 0.
	  if(ne>maxen) write(6,*) ' Only room for ',maxen,' points'
	  ne = min(ne,maxen)
	  ne = max(ne,2)
	  de = (emax-emin)/(ne-1)
	  write(6,1051)  emin,emax,ne,de
1051	  format(' EXCITATION function from ',f8.4,' to',f8.4,
     X		 ' in ',i5,' steps of',f8.5)
          is = num_energies; esave(:) = energy_list(:)
	  num_energies=ne
	  do i=1,ne
	  energy_list(i) = emin + (i-1)*de ! linear
           if(loge) then            ! log
             CF = (log(emax)-log(emin))/(emax-emin)
             CG =  log(emin) - CF * emin
             energy_list(i) = exp(CF*energy_list(i) + CG)
             endif

	  enddo
	  write(6,15) (energy_list(ien),ien=1,num_energies)
	  final = .true.
	  call fr
          num_energies=is; energy_list(:) = esave(:) ; final=.false.
	  write(6,*) ' File 71 has all phase shifts, and '
          write(6,*) ' file 40 the ',
     x      'fusion, reaction and non-elastic cross sections'
          write(6,*) ' file 35 the S-factors (CM energies)'
          write(6,*) ' file 75 the S-factors (LAB energies)'

       else if(cmd(1:4)=='SCAN'.or.cmd(1:4)=='scan') then
	  step=-1
	  read(cmd(5:100),*,err=980,end=981) ip,val1,val2,step
	  if(ip<1.or.ip>nvars) go to 982
	  if(step<=0) step=srch_step(ip)
	  write(6,1052) ip,srch_name(ip),val1,val2,step
1052	  format(/' SCAn',i3,'=',a10,' from ',g12.4,
     X       ' to ',g12.4,' in steps of ',g12.4)
	  oldval = srch_value(ip)
          ns = (val2-val1)/step+1
          do is=1,ns
	  val = val1 + (is-1)*step
	  srch_value(ip) = val
  	  call fr
	  totchisq = sum(data_chisq(1:datasets))
           if(datasets>1) then
	    if(abs(totchisq/ndof)<1e6) then
 	     write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	     else
 	     write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	     endif
	   else
	    if(abs(totchisq/ndof)<1e6) then
 	     write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof
	     else
 	     write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof
	     endif
	   endif
	    if(abs(totchisq/ndof)<1e6) then
  	     write(ko,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	     else
  	     write(ko,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	     endif
          enddo
	 srch_value(ip) = oldval
	 call fr
	 totchisq = sum(data_chisq(1:datasets))

        else  if(cmd(1:4)=='SHOW'.or.cmd(1:4)=='show') then
	 call fr
	 totchisq = sum(data_chisq(1:datasets))
!	 write(6,'('' Chisq='',10f10.3)') data_chisq(1:datasets)
	  do iof=ko,6,6-ko
	  do id=1,datasets
	  if(data_type(id)/=6) then 
	  datanorm=1.0
!		Adjust any datanorm search parameter!
	   do ip=1,nvars
	   if(srch_kind(ip)==5.and.srch_datanorm(ip)==id) 
     x	       datanorm = datanorm * srch_value(ip)
	   enddo

          write(iof,*) 
          write(iof,*) ' Dataset ',id,' <'//trim(dat_file(id))//'>'
	  if(data_type(id)==0) then
          write(iof,*) '   Angle   Datum      Abs. error  Theory',
     x               '         Chi'
	  else if(data_type(id)==1) then
          write(iof,*) '   Energy  Angle      Datum   Abs. error  ',
     x               '      Theory         Chi'
	  else if(data_type(id)==5.or.data_type(id)==9) then
          write(iof,*) '   Bound state    Datum      Abs. error  ',
     x               'Theory         Chi'
     	  else
          write(iof,*) '  Energy   Datum      Abs. error  Theory',
     x               '         Chi'
          endif
	  do ip=1,datalen(id)
            chi=(theoryvals(ip,id)-datanorm*datavals(ip,id))/
     x			(datanorm*dataerr(ip,id))
	  if(data_type(id)==0) then
	  write(iof,1062) datangles(ip,id),datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
     	  else if(data_type(id)==1) then
	  write(iof,1062) data_energies(ip,id),datangles(ip,id),
     x          datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
	  else if(data_type(id)==5.or.data_type(id)==9) then
	  write(iof,10621) bs_no(ip,id),datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
     	  else
	  write(iof,1062) data_energies(ip,id),datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
	  endif
1062 	   format(1x,f8.3,3g12.5,8f12.4)
10621 	   format(1x,i8,7x,3g12.5,f10.4)
	  enddo
	  endif ! type/=6
	  enddo
	if(abs(totchisq/ndof)<1e6) then
	write(iof,1040) totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	else
	write(iof,10401) totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	endif
	  enddo	

        else  if(cmd(1:4)=='PLOT'.or.cmd(1:4)=='plot'
     x       .or.cmd(1:4)=='LINE'.or.cmd(1:4)=='line') then
     	  nodat = cmd(1:4)=='LINE'.or.cmd(1:4)=='line'
	  tag = ' '
	  read(cmd,'(5x,a50)',err=980,end=981) tag
!	  write(6,*) ' tag = <'//trim(tag)//'>'
	  do k=0,min(2,maxrank)
	  plot_file = 'search.plot'
	  if(tag(1:1).ne.' ') plot_file = trim(tag)
	  if(k>0) plot_file = trim(plot_file)//char(k+ichar('0'))
	  open(304,file=plot_file)
	  do ip=1,nvars  ! remind parameters for fit
	   if(srch_error(ip)>1e-20.or.noerror) then
	    t = abs(srch_value(ip))
	    if(t<1e-3.or.t>1e3) then
	    write(304,10481) '#',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip)
	    else
	    write(304,1048) '#',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip)
	    endif
     	   else
	    write(304,1049) '#',ip,srch_name(ip),srch_value(ip)
	   endif
	  enddo
		if(k==0) then
		gettheoryplot = .true.; final = .true.
	  call fr
		gettheoryplot = .false.; final = .false.
          totchisq = sum(data_chisq(1:datasets))
		endif
	if(totchisq/ndof+
     x      maxval(data_chisq(1:datasets)/datalen(1:datasets))<1e6) then
 	  write(304,1042) '#',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	  else
 	  write(304,10421) '#',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)
	  endif
 	  write(6,1042) ' ',totchisq/ndof,
     x              (data_chisq(id)/datalen(id),id=1,datasets)

	    write(304,1063) 'Search file: '//trim(search_file),
     x 	                  'Fresco input: '//trim(input_file)
1063	    format('@subtitle "',a,'; ',a,'"',/,
     x             '@subtitle size 0.7'/'@legend ON')
	  is = -1; idp = -1
	  do id=1,datasets
	    if((data_rank_k(id)==k.or.(k==0.and.data_type(id)==5))
     x          .and.data_type(id)/=6) then
	    is=is+1; idp=idp+1
	   if(.not.nodat) then
            cmlab='CM '; if(data_lab(id)) cmlab='LAB'
	    if(is<=9) then
	    write(304,10631) is,cmlab,id,trim(dat_file(id)),is,is,is,
     x                     mod(idp+2,10),is,idp+1,is,is
10631	    format('@TYPE xydy'/
     x             '@legend string ',i2,' "',a3,' set ',i2,':',a,'"'/
     x             '@ s',i1,' linestyle 0',/
     x             '@ s',i1,' errorbar length 0.28'/
     x             '@ s',i1,' symbol ',i2/
     x             '@ s',i1,' color ',i2/
     x             '@ s',i1,' symbol fill 1'/
     x             '@ s',i1,' symbol size 0.5')
	    else
	    write(304,10632) is,cmlab,id,trim(dat_file(id)),is,is,is,
     x                     mod(idp+2,10),is,idp+1,is,is
10632	    format('@TYPE xydy'/
!    x             '@subtitle "',a,'; ',a,'"',/,
     x             '@subtitle size 0.7'/'@legend ON',/,
     x             '@legend string ',i2,' "',a3,' set ',i2,':',a,'"'/
     x             '@ s',i2,' linestyle 0',/
     x             '@ s',i2,' errorbar length 0.28'/
     x             '@ s',i2,' symbol ',i2/
     x             '@ s',i2,' color ',i2/
     x             '@ s',i2,' symbol fill 1'/
     x             '@ s',i2,' symbol size 0.5')
	    endif
	  datanorm=1.0
!		Adjust any datanorm search parameter!
	   do ip=1,nvars
	   if(srch_kind(ip)==5.and.srch_datanorm(ip)==id) 
     x	       datanorm = datanorm * srch_value(ip)
	   enddo

	   elast = data_energies(1,id)
	   do ip=1,datalen(id)
	    e = datanorm*dataerr(ip,id)
	    d = datanorm*datavals(ip,id)
	    if(e>abs(d)) e = .9*abs(d)
	      if(data_type(id)==1) then
              if(abs(elast-data_energies(ip,id))>1e-5) then
		write(304,*) 'NEW energy'
                elast = data_energies(ip,id)
 	      endif
 	      endif
	 
	    if(data_type(id)==0) write(304,1062) datangles(ip,id),d,e
	    if(data_type(id)==1) write(304,1062) datangles(ip,id),d,e
	    if(data_type(id)==5) write(304,1062) bs_no(ip,id)+0.,d,e
	    if(data_type(id)==9) write(304,1062) bs_no(ip,id)+0.,d,e
	    if(data_type(id)>1.and.data_type(id)<5) then
	      sfa = 1.
	      if(data_idir1(id)==-1)  then
		  ecmi = data_energies(ip,id) * ecmrat
		  etai = etarat / SQRT(ecmi)
		  sfa = exp(2d0*PI*etai) *  ecmi
		   endif
    		 write(304,1062) data_energies(ip,id),d*sfa,e*sfa
	      endif
	   enddo
	  write(304,*) '&'
	  is=is+1
	  endif
            cmlab='  '; 
	    if(data_type(id)<=1) then  ! theoryplot only available in cm..
	       cmlab='CM'
	       if(data_lab(id)) then
	        write(6,*)  '  Note: in search.plot, data is LAB,',
     x             ' continuous theory plot is CM, so not given!'
	       endif
	    endif

	    if(is<10) write(304,1064)   cmlab,id,is,idp+1
	    if(is>9)  write(304,10641)   cmlab,id,is,idp+1
1064	    format('@TYPE xy',/,'#     ',a2,' Dataset ',i3/
     x             '@ s',i1,' color ',i2)
10641	    format('@TYPE xy',/,'#     ',a2,' Dataset ',i3/
     x             '@ s',i2,' color ',i2)
	  if(data_type(id)==0) then
	   if(.not.data_lab(id)) then
	     do ia=1,NANGL
	      write(304,1062) thmin+(ia-1)*abs(THINC),theoryplot(ia,id)
	     enddo
	   else
	     do ia=1,datalen(id)
	      write(304,1062) datangles(ia,id),theoryvals(ia,id)
	     enddo
	   endif
	  else if(data_type(id)==1) then
	   do ia=1,NANGL
!	    write(304,1062) thmin+(ia-1)*abs(THINC),theoryplot(ia,id)
	   enddo
	  else if(data_type(id)==5.or.data_type(id)==9) then
	   do ia=1,datalen(id)
	    write(304,10621) bs_no(ia,id),theoryvals(ia,id)
	   enddo
	  else if(data_type(id)==6) then
c		 no nothing here
	  else  
	   do ia=1,datalen(id)
	      sfa = 1.
	      if(data_idir1(id)==-1)  then
		  ecmi = data_energies(ia,id) * ecmrat
		  etai = etarat / SQRT(ecmi)

!		  ecmi = data_energies(ia,id) *
!     x               MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
!		  etai = ETACNS * MASS(2+1,pel) * MASS(2+2,pel)
!     X                          * SQRT(RMASS(pel)/ ecmi)

		  sfa = exp(2d0*PI*etai) *  ecmi
	          endif
	    write(304,1062) data_energies(ia,id),theoryvals(ia,id)*sfa
	   enddo
	  endif
	   write(304,*) '&'
	   endif  ! correct k
	  enddo	
	  close(304)
	  write(6,*) ' xmgr file written: ',plot_file
	  enddo  ! rank k
        else  if(cmd(1:3)=='MIN'.or.cmd(1:3)=='min') then  ! call minuit
	  interactive = .true.
!		call fcn(1,srch_error,fval,srch_value,1,futil)
!		write(6,*) 'fcn =',fval
	   inquire(105,opened=op)
	   if(.not.op) open(105,file=trim(output_file)//'-trace')
	   inquire(106,opened=op)
	   if(.not.op) open(106,file=trim(output_file)//'-snap')
	  plist(1) = 1./ndof
	  write(6,1070) plist(1)
1070	  format(/'   Call MINUIT',/
     x            '   with NOGradient, STRategy=0, ERRordef=',f8.4)
	  call MNCOMD(fcn,'set nogradient',ICONDN,futil)
	  call MNCOMD(fcn,'set strat 0',ICONDN,futil)
	  call MNEXCM(fcn,'set errordef',plist(1),1,ICONDN,futil)
	  call MNINTR(fcn,futil)
	  do ip=1,nvars
!			Get variable values back from Minuit
	   call MNPOUT(ip,srch_name(ip),srch_value(ip),srch_error(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ivarbl)
	  enddo
	  noerror = .false.
        else
	  go to 980
	endif
	go to 100
980	  write(6,*) 'Unrecognised command : ',cmd
	go to 100
981	  write(6,*) 'Incomplete input command : ', cmd
	go to 100
982	  write(6,*) 'Parameter number',ip,' outside defined range 1:'
     x                   ,nvars
	go to 100
983	  write(6,*) 'Set number',id,' outside defined range'
	go to 100
999	stop
	end
	logical function intrac()
	use searchpar, only: interactive
	intrac = interactive
	return
	end
	subroutine futil
	end
	subroutine fcn(npar,grad,fval,xval,iflag,futil)
	use searchpar
	use searchdata
	implicit real*8(a-h,o-z)
	real*8 fval,grad(*),xval(*)
	external futil
	if(iflag==1) then
	  number_calls = 0
!					Initialise
	else if(iflag==2) then
!					Gradients -> grad
	endif
	  number_calls = number_calls+1
	  do ip=1,nvars
	   srch_value(ip) = xval(ip)
	  enddo
	  penalty = 0.

	  call fr

	  totchisq = sum(data_chisq(1:datasets)) + penalty
	  fval = totchisq/ndof
	  ip = min(6,datasets)
	  if(datasets==1) ip=0
	  write(105,2) number_calls,fval, 
     x           (data_chisq(id)/datalen(id),id=1,ip)
1	 format(1p,5e12.5)
2	 format(i6,1p,e12.4,1x,0p,6f10.5)
3	 format(i6,1p,e12.4,1x,1p,10e10.2)
	  ip = min(10,datasets)
	  write(106,3) number_calls,fval,
     x           (data_chisq(id)/datalen(id),id=1,ip)
	  write(106,1) xval(1:nvars)
	  call flush(105)
	  call flush(106)
	
	if(iflag==3) then
!					Wrapup
	endif
	return
	end
	

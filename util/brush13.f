*---------------------------------------------------------------------
      program interactive_data_set_create
*---------------------------------------------------------------------
*     JAT Jan 2005 [Version zzzz] Incorporates JLM potential as option
*     and can read BAB HF matter densities or other external file.
*     MSU/Betty Tsang July 2004 visit - for above. 
*---------------------------------------------------------------------
*     Revised March 2005 to include n and p densities in determining
*     alpha (isovector) rather than (N-Z)/A.
*
*     Modified: FSU (Kemper/Roeder) to include (d,n) and (n,d) ('05)
*     Modified: MSU (Tsang) to include (3He,d) (June '05)
*
*     front5z.f is the front end - as per normal - but, if Vso=/0,
*     then this version requests input also for the bound state 
*     potential's spin-orbit geometry (at the end). You thus
*     can produce data sets with a specified rso and and aso
*     that are different to the central r0 and a0 - as previous. 
*
*     Numbering of reactions has reverted to earlier conventions
*     for (p,d) and (d,p) to allow earlier data set usage - MSU 2005
*---------------------------------------------------------------------
*     Version 6z - allows negative ireac values, for rotation
*     of amplitudes and then m-dependent cross sections - WNC 2006
*---------------------------------------------------------------------
*     Version 7z - uses the Bauge JLM parameterisation - BT 2006
*---------------------------------------------------------------------
*     Version 8z - has option to read the Sao Paulo potential for 3He
*     Helio Dias and Betty Tsang - got postponed until early 2008
*     Requires no changes to twofnr itself so can use twofnr7.f
*     Also a number of additional input error traps - courtesy of
*     United Flight 929 - 15 March 2008
*---------------------------------------------------------------------
*     Version 9z: Having found an old code for the Watanabe and finite
*     range adiabatic potentials from the Reid soft core potential
*     and deuteron - these options have been introduced -  March 2008.
*---------------------------------------------------------------------
*     Version 10: Includes the triton global potential of Li et al. 
*     For Jeff Thomas (Surrey) April/May 2009 and also corrections to
*     the Bauge JLM parameterisation (due to exchange with Pang and
*     he with Bauge, January 2011 - e-mail records these changes)
*---------------------------------------------------------------------
*     Version 11: James Benstead (Surrey/AWE)/JAT. Version includes 
*     (p,t) and (t,p) reactions as di-neutron transfer - May 2011.
*---------------------------------------------------------------------
*     Version 12: Global A=3 potential GDP08 -- Pang et al. Jan 2012
*---------------------------------------------------------------------
*     Version 'Brush1': Ian Thompson (LLNL) - September 2011.
*     Also output Fresco input file as  *.fres
*---------------------------------------------------------------------
*     To keep step with front version being used - renamed brush12
*     JAT 2012. This version has corrected spin orbit weight for
*     read real spin-orbit formactors.
*---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character fname*12,titf*46,title*58,cmd*80
      character*12 form,potin,potout,potkind
      logical write15,pnlocs
      real*8 a(8),b(4)
      character*5 nuclei(4)
      character*2 NUCNAME(110)
      DATA NUCNAME/'NN','H ','He',
     +'Li','Be','B ','C ','N ','O ','F ','Ne',
     +'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     +'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co',
     +'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     +'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh',
     +'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',
     +'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     +'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +'Lu','Hf','Ta','W ','Re','Os','Ir','Pt',
     +'Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     +'Fr','Ra','Ac','Th','Pa','U ','Np','Pu',
     +'Am','Cm','Bk','Cf','Es','Fm','Md','No',
     +'Lr','Rf','Db','Sg','Bh','Hs','Mt'/

*---------------------------------------------------------------------
*     common blocks with channel potential parameters
*---------------------------------------------------------------------
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/cjlm/potr(900),poti(900),ijlm,ijlmp
      common/che3/pspr(900),pspi(900),psps(900),ihesp
      common/djlm/rlr,rli
      common/deut/iadia,iwat
      real*8 vreal(900),vimag(900),vspin(900)
!      common/fromfold/vreal(900),vimag(900),vspin(900),potin,potout
      common/fromfold/vreal,vimag,vspin
*---------------------------------------------------------------------
      print*,'====================================================== '
      print*,' BRUSH     (adapted from FRONT11)                      '
      print*,' Front end for FRESCO & TWOFNR zr-transfer data sets   '
      print*,'====================================================== '
      print*,' Data set identifier (xx in xx.fres: max 12 chars)    '
      read '(a)',fname
      open(16,file=trim(fname)//'.fres',status='unknown')
      open(17,file=trim(fname)//'.tran',status='unknown')
      open(18,file=trim(fname)//'.brin',status='new')
      open(15,file=trim(fname)//'.pots',status='unknown')
      write(18,'(a)') fname
      print '(a,a)','  >>>> ',fname
      print '(a,a)','  Output for twofnr, file: ',trim(fname)//'.tran'
      print '(a,a)','  Output for fresco, file: ',trim(fname)//'.fres'
      print '(a,a)','  Potentials > fresco: ',trim(fname)//'.pots'
      print '(a,a)','  Input saved to file: ',trim(fname)//'.brin'
      print*,'-------------------------------------------------   '
*---------------------------------------------------------------------
      nil=0
      unit=1.d0
*---------------------------------------------------------------------
*     version date
      iday=9; imon=09; iyear=2011
*---------------------------------------------------------------------
      iadia=0
      iwat=0
      ihesp=0
      write15 = .false.; pnlocs=.false.
      form='(5e14.7)'
*---------------------------------------------------------------------
*     the following used for ranges with jlm/bauge potentials
*---------------------------------------------------------------------
      ijlm=0
      rlr=0.d0
      rli=0.d0
*---------------------------------------------------------------------
*     title line
*---------------------------------------------------------------------
      print*,' Enter title information (for info only < 46 chars)'
      read '(a)',titf
      write(18,'(a)') titf
      print '(a,a)','  >>>> ',titf
      title=titf//fname
*---------------------------------------------------------------------
  932 print*,' Reaction type: [1] (p,d)   '
      print*,'                [2] (d,p)   '
      print*,'                [3] (n,d)   '
      print*,'                [4] (d,n)   '
      print*,'                [5] (d,t)   '
      print*,'                [6] (d,3He) '
      print*,'                [7] (3He,d) '
      print*,'                [8] (p,t)   '
      print*,'                [9] (t,p)   '
      read*,ireac
      ireaco=abs(ireac)
      if(ireaco.lt.1.or.ireaco.gt.9) goto 932
      write(18,*) ireac
      print*,' >>>> ',ireac
      ireaco=0
*---------------------------------------------------------------------
      if(ireac.lt.0) then
       ireaco=8
       ireac=abs(ireac)
       print*,' You have requested cross sections in a rotated '
       print*,' coordinate system. You must therefore specify '
       print*,' Euler angles (alfa,beta,gama) in units of pi '
       read*,angal,angbe,angga
       print*,' >>>> ',real(angal),real(angbe),real(angga)
       write(18,*) real(angal),real(angbe),real(angga)
       print '(a,a)','  Angles written to file: ','rot.'//fname
       open(19,file=trim(fname)//'.rot',status='unknown')
       write(19,*) angal,angbe,angga

	write(6,*) ' COORDINATE ROTATION NOT IMPLEMENTED FOR FRESCO'
	           stop
      endif
*---------------------------------------------------------------------
      write(16,1011) ireac,title,iday,imon,iyear
      write(17,101) ireaco,nil,iday,imon,iyear,title
  101 format(i1,i11,i2,i2,i4,2x,a)
 1011 format(i1,':',a,' with V',i2,'/',i2,'/',i5/'NAMELIST')
  102 format(0p8f10.4)
*---------------------------------------------------------------------
*     reaction masses
      m1=1
      if(ireac.gt.1) m1=2
      if(ireac.eq.3.or.ireac.eq.8) m1=1
      if(ireac.eq.7.or.ireac.eq.9) m1=3
      print*,' Laboratory incident energy per nucleon (MeV) '
      read*,energy
      write(18,*) real(energy)
      print*,' >>>> ',real(energy)
      energy=m1*energy
      print*, 'Total projectile energy = '
      print*, real(energy)
*---------------------------------------------------------------------
*     separation/binding energies for A=2 and A=3 systems (MeV)
*---------------------------------------------------------------------
      edeut=2.224573d0
      etrit=8.481821d0-edeut
      eheli=7.718058d0-edeut
*     binding energy of the triton
      etwo=8.481821d0
*---------------------------------------------------------------------
  933 print*,' Target mass (a1) and charge (z1)    '
      read*,a1,z1
      if(z1.gt.a1.or.a1.le.0.d0) go to 933
      write(18,*) real(a1),real(z1)
      print*,' >>>> ',real(a1),real(z1)
*---------------------------------------------------------------------
*     spins and charges - m1 has already been established
      s1=1.0
      zp1=1
      s2=0.5
      z2=z1
      zp2=1
      if(ireac.eq.1) then
       s1=0.5 
       a2=a1-1.0
       m2=2
       s2=1.0
      else if(ireac.eq.3) then
       zp1=0
       s1=0.5 
       a2=a1-1.0
       m2=2
       s2=1.0
       z2=z1-1.0
      else if(ireac.eq.2) then
       a2=a1+1.0
       m2=1
      else if(ireac.eq.4) then
       a2=a1+1.0
       m2=1
       zp2=0.0
       z2=z1+1.0
      else if(ireac.eq.5) then
       a2=a1-1.0
       m2=3
      else if(ireac.eq.6) then
       a2=a1-1.0
       m2=3
       z2=z1-1.0
       zp2=2.0
      else if(ireac.eq.7) then
       zp1=2
       s1=0.5 
       a2=a1+1.0
       m2=2
       s2=1.0
       z2=z1+1.0
       else if(ireac.eq.8) then
       s1=0.5
       a2=a1-2.0
       m2=3
       else if(ireac.eq.9) then
       s1=0.5
       a2=a1+2.0
       m2=1
      endif
*---------------------------------------------------------------------
*     energy and integration ranges line
*---------------------------------------------------------------------
      a(1)=unit
      a(2)=nil
      a(4)=nil
      a(6)=energy
      a(7)=nil
      a(8)=nil
      print*,'------------------------------------------------- '
 873  print*,' Integration ranges: [1] use defaults   '
      print*,' (defaults: 0-30 fm in 0.10 fm steps)   '
      print*,'                     [2] specify values '
      read*,iiran
      if(iiran.lt.1.or.iiran.gt.2) goto 873
      write(18,*) iiran
      print*,' >>>> ',iiran
      if(iiran.eq.1) then
       rmax=30.d0
       step1=0.10d0
      else
       print*,' max integration radius and step '
       read*,rmax,step1
       write(18,*) real(rmax),real(step1)
       print*,' >>>> ',real(rmax),real(step1)
      endif
      nrmax=nint(rmax/step1)
      step2=step1*(a1/a2)
      a(3)=rmax
      a(5)=nrmax
      print*,' integrations from 0 to',real(rmax),' fm'
      print*,' in steps of',real(step1),' fm'
      print*,' step in outgoing channel =',real(step2)
      print*  
      write(17,102) (a(i),i=1,7) 
      nr3max=nrmax+3
*---------------------------------------------------------------------
 874  print*,' number of partial waves [1] default (70) '
      print*,'                         [2] specify      '
      read*,ipw
      if(ipw.lt.1.or.ipw.gt.2) goto 874
      write(18,*) ipw
      print*,' >>>> ',ipw
      npw=70
      if(ipw.eq.2) then
       print*,' input number of partial waves (<90) '
       read*,npw
       write(18,*) npw
       print*,' >>>> ',npw
      endif
*---------------------------------------------------------------------
      print*,' Input the required centre of mass angles info:   '
      print*,' number of angles: step (degrees): starting value '
      print*,' (entering 0 0 0 will use 181  1.0  0.0 )          '
      read*,b(2),b(3),b(4)
      write(18,*) real(b(2)),real(b(3)),real(b(4))
      if (b(2).lt.0.01) then
       b(2)=181
       b(3)=1.0
       b(4)=0.0
      endif
      print*,' >>>> ',real(b(2)),real(b(3)),real(b(4))
*------------- FRESCO output so far
      rintp=0.3
      absend=-1  ! all partial waves up to npw!
      absend=1e-4  ! all partial waves until reaction < 1e-4 mb
	thmax =b(4)+b(2)*b(3)
	iter = 1   !  first-order DWBA !
	inh = 2    !  longitudinal recoil
      write(16,83) step1,rmax,rintp,npw,absend,inh,b(4),thmax,b(3),
     x       	 iter,energy
83    format(' &FRESCO   hcm=',f8.4,'  rmatch=',f8.3,' rintp=',f5.2,/
     x       '   jtmin=0  jtmax=',i5,' absend=',f8.4,' inh=',i1,/
     x       '   thmin=',f8.3,' thmax=',f8.3,' thinc=',f8.3,/
     x       '   iter=',i1,'   chans=1 listcc=0 smats=2 xstabl=1 ',/
     x       '   elab(1) =',f9.3,' /')
*---------------------------------------------------------------------
*     transferred angular momenta line
*---------------------------------------------------------------------
      a(1)=2.2
      print*,'------------------------------------------------- '
 911  if(ireac.eq.8.or.ireac.eq.9) then
       a(2)=0.0
       print*,' Enter quantum numbers L and J of transferred cluster'
       print*,' Uses simple di-neutron model, so S = 0 and L = J '
      else
       a(2)=0.5
       print*,' sp quantum numbers L and J of transferred nucleon'
      endif
      read*,ltr,rjtr
	str = a(2)  ! for FRESCO output
      if(abs(abs(ltr-rjtr)-a(2)).gt.0.1) then
       print*,' need |J-S| < L < J+S   '
       go to 911
      endif
      write(18,*) ltr,real(rjtr)
      print*,' >>>> ',ltr,real(rjtr)
 437  if(ireac.eq.8.or.ireac.eq.9)then
       print*,' number of nodes in di-neutron radial wave function '
      else
       print*,' number of nodes in nucleon sp radial wave function '
      endif
      print*,' (convention here: the lowest state has zero nodes)' 
      read*,nodes
      if(nodes.lt.0) go to 437
      write(18,*) nodes
      print*,' >>>> ',nodes
      a(3)=ltr
      a(4)=rjtr
      write(17,102) (a(i),i=1,4) 
      print*
*---------------------------------------------------------------------
*     sort out separation energy and/or Q-value
*---------------------------------------------------------------------
  935 if(ireac.eq.1.or.ireac.eq.2.or.ireac.eq.5) then
       print*,' specify : [1] neutron separation energy (>0 MeV)'
      else if(ireac.eq.8.or.ireac.eq.9) then
       print*,' specify : [1] two neutron separation energy (>0 MeV)'
      else
       print*,' specify : [1] proton separation energy (>0 MeV)'
      endif
      print*,' or        [2] reaction Q-value (MeV)            '
      read*,ietyp
      if(ietyp.lt.1.or.ietyp.gt.2) go to 935
      write(18,*) ietyp
      print*,' >>>> ',ietyp
      if(ietyp.eq.1) then
       print*,' transferred particle separation energy (MeV: >0) '
       read*,sn
       write(18,*) real(sn)
       print*,' >>>> ',real(sn)
*---------------------------------------------------------------------
*     compute Q-value
*---------------------------------------------------------------------
       if(ireac.eq.1.or.ireac.eq.3) then
        qval=-sn+edeut
       else if(ireac.eq.2.or.ireac.eq.4) then
        qval=sn-edeut
       else if(ireac.eq.5) then
        qval=-sn+etrit
       else if(ireac.eq.6) then
        qval=-sn+eheli
       else if(ireac.eq.7) then
        qval=sn-eheli
       else if(ireac.eq.8) then
        qval=-sn+etwo
       else if(ireac.eq.9) then
        qval=sn-etwo
       endif
       print*,'  Q-value is',real(qval),' MeV'
      else if(ietyp.eq.2) then
       print*,' reaction Q-value (MeV) '
       read*,qval
       write(18,*) real(qval)
       print*,' >>>> ',real(qval)
*---------------------------------------------------------------------
*     compute separation energy 
*---------------------------------------------------------------------
       if(ireac.eq.1.or.ireac.eq.3) then 
        sn=edeut-qval
       else if(ireac.eq.2.or.ireac.eq.4) then
        sn=edeut+qval
       else if(ireac.eq.5) then
        sn=etrit-qval
       else if(ireac.eq.6) then
        sn=eheli-qval
       else if(ireac.eq.7) then
        sn=eheli+qval
       else if(ireac.eq.8) then
        sn=etwo-qval
       else if(ireac.eq.9) then
        sn=etwo+qval
       endif       
       print*,'  Separation energy is',real(sn),' MeV'
       if(sn.lt.0.d0) then
        print*,' state is particle unbound - so cannot proceed '
        goto 935
       endif 
      endif
*---------------------------------------------------------------------
*     can now compute lab energy for final state potential (energy2)
*     print wavenumbers and look at likely mismatch of reaction
*---------------------------------------------------------------------
      print*,'================================================= '
      ecm1=energy*a1/(a1+m1)
      print*,' entrance channel cm energy ',real(ecm1)
      fmu1=a1*m1/(a1+m1)
      fkay1=0.2195376d0*sqrt(fmu1*ecm1) 
      ecm2=ecm1+qval
      print*,' exit     channel cm energy ',real(ecm2)
      if(ecm2.lt.0.d0) then
       print*,' reaction is below threshold'
       stop
      else if(ecm2.lt.10.d0) then
       print*,' reaction is near threshold'
      endif
      fmu2=a2*m2/(a2+m2)
      fkay2=0.2195376d0*sqrt(fmu2*ecm2) 
      rad=1.2d0*(a1**0.3333333333d0)
      rl1=fkay1*rad
      rl2=fkay2*rad
      print*,' wavenumbers and grazing angular momenta '
      print*,' kin  = ',real(fkay1),'  L(in ) = ',real(rl1)
      print*,' kout = ',real(fkay2),'  L(out) = ',real(rl2)
      rlmis=abs(rl1-rl2)
      print*,' so L mismatch of ',real(rlmis),' hbar'
      print*,' for an estimated radius of ',real(rad),' fm'
      print*,' and an L transfer of ',ltr,' hbar'
      energy2=ecm2*(a2+m2)/a2
*---------------------------------------------------------------------
*     entrance channel partial waves/nonlocality line
*---------------------------------------------------------------------
      print*,'================================================= '
      if(ireac.eq.1.or.ireac.eq.8) then
       print*,' incident (proton) channel information           '
      else if(ireac.eq.3) then
       print*,' incident (neutron) channel information          '
      else if(ireac.eq.7) then
       print*,' incident (3He) channel information              '
      else if(ireac.eq.9) then
       print*,' incident (triton) channel information           '
      else
       print*,' incident (deuteron) channel information         '
      endif
 876  print*,' nonlocality in incident channel [1] no  '
      print*,'                                 [2] yes '      
      read*,inonloc
      if(inonloc.lt.1.or.inonloc.gt.2) go to 876
      write(18,*) inonloc
      print*,' >>>> ',inonloc
      a(5)=0.d0 
      if(inonloc.eq.2) then
       if(ireac.eq.1.or.ireac.eq.8) then
        print*,' input proton nonlocality range (~0.85 fm) '
       else if(ireac.eq.3) then
        print*,' input neutron nonlocality range (~0.85 fm) '
       else if(ireac.eq.7) then
        print*,' input 3He nonlocality range (~0.20 fm) '
       else if(ireac.eq.9) then
        print*,' input triton nonlocality range (~0.20 fm) '
       else
        print*,' input deuteron nonlocality range (~0.54 fm) '
       endif
       read*,a(5)
       write(18,*) real(a(5))
       print*,' >>>> ',real(a(5))
      endif
      a(1)=3.1
      a(2)=nil
      a(3)=npw
      a(4)=nil
      write(17,102) (a(i),i=1,5) 
	pnlocs = pnlocs .or. abs(a(5))>1e-10
*---------------------------------------------------------------------
*     entrance channel masses/charges line
*---------------------------------------------------------------------
      a(1)=4.1
      a(2)=m1
      a(3)=a1
      a(4)=zp1
      a(5)=z1
      a(6)=s1
      print*,' target spin in incident channel '
      read*,spin1
      write(18,*) real(spin1)
      print*,' >>>> ',real(spin1)
      a(7)=spin1
      a(8)=0.0
      write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
	potin = ' Numerical  '
	potout= ' Numerical  '
 877  print*,' incident channel potential   '
      print*,'            [1] from those built in   '
      print*,'            [2] specify potential parameters    '
*---------------------------------------------------------------------
      read*,iopti
      if(iopti.lt.1.or.iopti.gt.2) go to 877
      write(18,*) iopti
      print*,' >>>> ',iopti
      print*,'-------------------------------------------------'
      print*,' initial state potential at Elab=',real(energy),' MeV'
      if(iopti.eq.1) then
       if(ireac.eq.1.or.ireac.eq.8) then 
        call proton(energy,a1,z1,step1,nr3max)
       else if(ireac.eq.3) then 
        call neutron(energy,a1,z1,step1,nr3max)
       else if(ireac.eq.7) then
        call helium(energy,a1,z1,step1,nr3max)
       else if(ireac.eq.9) then
        call triton(energy,a1,z1,step1,nr3max)
       else
        call deuteron(energy,a1,z1,ireac)
       endif
      else
       print*,' Central terms '
       print*,' ------------- '
       print*,' Coulomb radius parameter '
       read*,rcd      
       write(18,*) real(rcd)
       print*,' >>>> ',real(rcd)      
       print*,' Real volume  : depth(>0), radius, diffuseness '
       read*,vd,rrd,ard
       write(18,*) real(vd),real(rrd),real(ard)
       print*,' >>>> ',real(vd),real(rrd),real(ard)
       print*,' Imag volume  : depth(>0), radius, diffuseness '
       read*,wd,rid,aid  
       write(18,*) real(wd),real(rid),real(aid)
       print*,' >>>> ',real(wd),real(rid),real(aid)  
*---------------------------------------------------------------------
*      print*,' Imag surface : depth(>0), radius, diffuseness '
*      read*,wisd,risid,aisid   
*---------------------------------------------------------------------
       print*,' Imag surface : depth(>0) '
       read*,wisd    
       write(18,*) real(wisd)
       print*,' >>>> ',real(wisd)    
       print*,' Spin-orbit terms '
       print*,' ---------------- '
       if(ireac.eq.1.or.ireac.eq.8) then
        print*,' proton: coefficients are of L.sigma (~6.0 MeV) '
       else if(ireac.eq.3) then
        print*,' neutron: coefficients are of L.sigma (~6.0 MeV) '
       else if(ireac.eq.7) then
        print*,' 3He: coefficients are of L.sigma (~6.0 MeV) '
       else if(ireac.eq.9) then
        print*,' triton: coefficients are of L.sigma (~6.0 MeV) '
       else 
        print*,' Careful of convention here: non-standard strength '
        print*,' deuteron: half coefficient of L.S (~3.0 MeV)  '
       endif
       print*,' Real s/orbit : depth(>0), radius, diffuseness '
       read*,vsod,rsord,asord      
       write(18,*) real(vsod),real(rsord),real(asord)
       print*,' >>>> ',real(vsod),real(rsord),real(asord)      
       print*,' Imag s/orbit : depth(>0), radius, diffuseness '
       read*,wsod,rsoid,asoid                                   
       write(18,*) real(wsod),real(rsoid),real(asoid)
       print*,' >>>> ',real(wsod),real(rsoid),real(asoid)               
      endif                               
*---------------------------------------------------------------------
*     potential line 1
*---------------------------------------------------------------------
      a(1)=5.1
      a(2)=vd
      a(3)=(wd+wisd)
      a(4)=vsod
      a(5)=wsod
      a(6)=rrd
      a(7)=ard
      a(8)=rcd
      write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     potential line 2
*---------------------------------------------------------------------
      a(1)=6.1
      a(2)=rsord
      a(3)=asord
      a(4)=rsoid
      a(5)=asoid
      write(17,102) (a(i),i=1,5)
*---------------------------------------------------------------------
*     potential line 3
*---------------------------------------------------------------------
      a(1)=7.1
      if(abs(wd+wisd).gt.1.d-10) then
       a(2)=wisd/(wd+wisd)
      else
       a(2)=1.d0
      endif
      a(3)=rid
      a(4)=aid
      a(5)=nil
      a(6)=nil
      write(17,102) (a(i),i=1,6)
*---------------------------------------------------------------------
*     a(1)=8.1
*     a(2)=2.0
*     a(3)=nil
*     a(4)=nil
*     a(5)=nil
*     a(6)=wisd
*     a(7)=risid
*     a(8)=aisid
*     write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     exit channel partial waves/nonlocality line
*---------------------------------------------------------------------
      print*,'================================================= '
      if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then
       print*,' outgoing (deuteron) channel information '
      else if(ireac.eq.2.or.ireac.eq.9) then
       print*,' outgoing (proton) channel information '
      else if(ireac.eq.4) then
       print*,' outgoing (neutron) channel information '
      else if(ireac.eq.5.or.ireac.eq.8) then
       print*,' outgoing (triton) channel information '
      else if(ireac.eq.6) then
       print*,' outgoing (3He) channel information '
      endif
 878  print*,' nonlocality in outgoing channel [1] no  '
      print*,'                                 [2] yes '      
      read*,inonloc
      if(inonloc.lt.1.or.inonloc.gt.2) go to 878
      write(18,*) inonloc
      print*,' >>>> ',inonloc
      a(5)=0.d0 
      if(inonloc.eq.2) then
       if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then
        print*,' input deuteron nonlocality range (~0.54 fm) '
       else if(ireac.eq.2.or.ireac.eq.9) then
        print*,' input proton nonlocality range (~0.85 fm) '
       else if(ireac.eq.4) then
        print*,' input neutron nonlocality range (~0.85 fm) '
       else if(ireac.eq.5.or.ireac.eq.8) then
        print*,' input triton nonlocality range (~0.20 fm) '
       else if(ireac.eq.6) then
        print*,' input 3He nonlocality range (~0.20 fm) '
       endif
       read*,a(5)
       write(18,*) real(a(5))
       print*,' >>>> ',real(a(5))
      endif
      a(1)=3.2
      a(2)=nil
      a(3)=npw
      a(4)=unit
      write(17,102) (a(i),i=1,5) 
	pnlocs = pnlocs .or. abs(a(5))>1e-10
*---------------------------------------------------------------------
*     exit channel masses/charges line
*---------------------------------------------------------------------
      a(1)=4.2
      a(2)=m2
      a(3)=a2
      a(4)=zp2
      a(5)=z2
      a(6)=s2
 879  print*,' target spin in outgoing channel '
      read*,spin2
      ispierr=0
*---------------------------------------------------------------------
*     check consistency of target and transferred angular momenta
*---------------------------------------------------------------------
      if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.5.or.ireac.eq.6.or.
     1  ireac.eq.8) then
       big=spin2+rjtr+0.1
       sma=abs(spin2-rjtr)-0.1
       if(spin1.gt.big.or.spin1.lt.sma) then
	ispierr=1
        print*,' =============================================='
        print*,' input angular momenta are inconsistent:      '
        print*,real(spin2),' +',real(rjtr),' =',real(spin1),' !!'
        print*,' one or more of the target and/or transferred '
        print*,' nucleon angular momenta are wrong.           '
        print*,' =============================================='
       endif
      else
       big=spin1+rjtr+0.1
       sma=abs(spin1-rjtr)-0.1
       if(spin2.gt.big.or.spin2.lt.sma) then
	ispierr=1
        print*,' =============================================='
        print*,' input angular momenta are inconsistent:      '
        print*,real(spin1),' +',real(rjtr),' =',real(spin2),' !!'
        print*,' one or more of the target and/or transferred '
        print*,' nucleon angular momenta are wrong.           '
        print*,' =============================================='
       endif
      endif
      if(ispierr.gt.0) then
      print*,' problem with spins: [1] re-enter outgoing target spin'
      print*,'                     [2] abort and start again      '
      read*,ispierr
      if(ispierr.eq.1) go to 879
      if(ispierr.eq.2) stop
      endif
*---------------------------------------------------------------------
      write(18,*) real(spin2)
      print*,' >>>> ',real(spin2)
      a(7)=spin2
      a(8)=qval
      write(17,102) (a(i),i=1,8) 
      
*------------- FRESCO output of partitions, levels, and entrance potential
	nuclei(1)(4:5)=nucname(nint(zp1)+1)  ! target
	write(nuclei(1)(1:3),'(i3)') m1
	nuclei(2)(4:5)=nucname(nint(z1)+1)  ! target
	write(nuclei(2)(1:3),'(i3)') nint(a1)
* target parity assumed +1 (this choice does not have any effect, as only relative parity important)
        iptyt1 = 1  
        iptyt2 = iptyt1 * (-1)**ltr

	write(16,84)    nuclei(1),real(m1),nint(zp1),
     x                nuclei(2),real(a1),nint(z1),0.0
84	format(/'  &PARTITION namep=''',A5,''' massp=',f8.4,' zp =',i3,
     x        '       nex=1 namet=''',A5,''' masst=',f8.4,' zt =',i3,
     x             ' qval =',f9.4,'/')
      write(16,85)  s1,1,0.0, 1,  spin1,iptyt1,0.0
85	format('  &STATES jp=',f4.1,' ptyp=',i2,' ep=',f8.4,' cpot=',i2,
     x                ' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' /')

	nuclei(3)(4:5)=nucname(nint(zp2)+1)  ! target
	write(nuclei(3)(1:3),'(i3)') m2
	nuclei(4)(4:5)=nucname(nint(z2)+1)  ! target
	write(nuclei(4)(1:3),'(i3)') nint(a2)

 	write(16,84) nuclei(3),real(m2),nint(zp2),
     x             nuclei(4),real(a2),nint(z2),qval
      write(16,85)  s2,1,0.0, 2,  spin2,iptyt2,0.0

	write(16,'(''  &partition /'',/)')

	kp=1
	write(16,86) kp,0,a1,0.0,rcd
86	format('  &pot kp=',i1,' type=',i1,' p(1:3) =',3f9.4,' /')
87	format('  &pot kp=',i1,' type=',i1,' p(1:6) =',6f9.4,' /')
88	format('  &pot kp=',i1,' type=',i1,' shape=',i1,
     x              ' p(1:2) =',2f9.4,' /')
	if(rrd<50..and.rid<50) then
	  write(16,87) kp,1,vd,rrd,ard,wd,rid,aid
	  write(16,87) kp,2,0.0,rrd,ard,wisd,rid,aid  !  risid,aisid
	else if(rid>50.and.rrd<50.) then  ! WS real, numerical imag
	  write(16,86) kp,1,vd,rrd,ard
	  write(16,88) kp,1,8,0.0,1.0  ! read imag central
	  write(15,'(a)') potin//' imag central entrance potential'
	  write(15,*) nr3max,real(step1),real(step1)
	  write(15,'(1p,5e14.6)')(vimag(ii),ii=1,nr3max)
	  write15 = .true.
	else if(rid<50.and.rrd>50.) then  ! WS imag, numerical real
	  write(16,87) kp,1,0.,0.,0.,vd,rrd,ard
	  write(16,88) kp,1,7,1.0,0.0  ! read real central
	  write(15,'(a)') potin//' real central entrance potential'
	  write(15,*) nr3max,real(step1),real(step1)
	  write(15,'(1p,5e14.6)')(vreal(ii),ii=1,nr3max)
	  write15 = .true.
	else  ! both numerical
	  write(16,88) kp,1,9,1.0,1.0  ! read complex central
         write(0,*)'AMM: writing pots'
! Commented by AMoro: adiabatic potentials have not been already calculated here!
!	  write(15,'(a)') potin//' complex central entrance potential'
!	  write(15,*) nr3max,real(step1),real(step1)
!	  write(15,'(1p,6e14.6)')(vreal(ii),vimag(ii),ii=1,nr3max)
!	  write15 = .true.
	endif

	if(rsord<50.) then
	  write(16,87) kp,3,vsod,rsord,asord,wsod,rsoid,asoid
	else  if(rsoid<50.) then
*---------------------------------------------------------------------
*       Modified: 1.0 --> 0.5 as read spin-orbit formfactor is
*       the coefficient of L.S. JAT January 2012
*       Fresco formfactors are multiplied by 2L.S. twofnr reads
*       formfactor and multiplies by L.S
*---------------------------------------------------------------------
	  write(16,88) kp,3,7,0.5,0.0  ! read real 
! Commented by AMoro: adiabatic potentials have not been already calculated here!
!	  write(15,'(a)') potin//' real spin-orbit entrance potential'
!	  write(15,*) nr3max,real(step1),real(step1)
!	  write(15,'(1p,5e14.6)')(vspin(ii),ii=1,nr3max)
	  write15 = .true.
	else
	  write(6,*) ' Numerical imag spin-orbit not implemented'
	 stop
	endif

! AMM
!      stop !!!!!!!!!!!!!!!!!!!!

*---------------------------------------------------------------------
 880  print*,' outgoing channel potential   '
      print*,'            [1] from those built in   '
      print*,'            [2] specify potential parameters    '
*---------------------------------------------------------------------
      read*,iopti
      if(iopti.lt.1.or.iopti.gt.2) go to 880
      write(18,*) iopti
      print*,' >>>> ',iopti
      print*,'-------------------------------------------------   '
      print*,' final state potential at Elab=',real(energy2),' MeV'
      if(iopti.eq.1) then
       if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then 
        call deuteron(energy2,a2,z2,ireac)
       else if(ireac.eq.2.or.ireac.eq.9) then
        call proton(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.4) then
        call neutron(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.5.or.ireac.eq.8) then
        call triton(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.6) then
        call helium(energy2,a2,z2,step2,nr3max)
       endif
      else
       print*,' Central terms '
       print*,' ------------- '
       print*,' Coulomb radius parameter '
       read*,rcd      
       write(18,*) real(rcd)
       print*,' >>>> ',real(rcd)      
       print*,' Real volume  : depth(>0), radius, diffuseness '
       read*,vd,rrd,ard
       write(18,*) real(vd),real(rrd),real(ard)
       print*,' >>>> ',real(vd),real(rrd),real(ard)
       print*,' Imag volume  : depth(>0), radius, diffuseness '
       read*,wd,rid,aid  
       write(18,*) real(wd),real(rid),real(aid)
       print*,' >>>> ',real(wd),real(rid),real(aid)  
*---------------------------------------------------------------------
*      print*,' Imag surface : depth(>0), radius, diffuseness '
*      read*,wisd,risid,aisid   
*---------------------------------------------------------------------
       print*,' Imag surface : depth(>0) '
       read*,wisd   
       write(18,*) real(wisd)
       print*,' >>>> ',real(wisd)   
       print*,' Spin-orbit terms '
       print*,' ---------------- '
       if(ireac.gt.1.and.ireac.ne.3.and.ireac.ne.7) then
        print*,' spin 1/2: coefficients are of L.sigma (~6.0 MeV) '
       else 
        print*,' Careful of convention here: non-standard strength '
        print*,' deuteron: 1/2 coefficient of L.S (~3.0 MeV)'
       endif
       print*,' Real s/orbit : depth(>0), radius, diffuseness '
       read*,vsod,rsord,asord      
       write(18,*) real(vsod),real(rsord),real(asord)
       print*,' >>>> ',real(vsod),real(rsord),real(asord)      
       print*,' Imag s/orbit : depth(>0), radius, diffuseness '
       read*,wsod,rsoid,asoid                                   
       write(18,*) real(wsod),real(rsoid),real(asoid)
       print*,' >>>> ',real(wsod),real(rsoid),real(asoid)   
      endif                               
*---------------------------------------------------------------------
*     potential line 1
*---------------------------------------------------------------------
      a(1)=5.2
      a(2)=vd
      a(3)=(wd+wisd)
      a(4)=vsod
      a(5)=wsod
      a(6)=rrd
      a(7)=ard
      a(8)=rcd
      write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     potential line 2
*---------------------------------------------------------------------
      a(1)=6.2
      a(2)=rsord
      a(3)=asord
      a(4)=rsoid
      a(5)=asoid
      write(17,102) (a(i),i=1,5)
*---------------------------------------------------------------------
*     potential line 3
*---------------------------------------------------------------------
      a(1)=7.2
      if(abs(wd+wisd).gt.1.d-10) then
       a(2)=wisd/(wd+wisd)
      else
       a(2)=1.d0
      endif
      a(3)=rid
      a(4)=aid
      a(5)=nil
      a(6)=nil
      write(17,102) (a(i),i=1,6)
*---------------------------------------------------------------------
*     a(1)=8.2
*     a(2)=2.0
*     a(3)=nil
*     a(4)=nil
*     a(5)=nil
*     a(6)=wisd
*     a(7)=risid
*     a(8)=aisid
*     write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     angles line
*---------------------------------------------------------------------
      b(1)=9.0
      write(17,102) (b(i),i=1,4)
      write(17,102) 
*---------------------------------------------------------------------
*     write entrance channel nucleon JLM potentials if required
*     or write Sao Paulo potential in the 3He case
*     or write Li et al. triton potential in the triton case
*---------------------------------------------------------------------
      if(ijlmp.eq.1.and.(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.8))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(potr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(poti(ii),ii=1,nr3max)
       vreal=potr; vimag=poti; vspin=0.
      endif
      if(ihesp.eq.1.and.ireac.eq.7)then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
       vreal=pspr; vimag=pspi; vspin=0.
      endif
      if(ihesp.eq.1.and.ireac.eq.9)then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(psps(ii),ii=1,nr3max)
       vreal=pspr; vimag=pspi; vspin=psps
      endif
*---------------------------------------------------------------------
*     calculate and write adiabatic or watanabe deuteron potentials 
*     if required
*---------------------------------------------------------------------
      if(iadia.gt.0.or.iwat.gt.0) then
       print*,'-------------------------------------------------'
       if(iadia.gt.0) then
        print*,' Now construct the adiabatic potential: '
	potkind = ' Adiabatic  '
       else if (iwat.gt.0) then
        print*,' Now construct the Watanabe potential: '
	potkind = ' Watanabe   '
       endif
       if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then 
        call adiab(energy2,a2,z2,step2,nr3max)
	potout = potkind
       else if((-1)**ireac.gt.0.or.ireac.eq.5) then
        call adiab(energy ,a1,z1,step1,nr3max)
	potin  = potkind
       endif
      endif
*---------------------------------------------------------------------
*     write exit channel nucleon JLM potentials if required
*     or write Sao Paulo potential in the 3He case
*     or write Li et al. triton potential in the triton case
*---------------------------------------------------------------------
      if(ijlmp.eq.1.and.(ireac.eq.2.or.ireac.eq.4.or.ireac.eq.9))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(potr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(poti(ii),ii=1,nr3max)
       vreal=potr; vimag=poti; vspin=0.
      endif
      if(ihesp.eq.1.and.ireac.eq.6)then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
       vreal=pspr; vimag=pspi; vspin=0.
      endif
      if(ihesp.eq.1.and.(ireac.eq.5.or.ireac.eq.8))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(psps(ii),ii=1,nr3max)
       vreal=pspr; vimag=pspi; vspin=psps
      endif

*------------- FRESCO output of exit potential
	kp=2
	write(16,86) kp,0,a2,0.0,rcd
!	write(16,87) kp,1,vd,rrd,ard,wd,rid,aid
!	write(16,87) kp,2,0.0,rrd,ard,wisd,rid,aid  !  risid,aisid

        if(rrd<50..and.rid<50) then
          write(16,87) kp,1,vd,rrd,ard,wd,rid,aid
          write(16,87) kp,2,0.0,rrd,ard,wisd,rid,aid  !  risid,aisid
        else if(rid>50.and.rrd<50.) then  ! WS real, numerical imag
          write(16,86) kp,1,vd,rrd,ard
          write(16,88) kp,1,8,0.0,1.0  ! read imag central
          write(15,'(a)') potout//' imag central entrance potential'
          write(15,*) nr3max,real(step2),real(step2)
          write(15,'(1p,5e14.6)')(vimag(ii),ii=1,nr3max)
	  write15 = .true.
        else if(rid<50.and.rrd>50.) then  ! WS imag, numerical real
          write(16,87) kp,1,0.,0.,0.,vd,rrd,ard
          write(16,88) kp,1,7,1.0,0.0  ! read real central
          write(15,'(a)') potout//' real central entrance potential'
          write(15,*) nr3max,real(step2),real(step2)
          write(15,'(1p,5e14.6)')(vreal(ii),ii=1,nr3max)
	  write15 = .true.
        else  ! both numerical
          write(16,88) kp,1,9,1.0,1.0  ! read complex central
          write(15,'(a)') potout//' complex central entrance potential'
          write(15,*) nr3max,real(step2),real(step2)
          write(15,'(1p,6e14.6)')(vreal(ii),vimag(ii),ii=1,nr3max)
	  write15 = .true.
        endif


!	write(16,87) kp,3,vsod,rsord,asord,wsod,rsoid,asoid
        if(rsord<50.) then
          write(16,87) kp,3,vsod,rsord,asord,wsod,rsoid,asoid
        else  if(rsoid<50.) then
*---------------------------------------------------------------------
*       Modified: 1.0 --> 0.5 as read spin-orbit formfactor is
*       the coefficient of L.S. JAT January 2012
*       Fresco formfactors are multiplied by 2L.S. twofnr reads
*       formfactor and multiplies by L.S
*---------------------------------------------------------------------
	  write(16,88) kp,3,7,0.5,0.0  ! read real 
          write(15,'(a)') potout//' real spin-orbit potential'
          write(15,*) nr3max,real(step2),real(step2)
          write(15,'(1p,5e14.6)')(vspin(ii),ii=1,nr3max)
	  write15 = .true.
        else
          write(6,*) ' Numerical imag spin-orbit not implemented'
         stop
        endif



*---------------------------------------------------------------------
*     first formfactor line, D0^2 value needed
*---------------------------------------------------------------------
      a(1)=10.0
      a(2)=1.0
      a(3)=nil
      print*,'-------------------------------------------------'
      if(ireac.lt.5) then
       if(ireac.eq.1.or.ireac.eq.2) then
        print*,' <p|d>  vertex constant D0 = -122.5 MeV fm^3/2  '
       else 
        print*,' <n|d>  vertex constant D0 = -122.5 MeV fm^3/2  '
       endif
       print*,' (Reid SC)  e.g. Nucl. Phys. A241 (1975)  36  '
       a(4)=122.5d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3 '
 881   print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 881
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
         print*,' input D0^2 MeV^2 fm^3 '
         read*,a(4)
         write(18,*) a(4)
         print*,' >>>> ',a(4)
       endif
      else if(ireac.eq.5) then
       print*,' <d|t>  vertex constant D0 = -160.0 MeV fm^3/2  '
       print*,'            e.g. Phys. Rev. C 20 (1979) 1631  '
       a(4)=160.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3 '
 882   print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 882
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
         print*,' input D0^2 MeV^2 fm^3 '
         read*,a(4)
         write(18,*) a(4)
         print*,' >>>> ',a(4)
       endif
      else if(ireac.eq.8.or.ireac.eq.9) then
       print*,' <p|t> vertex constant D0 = -469.0 MeV fm^3/2 '
       print*,'            e.g. Phys. Rev. C 4 (1971) 196  '
       a(4)=469.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3 '
 8882  print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 8882
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
         print*,' input D0^2 MeV^2 fm^3 '
         read*,a(4)
         write(18,*) a(4)
         print*,' >>>> ',a(4)
       endif
      else if(ireac.eq.6.or.ireac.eq.7) then
       print*,' <d|3He> vertex constant D0 = -160.0 MeV fm^3/2 '
       print*,'            e.g. Phys. Rev. C 20 (1979) 1631  '
       a(4)=160.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3 '
 883   print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 883
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
         print*,' input D0^2 MeV^2 fm^3 '
         read*,a(4)
         write(18,*) a(4)
         print*,' >>>> ',a(4)
       endif
      endif
      a(5)=nil
	D0 = sqrt(A(4))  ! for FRESCO output later
      if(ireac.eq.8.or.ireac.eq.9)then
      write(17,33) (a(i),i=1,5)
  33  format(0p3f10.4,f10.2,3f10.4)
      goto 3337
      endif
      write(17,102) (a(i),i=1,5)
3337  continue
*---------------------------------------------------------------------
*     finite range correction for light particle vertices line
*---------------------------------------------------------------------
      a(1)=10.01
      print*,'-------------------------------------------------'
      if(ireac.eq.1) then
       print*,' Treatment of range (fnrng) of <p|d> vertex'
      else if(ireac.eq.3) then
       print*,' Treatment of range (fnrng) of <n|d> vertex'
      else if(ireac.eq.2) then
       print*,' Treatment of range (fnrng) of <d|p> vertex'
      else if(ireac.eq.4) then
       print*,' Treatment of range (fnrng) of <d|n> vertex'
      else if(ireac.eq.5) then
       print*,' Treatment of range (fnrng) of <d|t> vertex'
      else if(ireac.eq.6) then
       print*,' Treatment of range (fnrng) of <d|3He> vertex'
      else if(ireac.eq.7) then
       print*,' Treatment of range (fnrng) of <3He|d> vertex'
      else if(ireac.eq.8) then
       print*,' Treatment of range (fnrng) of <p|t> vertex'
       print*,'Glendenning: Ch 9, Nuclear Spectroscopy and'
       print*,'Reactions - recommends zero range for L<= 2'
       print*,'(p,t) or (t,p) reactions.                  '
      else if(ireac.eq.9) then
       print*,' Treatment of range (fnrng) of <t|p> vertex'
       print*,'Glendenning: Ch 9, Nuclear Spectroscopy and'
       print*,'Reactions - recommends zero range for L<= 2'
       print*,'(p,t) or (t,p) reactions.                  '
      endif
 884  print*,'            [1] zero-range (fnrng = 0)          '
      print*,'            [2] local-energy (default value) '
      print*,'            [3] local-energy (specify value) '
      read*,izr
      if(izr.lt.1.or.izr.gt.3) go to 884
      write(18,*) izr
      print*,' >>>> ',izr
      if(izr.eq.1) then
       a(2)=0.0
       print*,' Zero-range calculation (fnrng = 0) '
      else if(izr.eq.2) then
      if(ireac.lt.5) then
       a(2)=0.745712d0
       print*,' Hulthen finite range factor, fnrng=0.745712 '
       print*,' (Reid SC)  e.g. Nucl. Phys. A241 (1975)  36 '
      else if(ireac.eq.5) then
       a(2)=0.746269d0
       print*,' Hulthen finite range factor, fnrng=0.746269 '
       print*,'                 Nucl. Phys. A234 (1974) 301 '
      else if(ireac.eq.8.or.ireac.eq.9) then
       a(2)=0.746269d0
       print*,' Hulthen finite range factor, fnrng=0.746269 '
       print*,' Using same as (d,t),(d,3He),(3He,d) for now!'
      else if(ireac.eq.6.or.ireac.eq.7) then
       a(2)=0.746269d0
       print*,' Hulthen finite range factor, fnrng=0.746269 '
       print*,'                 Nucl. Phys. A234 (1974) 301 '
      endif
      else if(izr.eq.3) then       
       print*,' input local-energy vertex interaction range (fm) '
       print*,' for conventions:  DelVecchio and Daehnick PRC '
       print*,' Vol 6 (1972) p2095 (DVD)                      '
       print*,' [if +ve,  Hulthen, fnrng = R of DVD             ] ' 
       print*,' [fnrng = 1/beta of Nucl. Phys. A241 (1975)  36  ] '
       print*,' [if -ve, Gaussian, fnrng = 1/(2*epsilon) of DVD ] '
       read*,a(2)
       write(18,*) real(a(2))
       print*,' >>>> ',real(a(2))
      endif
	FNRNG = a(2)  ! for FRESCO output later
      write(17,102) (a(i),i=1,2)
*---------------------------------------------------------------------
*     formfactor line
*---------------------------------------------------------------------
      a(1)=10.41
      a(2)=nodes
      a(3)=nil
      if(ireac.eq.3.or.ireac.eq.6) a(3)=z2
      if(ireac.eq.4.or.ireac.eq.7) a(3)=z1
      a(4)=sn
      a(5)=1.0
      if(ireac.eq.8.or.ireac.eq.9) a(5)=2.0
*     set the core mass correctly for stripping/pickup cases
      if(ireac.eq.2.or.ireac.eq.4.or.ireac.eq.7.or.ireac.eq.9) then
       a(6)=a1
      else
       a(6)=a2
      endif
      a(7)=nil
      write(17,102) (a(i),i=1,7)
*---------------------------------------------------------------------
*     formfactor potential line
*---------------------------------------------------------------------
      a(1)=10.42
      print*,'-------------------------------------------------'
      if(ireac.eq.1.or.ireac.eq.2.or.ireac.eq.5) then
       print*,' neutron binding potential  '
      else if(ireac.eq.8.or.ireac.eq.9) then
       print*,' di-neutron binding potential  '
      else
       print*,' proton binding potential  '
      endif
      print*,' radius and diffuseness (e.g. 1.25  0.65 fm)'
      read*,a(2),a(4)
      radius=a(2)
      diffuse=a(4)
      write(18,*) real(a(2)),real(a(4))
      print*,' >>>> ',real(a(2)),real(a(4))
      a(3)=a(2)
      if(ireac.eq.8.or.ireac.eq.9) then
       print*,' Spin-orbit: di-neutron so input strength Vso = 0'
      else
       print*,' Spin-orbit: strength of l.sigma (~6.0 MeV) '
      endif
      read*,a(5)  
      vso=a(5)    
      write(18,*) real(a(5))
      print*,' >>>> ',real(a(5))      
      print*,' Bound state non-locality (0 else ~0.85 fm) '
      read*,a(6)      
      write(18,*) real(a(6))
      print*,' >>>> ',real(a(6))      
*     a(6)=nil                        
      write(17,102) (a(i),i=1,6)
	pnlocs = pnlocs .or. abs(a(6))>1e-10
*------------- FRESCO output of target binding potential
	kp=3
	vd=50.0
	core = min(a1,a2)
	write(16,86) kp,0,core,0.0,a(2)
	write(16,86) kp,1,vd,a(2),a(4)
*---------------------------------------------------------------------
*     formfactor line (for distinct spin-orbit geometry if vso>0)
*---------------------------------------------------------------------
      if(abs(vso).gt.1.d-3) then
       a(1)=10.43
       print*,' Bound state spin-orbit radius parameter  '
       print*,' (if 0 entered will use same geometry as  '
       print*,'  used for the real central interaction)  '
       read*,a(2)      
       write(18,*) real(a(2))
       print*,' >>>> ',real(a(2))      
       if(a(2).lt.0.01) then
	  a(2)=radius
        a(3)=diffuse
	goto 887
       endif
       print*,' Bound state spin-orbit diffuseness parameter  '
       read*,a(3)      
       write(18,*) real(a(3))
       print*,' >>>> ',real(a(3))      
 887   write(17,102) (a(i),i=1,3)
      endif
      write(17,102) 
      print*,'-------------------------------------------------'
      print*
      print*
      
*------------- FRESCO output of target binding potential
	write(16,86) kp,3,vso,a(2),a(3)
	write(16,'(''  &pot /'',/)')

*------------- FRESCO parameters for target overlap function
	nn = nodes+1  ! nn includes node at origin !

	write(16,91) 1,1,2,2,0,nn,ltr,str,rjtr,3,sn,1,0
91	format('  &OVERLAP kn1=',i1,' ic1=',i1,' ic2=',i1,' in=',i1,
     x        ' kind=',i1,' nn=',i1,' l=',i2,' sn=',f5.1,' j=',f6.1,
     x        ' kbpot=',i1,' be=',f9.4,' isc=',i1,' ipc=',i1,' /')
     	write(16,'(''  &overlap /'',/)')

	write(16,92) 2,1,5,D0,FNRNG
92	format('  &COUPLING icto=',i1,' icfrom=',i1,' kind=',i1,
     x        ' p1=',f9.4,' p2=',f9.4,' /')
      write(16,93) -2,1,1,1,1.0
93	format('  &CFP  in=',i2,' ib=',i1,' ia=',i1,' kn=',i1,
     x        ' a=',f9.4,' /')

     	write(16,'(''  &coupling /'',/)')
  
      if(write15) then
 	write(6,94) trim(fname)//'.pots'
94	format(//'  File ',a,' written with external potentials for ',
     x           'FRESCO input'/'   Rename or link as fort.4:'/)

	cmd =  'ln -sf '//trim(fname)//'.pots'//' fort.4'
	call system('echo '//cmd)
	call system(cmd)
	call system('echo link made')
       else
	close(15,status='delete')
       endif

	if(pnlocs) then
	 write(6,*) 
	 write(6,*) ' ****** WARNING: Perey-Buck ',
     x               'nonlocalities NOT implemented in FRESCO ******'
	endif

        write(6,95) trim(fname)//'.fres',trim(fname)//'.out',
     x             trim(fname)//'.tran',trim(fname)//'.brin'
95      format(/ '   Run fresco as:      fresco < ',a,' > ',a,/,
     x          '   Run twofnr7 with input file:  ',a,/,
     x          '   Rerun brush with input file:  ',a,
     x          '  (after changing file name on first line)'/)
      end
*---------------------------------------------------------------------
      subroutine proton(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      common/cjlm/potr(900),poti(900),ijlm,ijlmp
      character nucleon,form*12
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
*---------------------------------------------------------------------
 887  print*,' [1] Bechetti-Greenlees (A>40 20<E<50 MeV)    '
      print*,'            Phys Rev 182 (1969) 1190          '
      print*,' [2] Chapel-Hill 89 Global set (A>40 E>10 MeV)   '
      print*,'            Phys Rep 201 (1991) 57            '
      print*,' [3] Menet (30<E<60) see: ADNTD 17 (1976) p6  '
      print*,' [4] Perey (E<20MeV) see: ADNTD 17 (1976) p6  '
      print*,' [5] JLM microscopic optical potential        '
      print*,'-------------------------------------------------'
*---------------------------------------------------------------------
      ijlmp=0
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.5) go to 887
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       ed=energy
       e=z1
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=54.d0-0.32d0*ed+0.4d0*e/a13+24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.32d0
       api=0.51d0+0.7d0*an
       wpi=11.8d0-0.25d0*ed+12.d0*an
       if(wpi.lt.0.d0) wpi=0.d0
       wvpi=0.22d0*ed-2.7d0
       if(wvpi.lt.0.d0) wvpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV proton energy ')
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=6.2d0
       rsord=1.01d0
       asord=0.75d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
*---------------------------------------------------------------------
*      CH86 parameters
       v0=    52.9d0
       vt=    13.d0
       ve=   -0.3d0
       r0=    1.25d0
       r00=  -0.24d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   10.d0
       wve0=  35.d0
       wvew=  15.d0
       rw=    1.32d0
       rw0=  -0.41d0
       aw=    0.72d0
       ws0=   9.d0
       wst=   14.d0
       wse0=  29.d0
       wsew=  23.d0
*---------------------------------------------------------------------
*      CH89 parameters
       v0=    52.9d0
       vt=    13.1d0
       ve=   -0.299d0
       r0=    1.25d0
       r00=  -0.225d0
       a0=    0.69d0

       rc=    1.24d0
       rc0=   0.12d0

       wv0=   7.8d0 
       wve0=  35.d0
       wvew=  16.d0

       ws0=   10.d0
       wst=   18.d0
       wse0=  36.d0
       wsew=  37.d0

       rw=    1.33d0
       rw0=  -0.42d0
       aw=    0.69d0
*---------------------------------------------------------------------
       a13=a**(1.d0/3.d0)
       rrc=rc*a13+rc0
       rcn=rrc/a13
       ecpp=1.73d0*z/rrc
       erp=e-ecpp
       vrp=v0+vt*((n-z)/a)+erp*ve
       rp=r0*a13+r00
       rpn=rp/a13
       ap=a0
       wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
       if(wvp.lt.0.d0) wvp=0.d0
       rwp=rw*a13+rw0
       rwpn=rwp/a13
       awp=aw
       wsp=(ws0+wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
       if(wsp.lt.0.d0) wsp=0.d0
       print 20,a
       print 101,z,e
   20  format(1h ,' Chapel Hill 89 potentials for a = ',f5.1)
       rcd=rcn
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
       print 12,vrp,rpn,ap,wsp,rwpn,awp,wvp
       vd=vrp
       rrd=rpn
       ard=ap
       wd=wvp
       rid=rwpn
       aid=awp
       wisd=wsp
       risid=rwpn
       aisid=awp
*---------------------------------------------------------------------
*      CH86 parameters
       vsod=5.9d0
       rsord=(1.39d0*a13-1.43)/a13
       asord=0.65d0
*---------------------------------------------------------------------
*      CH89 parameters
       vsod=5.9d0
       rsord=(1.34d0*a13-1.20)/a13
       asord=0.63d0
*---------------------------------------------------------------------
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=49.9d0-0.22*e+26.4*(n-z)/a+0.4*z/a13
       wvp=1.2+0.09*e
       wsp=4.2-0.05*e+15.5*(n-z)/a
       if(wsp.lt.0.d0) wsp=0.d0
       awp=0.74d0-0.008*e+(n-z)/a
       rrd=1.16d0
       ard=0.75d0
       rid=1.37d0
       print 25,a
       print 101,z,e
   25  format(1h ,' Menet potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=6.04d0
       rsord=1.064d0
       asord=0.78d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.4) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=53.3d0-0.55*e+27.0*(n-z)/a+0.4*z/a13
       wvp=0.d0
       wsp=13.5d0
       rrd=1.25d0
       ard=0.65d0
       rid=1.25d0
       awp=0.47d0
       print 24,a
       print 101,z,e
   24  format(1h ,' Perey potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=7.50d0
       rsord=1.25d0
       asord=0.47d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.5) then
       ijlmp=1
       ijlm=1
       print 29,a1
   29  format(1h ,' JLM potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print*,' printout is in',nrmax,' steps of',real(step)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=0.d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       form='(5e14.7)'
       nucleon='p'
       call jlm(potr,poti,nucleon,a1,z1,energy,step,nrmax)
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine neutron(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      common/cjlm/potr(900),poti(900),ijlm,ijlmp
      character nucleon,form*12
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
*---------------------------------------------------------------------
 888  print*,' [1] Bechetti-Greenlees (A>40 20<E<50 MeV)    '
      print*,'            Phys Rev 182 (1969) 1190          '
      print*,' [2] Chapel-Hill 89 Global set (A>40 E>10 MeV)   '
      print*,'            Phys Rep 201 (1991) 57            '
*     print*,' [3] Menet (30<E<60) see: ADNTD 17 (1976) p6  '
*     print*,' [4] Perey (E<20MeV) see: ADNTD 17 (1976) p6  '
      print*,' [3] JLM microscopic optical potential        '
      print*,'-------------------------------------------------'
*---------------------------------------------------------------------
      ijlmp=0
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.3) go to 888
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       ed=energy
       e=z1
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=56.3d0-0.32d0*ed-24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.26d0
       api=0.58d0
       wpi=13.0d0-0.25d0*ed-12.d0*an
       if(wpi.lt.0.d0) wpi=0.d0
       wvpi=0.22d0*ed-1.56d0
       if(wvpi.lt.0.d0) wvpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV neutron energy ')
       rcd=1.25d0
*      print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=6.2d0
       rsord=1.01d0
       asord=0.75d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
*---------------------------------------------------------------------
*      CH86 parameters
       v0=    52.9d0
       vt=    13.d0
       ve=   -0.3d0
       r0=    1.25d0
       r00=  -0.24d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   10.d0
       wve0=  35.d0
       wvew=  15.d0
       rw=    1.32d0
       rw0=  -0.41d0
       aw=    0.72d0
       ws0=   9.d0
       wst=   14.d0
       wse0=  29.d0
       wsew=  23.d0
*---------------------------------------------------------------------
*      CH89 parameters
       v0=    52.9d0
       vt=    13.1d0
       ve=   -0.299d0
       r0=    1.25d0
       r00=  -0.225d0
       a0=    0.69d0

       rc=    1.24d0
       rc0=   0.12d0

       wv0=   7.8d0 
       wve0=  35.d0
       wvew=  16.d0

       ws0=   10.d0
       wst=   18.d0
       wse0=  36.d0
       wsew=  37.d0

       rw=    1.33d0
       rw0=  -0.42d0
       aw=    0.69d0
*---------------------------------------------------------------------
       a13=a**(1.d0/3.d0)
       rrc=rc*a13+rc0
       rcn=rrc/a13
       ecpp=0.d0
       erp=e-ecpp
       vrp=v0-vt*((n-z)/a)+erp*ve
       rp=r0*a13+r00
       rpn=rp/a13
       ap=a0
       wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
       if(wvp.lt.0.d0) wvp=0.d0
       rwp=rw*a13+rw0
       rwpn=rwp/a13
       awp=aw
       wsp=(ws0-wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
       if(wsp.lt.0.d0) wsp=0.d0
       print 20,a
       print 101,z,e
   20  format(1h ,' Chapel Hill 89 potentials for a = ',f5.1)
       rcd=rcn
*      print*,' Coulomb radius parameter = ',real(rcd)
       print 111
       print 12,vrp,rpn,ap,wsp,rwpn,awp,wvp
       vd=vrp
       rrd=rpn
       ard=ap
       wd=wvp
       rid=rwpn
       aid=awp
       wisd=wsp
       risid=rwpn
       aisid=awp
*---------------------------------------------------------------------
*      CH86 parameters
       vsod=5.9d0
       rsord=(1.39d0*a13-1.43)/a13
       asord=0.65d0
*---------------------------------------------------------------------
*      CH89 parameters
       vsod=5.9d0
       rsord=(1.34d0*a13-1.20)/a13
       asord=0.63d0
*---------------------------------------------------------------------
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.98) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=49.9d0-0.22*e+26.4*(n-z)/a+0.4*z/a13
       wvp=1.2+0.09*e
       wsp=4.2-0.05*e+15.5*(n-z)/a
       if(wsp.lt.0.d0) wsp=0.d0
       awp=0.74d0-0.008*e+(n-z)/a
       rrd=1.16d0
       ard=0.75d0
       rid=1.37d0
       print 25,a
       print 101,z,e
   25  format(1h ,' Menet potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=6.04d0
       rsord=1.064d0
       asord=0.78d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.99) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=53.3d0-0.55*e+27.0*(n-z)/a+0.4*z/a13
       wvp=0.d0
       wsp=13.5d0
       rrd=1.25d0
       ard=0.65d0
       rid=1.25d0
       awp=0.47d0
       print 24,a
       print 101,z,e
   24  format(1h ,' Perey potential for a = ',f5.1)
       rcd=1.25d0
*      print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=7.50d0
       rsord=1.25d0
       asord=0.47d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       ijlmp=1
       ijlm=1
       print 29,a1
   29  format(1h ,' JLM potential for a = ',f5.1)
       rcd=1.25d0
*      print*,' Coulomb radius parameter = ',real(rcd)
       print*,' printout is in',nrmax,' steps of',real(step)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=0.d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       form='(5e14.7)'
       nucleon='n'
       call jlm(potr,poti,nucleon,a1,z1,energy,step,nrmax)
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine deuteron(energy,a1,z1,ireac)
      implicit real*8(a-h,o-z) 
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/deut/iadia,iwat
      real*8 vreal(900),vimag(900),vspin(900) ! AMM
      common/fromfold/vreal,vimag,vspin  
*---------------------------------------------------------------------
 889  print*,' Optical potentials for DWBA                 '
      print*,' [1] Lohr-Haeberli (A>40 8<E<13 MeV)         '
      print*,'             see: ADNTD 17 (1976) p6         '
      print*,' [2] Perey-Perey (12<E<25 MeV) no spin-orbit '
      print*,'             see: ADNTD 17 (1976) p6         ' 
      print*,' [3] Daehnick Global (A>27 12<E<90 MeV)      '
      print*,'     Phys. Rev. C 21, 2253 (1980)            ' 
      print*,' [4] Watanabe folding model potential from   '
      print*,'     nucleon potentials and Reid SC deuteron ' 
      print*,'-------------------------------------------------'
      if(ireac.lt.5)then
       print*,' Adiabatic potentials for breakup           '
       print*,' [5] Zero range adiabatic potential         '
       print*,'     Johnson-Soper PRC 1 (1970) 976         '
       print*,' [6] Finite range adiabatic potential       '
       print*,'     Johnson-Tandy NPA 235 (1974) 56        '
       print*,'     using a Reid SC interaction/deuteron   '
       print*,'-------------------------------------------------'
      endif
*---------------------------------------------------------------------
      read*,inopt
      if(ireac.lt.5) then
       if(inopt.lt.1.or.inopt.gt.6) go to 889
      else
       if(inopt.lt.1.or.inopt.gt.4) go to 889
      endif
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=91.13+2.2*z/a13
       wvp=0.d0
       wsp=218.d0/a13/a13
       rrd=1.05d0
       ard=0.86d0
       rid=1.43d0
       awp=0.5d0+0.013*a13*a13
       print 20,a
   20  format(1h ,' LH deuteron potential for a = ',f5.1)
       print 101,z,e
  101  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV deuteron energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
   12  format(1h ,7f7.3)
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=7.0d0/2.d0
       rsord=0.75d0
       asord=0.50d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=81.1-0.22*e+2.0*z/a13
       wvp=0.d0
       wsp=14.4+0.24*e
       rrd=1.15d0
       ard=0.81d0
       rid=1.34d0
       awp=0.68d0
       print 21,a
   21  format(1h ,' P-P deuteron potential for a = ',f5.1)
       print 101,z,e
       rcd=1.15d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=0.0d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       a=a1
       e=energy
       bet=-(e/100.d0)**2
       bet=exp(bet)
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=88.5-0.26*e+0.88*z/a13
       wvp=(12.2+0.026*e)*(1.d0-bet)
       wsp=(12.2+0.026*e)*(bet)
       rrd=1.17d0
       ard=0.709d0+0.0017*e
       rid=1.325d0
       awp=0.53d0+0.07*a13
       awp=awp-0.04*exp(-((  8-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 20-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 28-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 50-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 82-n)/2.d0)**2)
       awp=awp-0.04*exp(-((126-n)/2.d0)**2)
       print 23,a
   23  format(1h ,' Daehnick deuteron potential for a = ',f5.1)
       print 101,z,e
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=Vrp
       wd=Wvp
       aid=awp
       wisd=Wsp
       risid=rid
       aisid=aid
       vsod=(7.33d0-0.029*e)/2.d0
       rsord=1.07d0
       asord=0.66d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.4) then
       iwat=1
       iadia=0
       print 27,a1
   33  format(1h ,' Watanabe deuteron potential for a = ',f5.1)
       print 101,z1,energy
       print*,' from one of the nucleon potentials:    '
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=1.d0
       rsord=99.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
      endif
*---------------------------------------------------------------------
      if(ireac.lt.5.and.inopt.gt.4) then
       iadia=inopt
       iwat=0
       if(iadia.eq.5) then
        print 27,a1
   27   format(1h ,' ZR Adiabatic deuteron potential for a = ',f5.1)
       else 
        print 57,a1
   57   format(1h ,' FR Adiabatic deuteron potential for a = ',f5.1)
       endif
       print 101,z1,energy
       print*,' from one of the nucleon potentials:    '
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=1.d0
       rsord=99.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine adiab(energy,a,z,step,nrmax)
      implicit real*8(a-h,o-z)
      real*8 vreal(900),vimag(900),vspin(900)
      real*8 prn(900),pin(900),prp(900),pip(900),pson(900),psop(900)
      real*8 pgrid(900)
      common/folded/prn,pin,pson,prp,pip,psop,pgrid
      common/fromfold/vreal,vimag,vspin
      character form*12,nucleon
      common/cjlm/potr(900),poti(900),ijlm,ijlmp
      common/djlm/rlr,rli
      common/deut/iadia,iwat
*---------------------------------------------------------------------
*     potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*---------------------------------------------------------------------
*     latter is coefficient of l.sigma for nucleons (ww~6MeV)
*---------------------------------------------------------------------
*     we are here if iadia > 0 or iwat > 0
*---------------------------------------------------------------------
      do ii=1,900
       prn(ii)= 0.d0
       pin(ii)= 0.d0
       pson(ii)=0.d0
       prp(ii)= 0.d0
       pip(ii)= 0.d0
       psop(ii)=0.d0
       vreal(ii)=0.d0
       vimag(ii)=0.d0
       vspin(ii)=0.d0
      enddo
      form='(5e14.7)'
*---------------------------------------------------------------------
 890  print*,' [1] Bechetti-Greenlees (A>40 20<E<50 MeV)    '
      print*,'            Phys Rev 182 (1969) 1190          '
      print*,' [2] Chapel-Hill 89 Global set (A>40 E>10 MeV)   '
      print*,'            Phys Rep 201 (1991) 57            '
      print*,' [3] JLM microscopic optical potentials       '
      print*,'-------------------------------------------------'
*---------------------------------------------------------------------
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.3) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
*---------------------------------------------------------------------
*      use half the deuteron energy
*---------------------------------------------------------------------
       ed=energy/2.d0
       e=z
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=54.d0-0.32d0*ed+0.4d0*e/a13+24.d0*an
       vnr=56.3d0-0.32d0*ed-24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.32d0
       rni=1.26d0
       api=0.51d0+0.7d0*an
       ani=0.58d0
       wpi=11.8d0-0.25d0*ed+12.d0*an
       wni=13.0d0-0.25d0*ed-12.d0*an
       if(wpi.lt.0.d0) wpi=0.d0
       if(wni.lt.0.d0) wni=0.d0
       wvpi=0.22d0*ed-2.7d0
       wvni=0.22d0*ed-1.56d0
       if(wvpi.lt.0.d0) wvpi=0.d0
       if(wvni.lt.0.d0) wvni=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV nucleon energy ')
       print 11
       print 111
   11  format(1h ,'   proton ')
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
       print 19
       print 191
   19  format(1h ,'   neutron')
  191  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vnr,rr,ar,wni,rni,ani,wvni
   12  format(1h ,7f7.3)
       vso=6.2d0
       rso=1.01d0
       aso=0.75d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vso,rso,aso
*---------------------------------------------------------------------
       rr =rr *a13
       rpi=rpi*a13
       rni=rni*a13
       rso=rso*a13
*      for /folded/prn,pin,pson,prp,pip,psop
       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
        prn(ii)=ws(r,vnr,rr,ar)
        pin(ii)=ws(r,wvni,rni,ani)+wsd(r,wni,rni,ani)
*       pson(ii)=(2.d0+a)/(4.d0*(1.d0+a))*wso(r,vso,rso,aso)
        pson(ii)=wso(r,vso,rso,aso)
        prp(ii)=ws(r,vpr,rr,ar)
        pip(ii)=ws(r,wvpi,rpi,api)+wsd(r,wpi,rpi,api)
        psop(ii)=pson(ii)
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
        endif 
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)

! AMoro
        write(15,'(a)') ' complex central entrance potential'
        write(15,*) nrmax,real(step1),real(step1)
        write(15,'(1p,6e14.6)')(vreal(ii),vimag(ii),ii=1,nrmax)

        write(15,'(a)') ' real spin-orbit entrance potential'
	 write(15,*) nrmax,real(step1),real(step1)
	 write(15,'(1p,5e14.6)')(vspin(ii),ii=1,nrmax)
!        write15 = .true.
       endif
      endif 
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       e=energy/2.d0
       n=nint(a-z)
*---------------------------------------------------------------------
*      CH86 parameters
       v0=    52.9d0
       vt=    13.d0
       ve=   -0.3d0
       r0=    1.25d0
       r00=  -0.24d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   10.d0
       wve0=  35.d0
       wvew=  15.d0
       rw=    1.32d0
       rw0=  -0.41d0
       aw=    0.72d0
       ws0=   9.d0
       wst=   14.d0
       wse0=  29.d0
       wsew=  23.d0
*---------------------------------------------------------------------
*      CH89 parameters
       v0=    52.9d0
       vt=    13.1d0
       ve=   -0.299d0
       r0=    1.25d0
       r00=  -0.225d0
       a0=    0.69d0

       rc=    1.24d0
       rc0=   0.12d0

       wv0=   7.8d0 
       wve0=  35.d0
       wvew=  16.d0

       ws0=   10.d0
       wst=   18.d0
       wse0=  36.d0
       wsew=  37.d0

       rw=    1.33d0
       rw0=  -0.42d0
       aw=    0.69d0
*---------------------------------------------------------------------
       a13=a**(1.d0/3.d0)
       rrc=rc*a13+rc0
       rcn=rrc/a13
       ecpp=1.73d0*z/rrc
       ecnn=0.d0
       erp=e-ecpp
       ern=e-ecnn
       vrp=v0+vt*((n-z)/a)+erp*ve
       vrn=v0-vt*((n-z)/a)+ern*ve
       rp=r0*a13+r00
       rpn=rp/a13
       rn=rpn
       ap=a0
       an=ap
       wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
       wvn=wv0/(1.d0+exp((wve0-ern)/wvew))
       if(wvp.lt.0.d0) wvp=0.d0
       if(wvn.lt.0.d0) wvn=0.d0
       rwp=rw*a13+rw0
       rwpn=rwp/a13
       rwn=rwpn
       awp=aw
       awn=awp
       wsp=(ws0+wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
       wsn=(ws0-wst*((n-z)/a))/(1.d0+exp((ern-wse0)/wsew))
       if(wsp.lt.0.d0) wsp=0.d0
       if(wsn.lt.0.d0) wsn=0.d0
       print 20,a
       print 101,z,e
   20  format(1h ,' Chapel Hill 89 potentials for a = ',f5.1)
       print 11
       print 111
       print 12,vrp,rpn,ap,wsp,rwpn,awp,wvp
       print 19
       print 191
       print 12,vrn,rn,an,wsn,rwn,awn,wvn
*---------------------------------------------------------------------
*      CH86 parameters
       vso=5.9d0
       rso=(1.39d0*a13-1.43)/a13
       aso=0.65d0
*---------------------------------------------------------------------
*      CH89 parameters
       vso=5.9d0
       rso=(1.34d0*a13-1.20)/a13
       aso=0.63d0
*---------------------------------------------------------------------
       print 112
       print 12,vso,rso,aso
*---------------------------------------------------------------------
       rr =rpn*a13
       rpi=rwpn*a13
       rni=rwn*a13
       rso=rso*a13
*      for /folded/prn,pin,pson,prp,pip,psop

       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
        prn(ii)=ws(r,vrn,rr,ap)
        pin(ii)=ws(r,wvn,rni,awn)+wsd(r,wsn,rni,awn)
*       pson(ii)=(2.d0+a)/(4.d0*(1.d0+a))*wso(r,vso,rso,aso)
        pson(ii)=wso(r,vso,rso,aso)
        prp(ii)=ws(r,vrp,rr,ap)
        pip(ii)=ws(r,wvp,rpi,awp)+wsd(r,wsp,rpi,awp)
        psop(ii)=pson(ii)
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
        endif 
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)

        do ii=1,nrmax
         r=ii*step
         write(35,1017) r,vreal(ii),vimag(ii),vspin(ii)
 1017    format(2x,f5.2,3(2x,d12.5))
        enddo
       endif
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       if(ijlm.eq.0) ijlm=1
       ed=energy/2.d0
       nucleon='n'
       call jlm(prn,pin,nucleon,a,z,ed,step,nrmax)
       nucleon='p'
       call jlm(prp,pip,nucleon,a,z,ed,step,nrmax)
*      for /folded/prn,pin,pson,prp,pip,psop
       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
        pson(ii)=0.d0
        psop(ii)=pson(ii)
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
        endif 
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)
       endif
      endif

! AMoro
      if(iadia.eq.5) then
        write(15,'(a)') ' complex central entrance potential'
        write(15,*) nrmax,real(step),real(step)
        write(15,'(1p,6e14.6)')(vreal(ii),vimag(ii),ii=1,nrmax)

        write(15,'(a)') ' real spin-orbit entrance potential'
	 write(15,*) nrmax,real(step),real(step)
	 write(15,'(1p,5e14.6)')(vspin(ii),ii=1,nrmax)
!        write15 = .true.
      endif

*---------------------------------------------------------------------
*     if finite range folding to be done. nucleon potentials to
*     be used have been calculated above
      if(iwat.eq.1.or.iadia.eq.6) then
       call folder(a,step,nrmax)
*      print the returned folded potentials
       print*,'printout is in',nrmax,' steps of',real(step)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)

! AMoro
        write(15,'(a)') ' complex central entrance potential'
        write(15,*) nrmax,real(step),real(step)
        write(15,'(1p,6e14.6)')(vreal(ii),vimag(ii),ii=1,nrmax)

        write(15,'(a)') ' real spin-orbit entrance potential'
	 write(15,*) nrmax,real(step),real(step)
	 write(15,'(1p,5e14.6)')(vspin(ii),ii=1,nrmax)
      endif 
*---------------------------------------------------------------------
      return
      end
*---------------------------------------------------------------------
      subroutine folder(amt,step,nrmax)
      implicit real*8(a-h,o-z)
      real*8 prn(900),pin(900),prp(900),pip(900),pson(900),psop(900)
      real*8 vreal(900),vimag(900),vspin(900),pgrid(900)
      dimension g(8,400),vtr(8,250),funct(8,200),gd(400)
      dimension wri(200),xri(200),sint(8),pint(8)
      common/folded/prn,pin,pson,prp,pip,psop,pgrid
      common/fromfold/vreal,vimag,vspin
c     ----------------------------------------------------------------
      ndim=900
      icalc=1
      if(iadia.gt.0) icalc=2
c     ----------------------------------------------------------------
c     icalc=1 watanabe    icalc=2 adiabatic    
c     ----------------------------------------------------------------
c     mass of target nucleus
c     write(6,23)amt
c     ----------------------------------------------------------------
      do ii=1,8
       do jj=1,250
        vtr(ii,jj)=0.d0
       enddo
      enddo
c     ----------------------------------------------------------------
      mmts=96
      call gauss(-1.d0,1.d0,mmts,xri,wri)
c-----------------------------------------------------------------------
c     read(5,*) ndata,sep,ndatar,sepr
c     small r - internal   big R - external 
c-----------------------------------------------------------------------
      ndata=201
      if(iwat.gt.0) then
       sep=0.1d0
      else if(iadia.gt.0) then
       sep=0.05d0
      endif
c-----------------------------------------------------------------------
      ndatar=200
      sepr=step
c-----------------------------------------------------------------------
      do 9999 jj=1,2
c     ----------------------------------------------------------------
c     loop over neutron followed by proton potentials
c     ----------------------------------------------------------------
      do 1001 k=1,ndatar
      rcap=k*sepr
      do ii=1,8
       sint(ii)=0.0
      enddo
      deno=0.0
      j21=2
c     ----------------------------------------------------------------
      do 1002 j=1,ndata
      r=(j-1)*sep
      do ii=1,8
       pint(ii)=0.0
      enddo
      do i=1,mmts
       um=xri(i)
       p2=0.5d0*(3.d0*um**2-1.0)
       r1=sqrt(rcap**2+0.25*r**2+rcap*r*um)
c      will use function terp(r,fun,rgrid,npts,ndim) to interpolate
c     ----------------------------------------------------------------
       if(jj.eq.1) then
        gr= terp(r1,prn ,pgrid,nrmax,ndim)
        gi= terp(r1,pin ,pgrid,nrmax,ndim)
        gso=terp(r1,pson,pgrid,nrmax,ndim)
       else if(jj.eq.2) then
        gr= terp(r1,prp ,pgrid,nrmax,ndim)
        gi= terp(r1,pip ,pgrid,nrmax,ndim)
        gso=terp(r1,psop,pgrid,nrmax,ndim)
       endif
c     ----------------------------------------------------------------
       funct(1,i)=gr
       funct(2,i)=gi
       funct(3,i)=0.d0
       funct(4,i)=p2*gr
       funct(5,i)=p2*gi
       funct(6,i)=gso
       funct(7,i)=gso*um
       funct(8,i)=p2*gso
c     ----------------------------------------------------------------
       do ii=1,8
        pint(ii)=pint(ii)+funct(ii,i)*wri(i)
       enddo
c     ----------------------------------------------------------------
      enddo
c     ----------------------------------------------------------------
      call fact(r,br,br2,br3,br4,icalc)
c     central  -------------------------------------------------------
      do ii=1,3
       g(ii,j)=br*pint(ii)/2.0
      enddo
c     tensor tr  -----------------------------------------------------
      do ii=4,5
       g(ii,j)=3.d0/sqrt(2.d0)*br2*pint(ii)/2.d0
      enddo
c     spin-orbit  ----------------------------------------------------
      g(6,j)=(br3+br4/3.0)*pint(6)/2.d0
      g(6,j)=g(6,j)+(r/(2.0*rcap))*br3*pint(7)/2.d0
      g(6,j)=g(6,j)-br4/3.0*pint(8)/2.0
c     overall factor of one-half in Keaton + Armstrong
      g(6,j)=g(6,j)/2.d0
c     denominator  ---------------------------------------------------
      gd(j)=br
c     ----------------------------------------------------------------
      if(j-j21.le.0) go to 1022
      j1=j-1
      j2=j-2
c     ----------------------------------------------------------------
      do ii=1,6
       sint(ii)=sint(ii)+sep*(g(ii,j2)+4.d0*g(ii,j1)+g(ii,j))/3.0
      enddo
c     ----------------------------------------------------------------
      deno=deno+sep*(gd(j2)+4.0*gd(j1)+gd(j))/3.0
      j21=j+1
 1022 continue
 1002 continue
c     ----------------------------------------------------------------
      do ii=1,6
       vtr(ii,k)=sint(ii)/deno+vtr(ii,k)
      enddo
c     ----------------------------------------------------------------
      if(jj.eq.2) then
       vreal(k)=vtr(1,k)
       vimag(k)=vtr(2,k)
       vspin(k)=vtr(6,k)
       write(36,1017) rcap,vtr(1,k),vtr(2,k),vtr(6,k),vtr(4,k),vtr(5,k)
      endif
 1017 format(2x,f5.2,5(2x,d12.5))
 1001 continue
 9999 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine saxon(r,rt,deet,g)
      implicit real*8 (a-h,o-z)
      g=1.0+exp((r-rt)*deet)
      g=1.0/g
      return
      end
c-----------------------------------------------------------------------
      subroutine dsaxon(r,rt,deet,dg)
      implicit real*8 (a-h,o-z)
      gg=exp((r-rt)*deet)
      g=1+gg
      g=1.0/g
      dg=4.0*gg*g*g
      return
      end
c-----------------------------------------------------------------------
      subroutine thomas(r,rt,deet,dg)
      implicit real*8 (a-h,o-z)
      dg=0.0
      if(r.eq.0.0) return
      gg=exp((r-rt)*deet)
      g=1+gg
      g=1.0/g
      dg=2.0*gg*g*g*deet/r
      return
      end
c-----------------------------------------------------------------------
      subroutine fact(r,br,br2,br3,br4,icalc)
      implicit real*8 (a-h,o-z)
      if(r.eq.0) go to 1
       u=us(r)
       w=ud(r)
      if(icalc.eq.1)then
       v0=u
       v2=w
      else
      v0=vs(r)
      v2=vd(r)
      endif
   18 br=(u*v0+w*v2)
      br2=(u-w/sqrt(2.d0))*v2+w*v0
      br3=u*v0-w*v2-0.5d0*(w*v0+u*v2)/sqrt(2.d0)
      br4=1.5d0*((w*v0+u*v2)/sqrt(2.d0)+w*v2)
      return
    1 continue
      br =0.d0
      br2=0.d0
      br3=0.d0
      br4=0.d0
      return
      end
c-----------------------------------------------------------------------
      function vs(r)
      implicit real*8 (a-h,o-z)
      if(r.eq.0.0) r=1.0d-10
      vs=(vc(r)*us(r)+2.82842712*vt(r)*ud(r))
      return
      end
c-----------------------------------------------------------------------
      function vd(r)
      implicit real*8 (a-h,o-z)
      if(r) 10,10,20
   10 vd=0.0
      return
   20 vd=((vc(r)-2.*vt(r)-3.*vls(r))*ud(r)+2.82842712*vt(r)*us(r))
      return
      end
c-----------------------------------------------------------------------
      function vc(r)
      implicit real*8 (a-h,o-z)
      x=0.7*r
      e=exp(-x)
      vc=e*(-10.463+e*(105.468+e*e*(-3187.8+9924.3*e*e)))/x
      return
      end
c-----------------------------------------------------------------------
      function vt(r)
      implicit real*8 (a-h,o-z)
      x=0.7*r
      y=1.0/x
      e=exp(-x)
      vt=e*y*(-10.463*(1.+3.*y*(1.+y)-3.*y*(4.+y)*e**3)+
     1e**3*(351.77-1673.5*e*e))
      return
      end
c-----------------------------------------------------------------------
      function vls(r)
      implicit real*8 (a-h,o-z)
      x=0.7*r
      e=exp(-x)
      vls=e**4*(708.91-2713.1*e*e)/x
      return
      end
c-----------------------------------------------------------------------
      function us(r)
      implicit real*8 (a-h,o-z)
      dimension x(33),y(33),v(33)
      data x,y/1.000d-2,4.125d-2,7.250d-2,1.350d-1,1.975d-1,
     1         2.600d-1,3.225d-1,3.850d-1,4.475d-1,5.100d-1,
     2         5.725d-1,6.350d-1,6.975d-1,7.600d-1,8.850d-1,
     3         1.010d-0,1.135d-0,1.260d-0,1.385d-0,1.510d-0,
     4         1.760d-0,2.010d-0,2.510d-0,3.010d-0,3.510d-0,
     5         4.010d-0,4.510d-0,5.010d-0,5.510d-0,6.010d-0,
     6         7.010d-0,8.010d-0,9.010d-0,
     7         0.0000d-0,3.3373d-5,2.3901d-4,2.7621d-3,1.2737d-2,
     8         3.6062d-2,7.5359d-2,1.2847d-1,1.8993d-1,2.5349d-1,
     9         3.1390d-1,3.6770d-1,4.1317d-1,4.4992d-1,4.9953d-1,
     a         5.2406d-1,5.3166d-1,5.2864d-1,5.1926d-1,5.0621d-1,
     b         4.7505d-1,4.4200d-1,3.7864d-1,3.2249d-1,2.7399d-1,
     c         2.3251d-1,1.9719d-1,1.6718d-1,1.4172d-1,1.2012d-1,
     d         8.6290d-2,6.1983d-2,4.4523d-2/
      data v/2.9751d-4,2.6127d-3,1.2335d-2,8.2951d-2,2.5326d-1,
     1  5.0052d-1,7.5072d-1,9.3349d-1,1.0162d-0,1.0034d-0,
     2  9.2042d-1,7.9674d-1,6.5744d-1,5.1985d-1,2.8494d-1,
     3  1.1865d-1,1.1397d-2,-5.4090d-2,-9.2523d-2,-1.1418d-1,
     4  -1.3097d-1,-1.3193d-1,-1.2004d-1,-1.0453d-1,-8.9709d-2,
     5  -7.6510d-2,-6.5054d-2,-5.5229d-2,-4.6850d-2,-3.9726d-2,
     6  -2.8547d-2,-2.0508d-2,-1.4731d-2/
      z=0.7*r
      if(z-0.01)10,10,20
   10 us=0.0d0
      return
   20 if(z-9.01) 40,40,30
   30 us=0.87758d0*exp(-.33087d0*z)
      return
   40 do 100 i=1,33
      if(x(i)-z) 70,80,100
   70 if(x(i+1)-z) 100,60,50
  100 continue
   80 us=y(i)
      return
   60 us=y(i+1)
      return
   50 h=x(i+1)-x(i)
      p=(z-x(i))/h
      y1=y(i)
      y2=y(i+1)
      v1=v(i)
      v2=v(i+1)
      a1=y1
      a2=h*v1
      a3=3.d0*(y2-y1)-h*(2.d0*v1+v2)
      a4=2.d0*(y1-y2)+h*(v1+v2)
      us=((a4*p+a3)*p+a2)*p+a1
      return
      end
c-----------------------------------------------------------------------
      function ud(r)
      implicit real*8 (a-h,o-z)
      dimension x(33),y(33),v(33)
      data x,y/1.000d-2,4.125d-2,7.250d-2,1.350d-1,1.975d-1,
     1         2.600d-1,3.225d-1,3.850d-1,4.475d-1,5.100d-1,
     2         5.725d-1,6.350d-1,6.975d-1,7.600d-1,8.850d-1,
     3         1.010d-0,1.135d-0,1.260d-0,1.385d-0,1.510d-0,
     4         1.760d-0,2.010d-0,2.510d-0,3.010d-0,3.510d-0,
     5         4.010d-0,4.510d-0,5.010d-0,5.510d-0,6.010d-0,
     6         7.010d-0,8.010d-0,9.010d-0,
     8         0.0000d-0,1.0850d-5,8.4073d-5,1.0369d-3,4.9642d-3,
     9         1.4446d-2,3.0795d-2,5.3157d-2,7.8995d-2,1.0525d-1,
     a         1.2933d-1,1.4958d-1,1.6529d-1,1.7645d-1,1.8710d-1,
     b         1.8654d-1,1.7946d-1,1.6910d-1,1.5742d-1,1.4553d-1,
     c         1.2314d-1,1.0373d-1,7.3859d-2,5.3293d-2,3.9077d-2,
     d         2.9115d-2,2.2016d-2,1.6871d-2,1.3079d-2,1.0243d-2,
     f         6.4412d-3,4.1575d-3,2.7363d-3/
      data v/7.4202d-5,8.9321d-4,4.4755d-3,3.1982d-2,1.0121d-1,
     1  2.0596d-1,3.1477d-1,3.9371d-1,4.2466d-1,4.0842d-1,
     2  3.5772d-1,2.8848d-1,2.1428d-1,1.4427d-1,3.3294d-2,
     3  -3.5899d-2,-7.3151d-2,-9.0091d-2,-9.5299d-2,-9.4158d-2,
     4  -8.3973d-2,-7.1282d-2,-4.9327d-2,-3.3940d-2,-2.3612d-2,
     5  -1.6691d-2,-1.2002d-2,-8.7768d-3,-6.5201d-3,-4.9136d-3
     6  ,-2.8956d-3,-1.7758d-3,-1.1227d-3/
      z=0.7*r
      if(z-0.01)10,10,20
   10 ud=0.0d0
      return
   20 if(z-9.01) 40,40,30
   30 z=0.33087d0*z
      w=1.d0/z
      ud=0.023013d0*exp(-z)*(1.d0+3.d0*w*(1.d0+w))
      return
   40 do 100 i=1,33
      if(x(i)-z) 70,80,100
   70 if(x(i+1)-z) 100,60,50
  100 continue
   80 ud=y(i)
      return
   60 ud=y(i+1)
      return
   50 h=x(i+1)-x(i)
      p=(z-x(i))/h
      y1=y(i)
      y2=y(i+1)
      v1=v(i)
      v2=v(i+1)
      a1=y1
      a2=h*v1
      a3=3.d0*(y2-y1)-h*(2.d0*v1+v2)
      a4=2.d0*(y1-y2)+h*(v1+v2)
      ud=((a4*p+a3)*p+a2)*p+a1
      return
      end
*---------------------------------------------------------------------
      subroutine triton(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/che3/pspr(900),pspi(900),psps(900),ihesp
*---------------------------------------------------------------------
*     potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*---------------------------------------------------------------------
*     latter is the coefficient of l.sigma 
*---------------------------------------------------------------------
 890  print*,' [1] Bechetti-Greenlees (not well determined)    '
      print*,'                      see: ADNTD 17 (1976) p6    ' 
      print*,' [2] X. Li et al. Global potential E<40 MeV      '    
      print*,'                      see: NPA 789 (2007) 103    '
      print*,' [3] D.Y. Pang et al. GDP08                      '
      print*,'                      see: PRC 79 (2009) 024615  '
      print*,'-------------------------------------------------'
*---------------------------------------------------------------------
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.3) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
      ihesp=0
      a=a1
      ed=energy
      e=z1
      an=(a-2.d0*e)/a
      a13=a**0.3333333333d0
      z13=e**0.3333333333d0
      if(inopt.eq.1) then
       vpr=165.d0-0.17d0*ed-6.4d0*an
       rr=1.20d0
       ar=0.72d0
       rpi=1.40d0
       api=0.84d0
       wvpi=46.0d0-0.33d0*ed-110.d0*an
       if(wvpi.lt.0.d0) wvpi=0.d0
       wpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees 3h potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV triton energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=2.5d0
       rsord=1.20d0
       asord=0.72d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
      if(inopt.eq.2) then
*---------------------------------------------------------------------
*     X Li et al. Global potential NPA 2007
*---------------------------------------------------------------------
       ihesp=1
       do i=1,900
        pspr(i)=0.d0
        pspi(i)=0.d0
        psps(i)=0.d0
       enddo
       print 29,a1
   29  format(1h ,' Li et al. triton potential for a = ',f5.1)
       rcd=1.4219d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=1.d0
       rsord=99.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
*---------------------------------------------------------------------
       vlir=137.6d0-0.1456d0*ed+0.0436d0*ed*ed+4.3751d0*an
       vlir=vlir+1.0474d0*e/a13
       wlis=37.06d0-0.6451d0*ed-47.19d0*an
       wliv=7.383d0+0.5025d0*ed-0.0097d0*ed*ed
       vliso=1.9029d0
       rlir= (1.12010d0-0.1504d0/a13)
       rliv= (1.32020d0-0.1776d0/a13)
       rlis= (1.25100d0-0.4622d0/a13)
       rliso=(0.46991d0+0.1294d0/a13)
       alir= (0.68330d0+0.01910d0*a13)
       aliv= (1.11900d0+0.01913d0*a13)
       alis= (0.81140d0+0.01159d0*a13)
       aliso=(0.35450d0-0.05220d0*a13)
       print 117
  117  format(1h ,'    vr     ro     ao     ws     ris     ais  ')
       print 12,vlir,rlir,alir,wlis,rlis,alis
       print 119
  119  format(1h ,'    wv     riv     aiv   ')
       print 12,wliv,rliv,aliv
       print 118
  118  format(1h ,'    vso    rso    aso ')
       print 12,vliso,rliso,aliso
       rlir= rlir*a13
       rliv= rliv*a13
       rlis= rlis*a13
       rliso=rliso*a13
       print*,' printout is in',nrmax,' steps of',real(step)
*---------------------------------------------------------------------
       do ii=1,nrmax
        r=ii*step
        pspr(ii)=ws(r,vlir,rlir,alir)
        pspi(ii)=ws(r,wliv,rliv,aliv)+wsd(r,wlis,rlis,alis)
        psps(ii)=wso(r,vliso,rliso,aliso)
       enddo
      endif
      if(inopt.eq.3) then
*---------------------------------------------------------------------
*     GDP08 potential for 3H  dypang
*---------------------------------------------------------------------
       ap=3.0d0
       zp=1.0d0
       call gdp08(ap,zp,a1,z1,energy)
      endif
*---------------------------------------------------------------------
      return
      end
*---------------------------------------------------------------------
      subroutine helium(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      character spfile*30,guff*10
      real*8 funr(1500),funi(1500),rgrid(1500)
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/che3/pspr(900),pspi(900),psps(900),ihesp
      character form*12
*---------------------------------------------------------------------
 890  print*,' [1] Bechetti-Greenlees (not well determined)    '
      print*,'                      see: ADNTD 17 (1976) p6    ' 
      print*,' [2] Read Sao Paulo potential from file          '
      print*,'                      see: Reference             ' 
      print*,' [3] D.Y. Pang et al. GDP08                      '
      print*,'                      see: PRC 79 (2009) 024615  '
      print*,'-------------------------------------------------'
      ndim=1500
*---------------------------------------------------------------------
      ihesp=0
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.3) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       ed=energy
       e=z1
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=151.9d0-0.17d0*ed+50.d0*an
       rr=1.20d0
       ar=0.72d0
       rpi=1.40d0
       api=0.88d0
       wvpi=41.7d0-0.33d0*ed+44.d0*an
       if(wvpi.lt.0.d0) wvpi=0.d0
       wpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees 3he potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV 3He  energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=2.5d0
       rsord=1.20d0
       asord=0.72d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
      if(inopt.eq.2) then
*---------------------------------------------------------------------
*     the externally read Sao Paulo potential option requested
*---------------------------------------------------------------------
       ihesp=1
       print 29,a1
   29  format(1h ,' Sao Paulo 3He potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print*,' printout is in',nrmax,' steps of',real(step)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=0.d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       form='(5e14.7)'
*---------------------------------------------------------------------
       print*,' filename with Sao Paulo potential ?  '
       read '(a)',spfile
       open(22,file=spfile,status='unknown')
       print '(a,a)','  >>>> ',spfile
       write(18,'(a)') spfile
*---------------------------------------------------------------------
*      read things from the external file on channel 22
*      format is to be provided here (JAT 14/3/2008)
*---------------------------------------------------------------------
       read(22,'(a)') guff
       read(22,*) npts,spstep,spstart
       do ii=1,npts
	rgrid(ii)=spstart+(ii-1)*spstep
        read(22,*) funr(ii),funi(ii)
*       write(31,*)rgrid(ii),funr(ii),funi(ii)
       enddo
       rmax=rgrid(npts)
*---------------------------------------------------------------------
*      interpolate the potential to the internal working grid
*      c/i/so (nrmax points with step) first radial point at step
*---------------------------------------------------------------------
       do ii=1,nrmax
        r=ii*step
	if(r.lt.rmax) then
         pspr(ii)=terp(r,funr,rgrid,npts,ndim)
         pspi(ii)=terp(r,funi,rgrid,npts,ndim)
	else
    	 pspr(ii)=0.d0
	 pspi(ii)=0.d0
	endif
*       diagnostic prints
*       write(30,*) r,pspr(ii),pspi(ii)
       enddo
      endif
      if(inopt.eq.3) then
*---------------------------------------------------------------------
*     GDP08 potential for 3He  dypang
*---------------------------------------------------------------------
       ap=3.0d0
       zp=2.0d0
       call gdp08(ap,zp,a1,z1,energy)
      endif
*---------------------------------------------------------------------
      return
      end
*-------------------------------------------------------------
      subroutine gdp08(ap,zp,at,zt,Einc)
      implicit real*8(a-h,o-z)
      real*8 V, rv, av
      real*8 Wv, rrwv, aawv, Ws, rrws, aaws
      real*8 VVso, rvso, avso, WWso, rwso, awso
      real*8 RC
      real*8 V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0
      real*8 WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT
      real*8 WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST
      real*8 VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP
      real*8 WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP
      real*8 RC0,RCA, RCAP
      real*8 AT, ZT, AP, ZP, Einc
      real*8 EC, Epot, Ecm
      real*8 varpsilon
      real*8 MDA,MDB,MDC,VMD
      real*8 PI
      real*8 AT13, zero
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/iunit/iunit
      common/GDP08Para/
     &        V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0,
     &        WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT,
     &        WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST,
     &        VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP,
     &        WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP,
     &        RC0, RCA, RCAP
      PI=4.0d0*datan(1.0d0)
      call initGDP08
c     real part
      V0=118.25d0
      VE=-0.12512d0
      R0=1.3007
      R0A =-0.4816
      A0 =0.8148
      VC=1.0
c     Wv
      WV0=38.481
      RWV=1.3120
      RWVA=-0.1290
      AWV=0.8399
      WVE0=156.09
      WVEW=52.442
c     Ws
      WS0=35.037
      WST=34.181
      RWS=1.3120
      RWSA=-0.1290
      AWS=0.8399
      WSE0=30.755
      WSEW=106.36
c     Coulomb
      RC0=1.238
      RCA=0.116
c     isovector parameter
      if(ZP.EQ.2.) then
        varpsilon=(AT-2.0d0*ZT)/AT
      elseif(ZP.EQ.1.) then
        varpsilon=(2.0d0*ZT-AT)/AT
      else
        varpsilon=0.0d0
      endif
c     Coulomb correction, asumming alpha=1.0 as in CH89
      if(RCAP.EQ.0.0) RCAP=-1.0d0/3.0d0
      Ecm=AT*Einc/(AT+AP)
      RC=RC0+RCA*AT**RCAP
      if(VC .EQ. 1.0) then
        EC =(1.73d0/RC)*ZT*ZP/AT**(1.d0/3.d0)
      else
        EC = 0.0d0
      endif
      Epot=Einc-EC
      AT13=AT**(1.0d0/3.0d0)
      zero=0.0d0
c     real part
      V  = V0 + VE*Epot + (VT+VTE*Epot)*varpsilon
      rv = R0 + (R0A+R0AE*Epot)/AT13
      av = A0
c     volume imag
      if(WVE0 .LE. 0.0) then
        Wv  = (WV0+WVE*Epot) + (WVT+WVTE*Epot)*varpsilon
      else
        Wv  = (WV0+WVT*varpsilon)/(1.0d0+dexp(-(Epot-WVE0)/WVEW))
      endif
      rrwv= RWV + RWVA/AT13
      aawv= AWV + AWVT*varpsilon
c     surface imag
      if(WSE0 .LE. 0.0) then
        Ws = (WS0+WSE*Epot) + (WST+WSTE*Epot)*varpsilon
      else
        Ws = (WS0+WST*varpsilon)/(1.0d0+dexp( (Epot-WSE0)/WSEW))
      endif
      rrws= RWS + RWSA/AT13
      aaws= AWS + AWST*varpsilon
c     spin-orbit real
      VVso = VSO + VSOE*Epot + VSOA*AT**VSAP
      rvso = RSO + RSOA*AT13
      avso = ASO
c     spin-orbit imag
      WWso = WSO + WSOE*Epot + WSOA*AT**WSAP
      rwso = RWO + RWOA*AT13
      awso = AWO
c     Coulomb radius
      RC = RC0 +RCA*AT**RCAP
      vd   = V
      rrd  = RV
      ard  = AV
      wd   = Wv
      rid  = rrwv
      aid  = aawv
      wisd = Ws
      risid= rrws
      aisid= aaws
      vsod = VVso
      rsord= rvso
      asord= avso
      wsod = WWso
      rsoid= rwso
      asoid= awso
      return
      end
*-------------------------------------------------------------
      subroutine initGDP08
      implicit none
      real*8 V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0
      real*8 WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT
      real*8 WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST
      real*8 VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP
      real*8 WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP
      real*8 RC0,RCA, RCAP
      common/GDP08Para/
     &        V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0,
     &        WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT,
     &        WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST,
     &        VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP,
     &        WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP,
     &        RC0, RCA, RCAP

      V0  = 0.0d0
      VE  = 0.0d0
      VT  = 0.0d0
      VTE = 0.0d0
      R0  = 0.0d0
      R0A = 0.0d0
      R0Ap= 0.0d0
      A0  = 0.0d0
      VC  = 1.0d0
      R0AE= 0.0d0
      MDAE= 0.0d0
      MDA0= 0.0d0

      WV0 = 0.0d0
      WVE = 0.0d0
      WVT = 0.0d0
      WVTE= 0.0d0
      RWV = 0.0d0
      RWVA= 0.0d0
      RWVP= 0.0d0
      AWV = 0.0d0
      WVC = 1.0d0
      WVE0= 0.0d0
      WVEW= 0.0d0
      AWVT= 0.0d0

      WS0 = 0.0d0
      WSE = 0.0d0
      WST = 0.0d0
      WSTE= 0.0d0
      RWS = 0.0d0
      RWSA= 0.0d0
      RWSP= 0.0d0
      AWS = 0.0d0
      WSC = 1.0d0
      WSE0= 0.0d0
      WSEW= 0.0d0
      AWST= 0.0d0

      VSO = 0.0d0
      VSOE= 0.0d0
      RSO = 1.2d0
      RSOA= 0.0d0
      RSOP= 0.0d0
      ASO = 0.6d0
      VSOA= 0.0d0
      VSAp= 0.0d0

      WSO = 0.0d0
      WSOE= 0.0d0
      RWO = 1.2d0
      RWOA= 0.0d0
      RWOP= 0.0d0
      AWO = 0.6d0
      WSOA= 0.0d0
      WSAp= 0.0d0

      RC0 = 0.0d0
      RCA = 0.0d0
      RCAP= 0.0d0
      return
      end
*---------------------------------------------------------------------
      subroutine jlm(potr,poti,nucleon,aa,zz,Elab,step,nrmax)
*     ------------------------------------------------------------
*     JLM local density approximation (LDA) N+A potentials          
*     JAT February 1999 - tested March 1999.
*     ------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/density/rad(500),rhon(500),rhop(500),rhom(500)
      real*8 a(3,3),b(3,3),c(3,3),d(4,4),f(4,4),nn
      real*8 xrt(200),wrt(200),xmu(200),wmu(200)
      real*8 potr(900),poti(900)
*     real*8 gbden(401)
      character dname*12,guff*4
*     character fmt*20
      character nucleon,ans,nucl*12
      common/arrays/a,b,c,d,f
      common /dcom/ rho0, rad0, diff, beta2, b2bar
      common/djlm/rlr,rli
      pi=4.d0*datan(1.d0)
      do i=1,900
       potr(i)=0.d0
       poti(i)=0.d0
      enddo
      print*,' ---------------------------------------------------'
      print*,' JLM local density approximation nucleon potentials'
      print*
*     ------------------------------------------------------------
*     print*,' Output file trailer for potentials:'
*     read '(a)',fname
*     open(17,file='jlm.'//fname,status='unknown')
*     open(22,file='density',status='unknown')
*     open(27,file='jlmf.'//fname,status='unknown')
*     print '(a,a)','  Output is to file: ','jlm.'//fname
*     ------------------------------------------------------------
*  16 print*,' neutron (n) or proton (p) optical potential'  
*     read '(a)',nucleon
*     if(nucleon.ne.'n'.and.nucleon.ne.'p') goto 16
      nuc=-1
      nucl='proton'
      if(nucleon.eq.'n') then
       nuc=1
       nucl='neutron'
      endif
      print '(a,2x,a)','  Optical potential for ',nucl
*     ------------------------------------------------------------
*     print*,' Target mass A and charge Z '
*     read*,aa,zz
*     print9,aa,zz
      nn=aa-zz
      aa13=aa**(1.d0/3.d0)
      if(nuc.eq.-1) rc=1.123d0*aa13+2.35d0/aa13-2.07d0/aa
*     print*,' Lab Energy (MeV) '
*     read*,E
*     read*,Elab
*     print9,Elab
*     compute cm energy
      E=Elab*aa/(aa+1.d0) 
*     print'(a,f10.4)',' incident cm energy (MeV) =  ',E
*     ------------------------------------------------------------
*     print*,' potential at radii: rmin,rmax,rstep '
*     read*,rmin,rmax,rstep
*     ------------------------------------------------------------
      rmin=0.d0
      rmax=nrmax*step
      rstep=step
   14 print*,' Uses the JLM parameterisation of Bauge '
*     ------------------------------------------------------------
*     print9,rmin,rmax,rstep
    9 format(3f10.4)
*  14 print*,' Choose JLM potential parameter set '
*     print*,'   1) High energy set (>15 MeV) '
*     print*,'   2) Low energy set  (<15 MeV) '
*     read*,iset
*     print*,iset
*     if(iset.ne.1.and.iset.ne.2) goto 14
      iset=1
*     if(Elab.lt.15.d0) iset=2 
*  11 print*,' Version of local density approximation '
*     print*,'   1) LDA (rx=rt)'
*     print*,'   2) LDA (rx=rp)'
*     print*,'   3) LDA (Mid-point) ******** '
*     read*,ilda
*     print*,ilda
*     if(ilda.lt.1.or.ilda.gt.3) goto 11
      print*,' Uses the mid-point LDA prescription '
      ilda=3
   17 print*,' Choose assumed target density     '
      print*,'   1) Negele (Fermi) form          '
      print*,'   2) Specify rms radius           '
      print*,'   3) Oscillator form              '
      print*,'   4) 3pF three parameter Fermi    '
      print*,'   5) read Alex Brown HF densities '
*     print*,'   6) read other density           '
      read*,icho
      if(icho.lt.1.or.icho.gt.5) go to 17
      print*,' >>>> ',icho
      write(18,*) icho
      if(abs(icho).gt.5) goto 17
      if(icho.eq.2) then
        print*,' ---------------------------------------------------'
        print*,' input rms radius (fm) '
        read*,rms
        print*,' >>>> ',rms
	write(18,*) rms
        print*,' ---------------------------------------------------'
   18   print*,'   1) Gaussian density '
        print*,'   2) Woods-Saxon density'
        read*,irms
        if(irms.lt.1.or.irms.gt.2) go to 18
        print*,' >>>> ',irms
	write(18,*)irms
        if(irms.ne.1.and.irms.ne.2) goto 18
        if(irms.eq.2) then
         diff=0.54d0
         print*,'  default diffuseness (0.54 fm) (y/n) '
         read '(a)',ans
         print '(a,a)',' >>>> ',ans
         write(18,'(a)') ans
         if(ans.eq.'n') then
          print*,'  enter diffuseness '
          read*,diff
          print*,' >>>> ',diff
	  write(18,*)diff
         endif
         call finder(aa,rms,rho0,rad0,diff)
         iden=3
         print*,' ---------------------------------------------------'
         print*,' iden = 3: WooSax:  rho0 = ',rho0
         print*,' rad0 = ',rad0,' diff = ',diff
        else
         iden=2
         gamma=sqrt(2.d0/3.d0)*rms
         rho0=aa/(sqrt(pi)*gamma)**3
         print*,' ---------------------------------------------------'
         print*,' iden = 2: Gauss :  gamma = ',gamma
         print*,' rho0 = ',rho0 
        endif       
       else if(icho.eq.1) then
*      Negele Woods-Saxon parameters
        iden=1
        diff=0.54d0
        rad0=(0.978d0+0.0206d0*aa13)*aa13
        rho0=3.d0*aa/(4.d0*pi*rad0**3*(1.d0+(pi*diff/rad0)**2))
	vol=aa
	call volrms(vol,rms)
        print*,' ---------------------------------------------------'
        print*,' iden = 1: Negele:  rho0 = ',rho0
        print*,' rad0 = ',rad0,' diff = ',diff
	print*,' rms mass radius =',rms
      else if(icho.eq.3) then
        print*,' Input a and alfa values in '
        print*,' (1+alfa*[r/a]**2)exp(-r**2/a**2) '
        read*,rad0,diff
        print*,' >>>> ',rad0,diff
        write(18,*) rad0,diff
        iden=4
        call finder2(aa,rms,rho0,rad0,diff)
        print*,' ---------------------------------------------------'
        print*,' iden = 4: Oscill:  rho0 = ',rho0
        print*,' a    = ',rad0,' alfa = ',diff
        print*,' rho0 = ',rho0,' rms  = ',rms 
      else if(icho.eq.4) then
        print*,' Input rad, diff and w values in '
        print*,' (1+w*[r/rad]**2)/(1+exp(r-rad)/diff)  '
        read*,rad0,diff,www
        print*,' >>>> ',rad0,diff,www
        write(18,*) rad0,diff,www
        iden=5
        call finder3(aa,rms,rho0,rad0,diff,www)
        print*,' ---------------------------------------------------'
        print*,' iden = 5: 3pFermi:  rho0 = ',rho0
        print*,' rad  = ',rad0,' diff = ',diff
        print*,'  w =   ',www
        print*,' rho0 = ',rho0,' rms  = ',rms 
*     ------------------------------------------------------------
*       print to file for reading by smat 
*       drin=0.05d0
*       nval=401
*       istart=0
*       fmt='(4d19.8)'
*       write(22,10) fmt,drin,nval,istart
*  10   format(a20,f10.4,2i10)
*       do iir=1,nval
*        rr=(iir-1)*drin
*        gbden(iir)=oscden(rr,rho0,rad0,diff)
*        write(19,*) rr,gbden(iir)
*       enddo
*       write(22,'(4(d19.8))') (gbden(iir),iir=1,nval)
*     ------------------------------------------------------------
      else if(icho.eq.5.or.icho.eq.6) then
*       density is read from file
        iden=6
        ijax=500
        print*,' using read matter density '
        print*,' ---------------------------------------------------'
        print*,' iden = 6: read matter density  '
        if(icho.eq.5) print*,' f(r)=Hartree-Fock density'
        if(icho.eq.6) print*,' f(r)=external density'
        print*
        print*,' file with required density is'
        read '(a)',dname
        print '(a)','  >>>> '//dname
        write(18,'(a)') dname
        open(19,file=dname,status='unknown')
        print*
*       following for reading of Hartree-Fock densities
        if(icho.eq.5)then
         rewind 19
         ico=0
         read(19,'(a)') guff
         read(19,'(a)') guff
         read(19,'(a)') guff
         do ijj=1,1000
          read(19,*,end=707) rad(ijj),rhop(ijj),rhon(ijj),rhom(ijj)
          ico=ico+1
         enddo
  707    continue
         close(19)
         ival=ico
         drin=rad(2)-rad(1)
         print*,' number of read HF radii ',ival
         print*,' with step ',drin
*     ------------------------------------------------------------
*        do ijj=1,ival
*         write(30,*) rad(ijj),rhom(ijj)
*        enddo
*     ------------------------------------------------------------
        endif
        if(icho.eq.6)then
         rewind 19
         ico=0
         do ijj=1,1000
          read(19,*,end=708) rad(ijj),rhom(ijj)
          ico=ico+1
         enddo
  708    continue
         close(19)
         ival=ico
         drin=rad(2)-rad(1)
         print*,' number of read radii ',ival
         print*,' with step ',drin
*     ------------------------------------------------------------
*        do ijj=1,ival
*         write(31,*) rad(ijj),rhom(ijj)
*        enddo
*     ------------------------------------------------------------
        endif
      endif
      print*,' ---------------------------------------------------'
*     ------------------------------------------------------------
*     look up potential parameters
      call assign(iset)
      alpha=nuc*(nn-zz)/aa
*     ------------------------------------------------------------
*     print*,' Maximum NN relative separation (rptmax) for folding'
*     read*,rptmax
*     ------------------------------------------------------------
      rptmax=4.d0
*     ------------------------------------------------------------
*     print9,rptmax
*     print*,' Quadrature points: radial and cos(theta) integrals'
*     read*,mqr,mqmu
*     ------------------------------------------------------------
      mqr=32
      mqmu=48
*     ------------------------------------------------------------
*     print*,mqr,mqmu
      zero=0.d0
      one=1.d0
*     for integration over radius of target density 
      call gauss(zero,rptmax,mqr,xrt,wrt)
*     for integration over cos(theta) 
      call gauss(-one,one,mqmu,xmu,wmu)
*     ------------------------------------------------------------
*     Gaussian effective interaction strength (for t=1.0 fm)
*     and for unit integrated strength
*     print*,' real and imaginary Gaussian folding ranges t'
*     read*,tr,ti
*     print9,tr,ti
*     ------------------------------------------------------------
      tr=1.d0
      ti=1.d0
*     ------------------------------------------------------------
      t2r=tr*tr
      t2i=ti*ti
      rmagr=1.d0/(sqrt(pi)*tr)**3
      rmagi=1.d0/(sqrt(pi)*ti)**3
*     ------------------------------------------------------------
      print*,' real and imaginary potential scalings lambda '
      if(rlr.lt.1.d-3.and.rli.lt.1.d-3) then
       print*,' systematics usually suggest  1.0  0.8 '
       read*,rlr,rli
       print*,' >>>> ',rlr,rli
       write(18,*) rlr,rli
       print*
      else
       print*,' using  ',rlr,rli
       print*
      endif 
*     ------------------------------------------------------------
      i=0
! AMM
!      do rp=rstep,rmax+0.01d0,rstep
!      i=i+1
      nrp=nint((rmax-rstep)/rstep)+1
      do i=1,nrp
      potr(i)=0.d0
      poti(i)=0.d0
      rp2=rp*rp
      do imqr=1,mqr
       rpt=xrt(imqr)
       rpt2=rpt*rpt
       gg3r=g3(rmagr,t2r,rpt2)*rpt2
       gg3i=g3(rmagi,t2i,rpt2)*rpt2
       summr=0.d0
       summi=0.d0
        do imqmu=1,mqmu      
        bval=rp*rpt*xmu(imqmu)
        rt2=rpt2+rp2-2.d0*bval
        rt=sqrt(rt2)     
        if(ilda.eq.1) then
         rx=rt
        elseif(ilda.eq.2) then
         rx=rp
        elseif(ilda.eq.3) then
         rx=sqrt(rpt2/4.d0+rp2-bval)
        endif
       if(iden.eq.1) then
        rhot=woosax(rt,rho0,rad0,diff)
        rhox=woosax(rx,rho0,rad0,diff)
       elseif(iden.eq.2) then
        rhot=gauden(rt,rho0,gamma)
        rhox=gauden(rx,rho0,gamma)
       elseif(iden.eq.3) then
        rhot=woosax(rt,rho0,rad0,diff)
        rhox=woosax(rx,rho0,rad0,diff)
       elseif(iden.eq.4) then
        rhot=oscden(rt,rho0,rad0,diff)
        rhox=oscden(rx,rho0,rad0,diff)
       elseif(iden.eq.5) then
        rhot=fermi3(rt,rho0,rad0,diff,www)
        rhox=fermi3(rx,rho0,rad0,diff,www)
       elseif(iden.eq.6) then
        rhot=1.d-15
        if(rt.lt.rad(ival)) rhot=terp(rt,rhom,rad,ival,ijax)
        rhox=1.d-15
        if(rx.lt.rad(ival)) then
           rhoxn=terp(rx,rhon,rad,ival,ijax)
           rhoxp=terp(rx,rhop,rad,ival,ijax)
           rhox=rhoxn+rhoxp
           alpha=nuc*(rhoxn-rhoxp)/rhox
*          alpha=nuc*(nn-zz)/aa
        endif
       endif
*       if proton take care of Coulomb interaction
        Ecal=E
        if(nuc.eq.-1) then
         if(rx.ge.rc) then
          vc=1.44d0*zz/rx
         else
          vc=0.72d0*zz/rc*(3.d0-(rx/rc)**2) 
         endif
         Ecal=E-vc
        endif
        call pots(Ecal,rhox,iset,V0,W0,V1,W1,rmtilde)
        summr=summr+wmu(imqmu)*rhot*(V0+alpha*V1)/rhox
        summi=summi+wmu(imqmu)*rhot*rmtilde*(W0+alpha*W1)/rhox
       enddo
       potr(i)=potr(i)+wrt(imqr)*summr*gg3r
       poti(i)=poti(i)+wrt(imqr)*summi*gg3i
      enddo
      potr(i)=2.d0*pi*potr(i)*rlr
      poti(i)=2.d0*pi*poti(i)*rli
*     print*,i
*     write(24,*) rp,potr(i),poti(i)
      enddo     
*     ------------------------------------------------------------
*     call print1(17,aa,zz,nuc,Elab)
*     call print2(17)
*     write(17,'(5(e14.7))') (potr(i),i=2,301) 
*     call print2(17)
*     write(17,'(5(e14.7))') (poti(i),i=2,301)
*     ------------------------------------------------------------
*     for cupid
*     ------------------------------------------------------------
*     write(18,'(a)')' 31.       1.0   1.20    0.61    0.       0.0  '
*     write(18,'(a)')'150.0   0.0'
*     write(18,'(5(e16.7))') (potr(i),i=2,151) 
*     write(18,'(a)')' 31.       0.0   1.20    0.61    0.       1.0  '
*     write(18,'(a)')'150.0   1.0'
*     write(18,'(5(e16.7))') (poti(i),i=2,151) 
*     write(18,'(a)') ' 0'
*     ------------------------------------------------------------
*     for fresco
*     ------------------------------------------------------------
*     write(27,'(a)')'201   0.1    0.0'
*     write(27,'(5(e16.7))') (potr(i),i=1,201) 
*     write(27,'(a)')'201   0.1    0.0'
*     write(27,'(5(e16.7))') (poti(i),i=1,201) 
*     ------------------------------------------------------------
      return
      end
*     ------------------------------------------------------------
      subroutine pots(E,rho,iset,V0,W0,V1,W1,rmtilde)
      implicit real*8(a-h,o-z)
      real*8 a(3,3),b(3,3),c(3,3),d(4,4),f(4,4),ImN
      common/arrays/a,b,c,d,f
      ef=fermi(E,rho,iset)  
      V0=0.d0
      ReN=0.d0
      rmtilde=0.d0
      do i=1,3
       do j=1,3
        con=rho**i*E**(j-1)
        V0=V0+a(i,j)*con
        ReN=ReN+b(i,j)*con
        rmtilde=rmtilde+c(i,j)*con
       enddo
      enddo
      rmtilde=1.d0-rmtilde
      V0P=0.d0
      do i=1,3
       do j=2,3
        con=(j-1)*rho**i*E**(j-2)
        V0P=V0P+a(i,j)*con
       enddo
      enddo
      rmstar=1.d0-V0P  
      W0=0.d0
      ImN=0.d0
      do i=1,4
       do j=1,4
        con=rho**i*E**(j-1)
        W0=W0+d(i,j)*con
        ImN=ImN+f(i,j)*con
       enddo
      enddo
*     dd=100.d0
*     if(iset.eq.1) dd=600.d0
*     dd=625.d0
*     changed according to Pang/Bauge 2011
      dd=126.25d0
      W0=W0/(1.d0+dd/(E-ef)**2)
      ImN=ImN/(1.d0+1.d0/(E-ef))
      rmbar=rmstar/rmtilde
      V1=rmtilde*ReN
      W1=ImN/rmbar
      return
      end
*     ------------------------------------------------------------
      subroutine assign(iset)
      implicit real*8(a-h,o-z)
      real*8 a(3,3),b(3,3),c(3,3),d(4,4),f(4,4)
      common/arrays/a,b,c,d,f
      a(1,1)=-0.9740d+3
      a(1,2)= 0.1126d+2
      a(1,3)=-0.4250d-1
      a(2,1)= 0.7097d+4
      a(2,2)=-0.1257d+3
      a(2,3)= 0.5853d+0
      a(3,1)=-0.1953d+5
      a(3,2)= 0.4180d+3
      a(3,3)=-0.2054d+1
     
      b(1,1)= 0.3601d+3
      b(1,2)=-0.5224d+1
      b(1,3)= 0.2051d-1
      b(2,1)=-0.2691d+4
      b(2,2)= 0.5130d+2
      b(2,3)=-0.2470d+0 
      b(3,1)= 0.7733d+4
      b(3,2)=-0.1717d+3
      b(3,3)= 0.8846d+0

      c(1,1)= 0.4557d+1
      c(1,2)=-0.5291d-2
      c(1,3)= 0.6108d-5
      c(2,1)=-0.2051d+1
      c(2,2)=-0.4906d+0
      c(2,3)= 0.1812d-2
      c(3,1)=-0.6509d+2
      c(3,2)= 0.3095d+1   
      c(3,3)=-0.1190d-1  

*     if(iset.eq.1) then

      d(1,1)=-0.1483d+4
      d(1,2)= 0.3718d+2
      d(1,3)=-0.3549d+0
      d(1,4)= 0.1119d-2
      d(2,1)= 0.2988d+5
      d(2,2)=-0.9318d+3
      d(2,3)= 0.9591d+1
      d(2,4)=-0.3160d-1
      d(3,1)=-0.2128d+6
      d(3,2)= 0.7209d+4
      d(3,3)=-0.7752d+2
      d(3,4)= 0.2611d+0
      d(4,1)= 0.5125d+6
      d(4,2)=-0.1796d+5
      d(4,3)= 0.1980d+3
      d(4,4)=-0.6753d+0
              
      f(1,1)= 0.5461d+3
      f(1,2)=-0.1120d+2
      f(1,3)= 0.1065d+0
      f(1,4)=-0.3541d-3
      f(2,1)=-0.8471d+4
      f(2,2)= 0.2300d+3
      f(2,3)=-0.2439d+1
      f(2,4)= 0.8544d-2
      f(3,1)= 0.5172d+5
      f(3,2)=-0.1520d+4
      f(3,3)= 0.1717d+2
      f(3,4)=-0.6211d-1
      f(4,1)=-0.1140d+6
      f(4,2)= 0.3543d+4
      f(4,3)=-0.4169d+2
      f(4,4)= 0.1537d+0

*     else
*     These from JLM

      d(1,1)=-0.5138d+3
      d(1,2)=-0.2985d+2
      d(1,3)= 0.1452d+1
      d(1,4)= 0.9428d-1
      d(2,1)= 0.9078d+4
      d(2,2)= 0.5757d+3
      d(2,3)=-0.3435d+2
      d(2,4)=-0.2310d+1
      d(3,1)=-0.6192d+5
      d(3,2)=-0.4155d+4
      d(3,3)= 0.2657d+3
      d(3,4)= 0.1882d+2
      d(4,1)= 0.1516d+6
      d(4,2)= 0.1037d+5
      d(4,3)=-0.6748d+3
      d(4,4)=-0.5014d+2
              
      f(1,1)= 0.6597d+3
      f(1,2)= 0.4509d+1
      f(1,3)=-0.2383d+1
      f(1,4)=-0.4324d-1
      f(2,1)=-0.1263d+5
      f(2,2)=-0.6572d+2
      f(2,3)= 0.5866d+2
      f(2,4)= 0.1348d+1
      f(3,1)= 0.9428d+5
      f(3,2)= 0.5972d+3
      f(3,3)=-0.4923d+3
      f(3,4)=-0.1295d+2
      f(4,1)=-0.2453d+6
      f(4,2)=-0.1800d+4
      f(4,3)= 0.1358d+4
      f(4,4)= 0.3836d+2

*     These from Bauge (old, Bauge paper)

      d(1,1)=-0.6599d+3
      d(1,2)= 0.1077d+2
      d(1,3)=-0.7886d-1
      d(1,4)= 0.1875d-3
      d(2,1)= 0.1144d+5
      d(2,2)=-0.2908d+3
      d(2,3)= 0.2443d+1
      d(2,4)=-0.6203d-2
      d(3,1)=-0.7451d+5
      d(3,2)= 0.2207d+4
      d(3,3)=-0.1993d+2
      d(3,4)= 0.5175d-1
      d(4,1)= 0.1761d+6
      d(4,2)=-0.5458d+4
      d(4,3)= 0.5113d+2
      d(4,4)=-0.1339d+0
              
      f(1,1)= 0.4596d+3
      f(1,2)=-0.6440d+1
      f(1,3)= 0.4040d-1
      f(1,4)=-0.9009d-4
      f(2,1)=-0.7693d+4
      f(2,2)= 0.1464d+3
      f(2,3)=-0.1025d+1
      f(2,4)= 0.2337d-2
      f(3,1)= 0.5525d+5
      f(3,2)=-0.1112d+4
      f(3,3)= 0.7967d+1
      f(3,4)=-0.1802d-1
      f(4,1)=-0.1437d+6
      f(4,2)= 0.3038d+4
      f(4,3)=-0.2220d+2
      f(4,4)= 0.5026d-1

*     These from Bauge (from MOM code, Pang 2011)

      d(1,1)=-0.65986d+3
      d(1,2)= 0.10768d+2
      d(1,3)=-0.78863d-1
      d(1,4)= 0.18755d-3
      d(2,1)= 0.11437d+5
      d(2,2)=-0.29076d+3
      d(2,3)= 0.24430d+1
      d(2,4)=-0.62028d-2
      d(3,1)=-0.74505d+5
      d(3,2)= 0.22068d+4
      d(3,3)=-0.19926d+2
      d(3,4)= 0.51754d-1
      d(4,1)= 0.17609d+6
      d(4,2)=-0.54579d+4
      d(4,3)= 0.51127d+2
      d(4,4)=-0.13386d+0

      f(1,1)= 0.45959d+3
      f(1,2)=-0.64399d+1
      f(1,3)= 0.40403d-1
      f(1,4)=-0.90086d-4
      f(2,1)=-0.76929d+4
      f(2,2)= 0.14639d+3
      f(2,3)=-0.10244d+1
      f(2,4)= 0.23367d-2
      f(3,1)= 0.55250d+5
      f(3,2)=-0.11121d+4
      f(3,3)= 0.79667d+1
      f(3,4)=-0.18008d-1
      f(4,1)=-0.14373d+6
      f(4,2)= 0.30382d+4
      f(4,3)=-0.22202d+2
      f(4,4)= 0.50258d-1

*     endif
      return  
      end       
c--------------------------------------------------------------- 
      real*8 function fermi(E,rho,iset)
      implicit real*8(a-h,o-z)
*     Bauge modification PRC 58 page 1120
*     E0=10.d0
*     change accirding to Pang/Bauge 2011
      E0=9.d0
      ae=2.d0
      fermih=rho*(-510.8d0+3222.d0*rho-6250.d0*rho*rho)
      fermil=-22.d0-rho*(298.52d0-3760.23d0*rho+
     +         12345.82d0*rho*rho)
      fwt=1.d0/(1.d0+exp((E-E0)/ae))
      fermi=fwt*fermil+(1.d0-fwt)*fermih
      return
      end
c--------------------------------------------------------------- 
      real*8 function woosax(r,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      con=(r-rad0)/diff
      woosax=rho0/(1.d0+exp(con))
      return
      end
c--------------------------------------------------------------- 
      real*8 function fermi3(r,rho0,rad0,diff,www)
      implicit real*8(a-h,o-z)
      con=(r-rad0)/diff
      con1=(r/rad0)**2
      fermi3=rho0*(1+www*con1)/(1.d0+exp(con))
      return
      end
c--------------------------------------------------------------- 
      subroutine finder3(aa,rms,rho0,rad0,diff,www)
      implicit real*8 (a-h,o-z)
      rho0=1.d0
      call vol3pf(vol,rms,rho0,rad0,diff,www)
      rho0=aa/vol
      call vol3pf(vol,rms,rho0,rad0,diff,www)
      print 24, vol, rms
 24   format('  Volume: ',f12.7,' rms ',f12.7)
      return
      end 
c--------------------------------------------------------------- 
      real*8 function oscden(r,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      con=(r/rad0)**2
      oscden=rho0*(1+diff*con)*exp(-con)
      return
      end
c--------------------------------------------------------------- 
      real*8 function gauden(r,rho0,gamma)
      implicit real*8(a-h,o-z)
      gauden=rho0*exp(-r*r/gamma/gamma)
      return
      end
c--------------------------------------------------------------- 
      real*8 function g3(rmag,t2,r2)
      implicit real*8(a-h,o-z)
      g3=rmag*exp(-r2/t2)
      return
      end
c--------------------------------------------------------------- 
      subroutine print1(icc,aa,zz,nuc,E)
      implicit real*8(a-h,o-z)
      zn=1.0
      if(nuc.eq.1) zn=0.0
      write(icc,'(f4.1,a,2f5.1)') zn,'  1.0',zz,aa
      write(icc,'(f5.1)') E
      write(icc,1) '0.005  16000                                      '
      write(icc,1) '0.  30.  0.1                                      '   
      write(icc,1) '1.25                                              '
    1 format(a,a,a)
      return
      end
c--------------------------------------------------------------- 
      subroutine print2(icc)
      implicit real*8(a-h,o-z)
      write(icc,1) 'read                                              '  
      write(icc,1) '  300       0.100    ((5e14.7))                   '
    1 format(a,a,a)
      return
      end
c--------------------------------------------------------------- 
      subroutine finder(aa,rms,rho0,rad0,diff)
      implicit real*8 (a-h,o-z)
      common /dcom/ v0, r0, a0, beta2, b2bar
      pi=4.d0*datan(1.d0) 
c     Woods-Saxon for quadrupole deformed nucleus (beta4 = 0)
           a0=diff
           beta2=0.d0
           rwant=rms
           v0 = 0.13269d0
           r0 = 1.1d0*aa**(1.d0/3.d0) 
           b2bar = r0 * sqrt(5.d0/(4.d0*pi)) * beta2
           call volrms(vol,rms)
           step = 0.05d0
           r0 = r0 + step
           err = rwant -rms
           icount =0
 10           call volrms(vol,rms2)
              err2 = rwant -rms2 
              if (abs(err2).gt.1.d-5) then
                  icount = icount+1
                  write(6,*) 'radius',icount, rms2
                  rms =  rms2
                  grad = (err2-err)/step
                  err = err2
                  step = - err/grad
                  r0 = r0 + step
                  go to 10
              end if
           call volrms(vol,rms)
           step = 0.05d0
           v0 = v0 + step
           err = aa - vol
           icount =0
 20           call volrms(vol2,rms)
              err2 = aa - vol2
              if (abs(err2).gt.1.d-8) then
                  icount = icount+1
                  write(6,*) 'depth',icount, vol2
                  vol = vol2 
                  grad = (err2-err)/step
                  err = err2
                  step = - err/grad
                  v0 = v0 + step
                  go to 20
              end if
           call volrms(vol,rms)
           print*
           print 24, vol, rms
 24        format('Volume: ',f12.7,' rms ',f12.7)
           print 25, v0, r0, a0, beta2
 25        format('depth: ',f12.7,' radius ',f12.7,'\n',
     +         'diffuse: ',f12.7,' defm ',f12.7)
      print*
      rho0=v0
      rad0=r0
      return
      end 
c--------------------------------------------------------------- 
      real*8 function den0(r)
      implicit real*8 (a-h,o-z)
      common /dcom/ v0, r0, a0, beta2, b2bar
           e0 = 1.d0/(1.d0+exp((r-r0)/a0))
           e1 = e0*(1.d0-e0)/a0
           e2 = e1*(1.d0 - 2.d0*e0)/a0
           e3 = e1*(1.d0-6.d0*e0+6.d0*e0*e0)/(a0*a0)
            den0 = v0*( e0 + 
     +                 b2bar*b2bar*e2/10.d0 +
     +           b2bar*b2bar*b2bar*e3/105.d0 )
           return
      end 
c-------------------------------------------------------------
	 subroutine volrms(vol,rms)
	 implicit real*8 (a-h,o-z)
         pi=4.d0*datan(1.d0)
         h = 0.1d0
         rmax = 20.d0
         nval = int(rmax/h) + 1
         if (mod(nval,3).ne.0) nval = nval+3 - mod(nval,3)
         sum = 0.d0 
         sumr2 = 0.d0 
         do i=0, nval
             r = i*h
             r2 = r*r
             v = den0(r)
             r2vs = r2*v*simpfac(i,nval)
             sum = sum + r2vs
             sumr2 = sumr2 + r2*r2vs
          end do  
          sum = 4.d0 * pi * sum * (3.d0*h/8.d0)
         sumr2 = 4.d0 * pi * sumr2 * (3.d0*h/8.d0)
         vol = sum
         rms = sqrt(sumr2/sum)
         return
      end 
c-----------------------------------------------------------
      real*8 function simpfac(i,nstep) 
c     convert degrees to radians
         implicit none 
         integer i,nstep
	 real*8 sfact 
         if ((i.eq.0).or.(i.eq.nstep)) then  
		sfact= 1.d0 
	 else if (mod(i,3).eq.0) then
                sfact= 2.d0
 	 else
                sfact=3.d0
         end if
         simpfac =sfact
         return
      end 
      subroutine gauss(a,b,npoint,xri,wri)
      implicit real*8(a-h,o-z)
      real*8 xg(200),wg(200),xri(200),wri(200)
      call setmgl(npoint,xg,wg)
      do 20 j=1,npoint
      xri(j) = (a+b)/2.d0 + (b-a)/2.d0*xg(j)
      wri(j) = (b-a)/2.d0*wg(j)
   20 continue
      return
      end
      subroutine setmgl( n, points, weight )
      implicit real*8  ( a-h, o-z )
      real*8 points(200), weight(200)
      real*8 poin16(300)
      pi = 4.d0*atan(1.d0)
      if (n.gt.200)  write (1,50)
50    format(' setmlg call with too many points')
      m = ( n + 1 ) / 2
      e1 = n * ( n + 1 )
      do 1 i = 1, m
      t = ( 4*i - 1 ) * pi / ( 4*n + 2 )
      x0 = ( 1.d0 - ( 1.d0 - 1.d0/n ) / ( 8.d0*n*n ) ) * cos(t)
      pkm1 = 1.d0
      pk = x0
      do 3 k = 2, n
      t1 = x0 * pk
      pkp1 = t1 - pkm1 - ( t1-pkm1 )/k + t1
      pkm1 = pk
      pk = pkp1
3     continue
      den = 1.d0 - x0*x0
      d1 = n * ( pkm1 - x0*pk )
      dpn = d1 / den
      d2pn = ( 2.d0*x0*dpn - e1*pk ) / den
      d3pn = ( 4.d0*x0*d2pn + (2.d0-e1)*dpn ) / den
      d4pn = ( 6.d0*x0*d3pn + (6.d0-e1)*d2pn ) / den
      u = pk / dpn
      v = d2pn / dpn
      h = -u * ( 1.d0 + 0.5d0*u*(v+u*(v*v-u*d3pn/(3.d0*dpn))))
      p = pk + h*(dpn+0.5d0*h*(d2pn+h/3.d0*(d3pn+0.25d0*h*d4pn)))
      dp = dpn + h*(d2pn+0.5d0*h*(d3pn+h*d4pn/3.d0))
      h = h - p / dp
      poin16(i) = x0 + h
      fx = d1 - h*e1*(pk+0.5d0*h*(dpn+h/3.d0* 
     1     (d2pn+0.25d0*h*(d3pn+0.2d0*h*d4pn))))
      weight(i) = 2.d0 * ( 1.d0 - poin16(i)*poin16(i)) / (fx*fx)
1     continue
      if ( m + m .gt. n ) poin16(m) = 0.d0
      do 10 i = n/2 + 1, n
      poin16(i) = poin16( n + 1 - i )
      weight(i) = weight( n + 1 - i )
      poin16( n + 1 - i ) = -poin16( n + 1 - i )
10    continue
      do 30 i=1,n
 30   points(i)=poin16(i)
      return
      end
c--------------------------------------------------------------- 
      subroutine finder2(aa,rms,rho0,rad0,diff)
      implicit real*8 (a-h,o-z)
      pi=4.d0*datan(1.d0) 
      rho0=1.d0
      call volosc(vol,rms,rho0,rad0,diff)
      rho0=aa/vol
      call volosc(vol,rms,rho0,rad0,diff)
      print 24, vol, rms
 24   format('  Volume: ',f12.7,' rms ',f12.7)
      return
      end 
c-------------------------------------------------------------
      subroutine volosc(vol,rms,rho0,rad0,diff)
      implicit real*8 (a-h,o-z)
      pi=4.d0*datan(1.d0)
      h = 0.1d0
      rmax = 20.d0
      nval = int(rmax/h) + 1
      if (mod(nval,3).ne.0) nval = nval+3 - mod(nval,3)
      sum = 0.d0 
      sumr2 = 0.d0 
      do i=0, nval
       r = i*h
       r2 = r*r
       v=oscden(r,rho0,rad0,diff)
       r2vs = r2*v*simpfac(i,nval)
       sum = sum + r2vs
       sumr2 = sumr2 + r2*r2vs
      end do  
      sum = 4.d0 * pi * sum * (3.d0*h/8.d0)
      sumr2 = 4.d0 * pi * sumr2 * (3.d0*h/8.d0)
      vol = sum
      rms = sqrt(sumr2/sum)
      return
      end 
c-------------------------------------------------------------
      subroutine vol3pf(vol,rms,rho0,rad0,diff,www)
      implicit real*8 (a-h,o-z)
      pi=4.d0*datan(1.d0)
      h = 0.1d0
      rmax = 20.d0
      nval = int(rmax/h) + 1
      if (mod(nval,3).ne.0) nval = nval+3 - mod(nval,3)
      sum = 0.d0 
      sumr2 = 0.d0 
      do i=0, nval
       r = i*h
       r2 = r*r
       v=fermi3(r,rho0,rad0,diff,www)
       r2vs = r2*v*simpfac(i,nval)
       sum = sum + r2vs
       sumr2 = sumr2 + r2*r2vs
      end do  
      sum = 4.d0 * pi * sum * (3.d0*h/8.d0)
      sumr2 = 4.d0 * pi * sumr2 * (3.d0*h/8.d0)
      vol = sum
      rms = sqrt(sumr2/sum)
      return
      end 
c-------------------------------------------------------------
      real*8 function terp(r,fun,rgrid,npts,ndim)
c------------------------------------------------------------------------------
c     this function calculates, by interpolation, the value of a real
c     function at an arbitrary point r, when the value of the function 
c     (stored in array fun) is known on a grid of points. The values of the
c     npts points at which the function is known is stored in array rgrid.
c     ndim is the externally defined dimensions of the arrays fun and rgrid
c     JAT routine of some vintage.
c------------------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      real*8 fun(ndim),y1,y2,y3,y4,y5,y6
      double precision rgrid(ndim)
      do 30 k=1,npts
      nst=0
      if(rgrid(k).lt.r) goto 30
      nst=max0(k-3,1)
      goto 33
   30 continue
   33 if(nst.gt.npts-5) nst=npts-5
      x1=rgrid(nst+0)
      x2=rgrid(nst+1)
      x3=rgrid(nst+2)
      x4=rgrid(nst+3)
      x5=rgrid(nst+4)
      x6=rgrid(nst+5)
      y1=fun(nst+0)
      y2=fun(nst+1)
      y3=fun(nst+2)
      y4=fun(nst+3)
      y5=fun(nst+4)
      y6=fun(nst+5)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)*(x3-x6)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)*(x4-x6)
      pii5=(x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)*(x5-x6)
      pii6=(x6-x1)*(x6-x2)*(x6-x3)*(x6-x4)*(x6-x5)
  777 xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      xd5=r-x5
      xd6=r-x6
      pi1= xd2*xd3*xd4*xd5*xd6
      pi2= xd1*xd3*xd4*xd5*xd6
      pi3= xd1*xd2*xd4*xd5*xd6
      pi4= xd1*xd2*xd3*xd5*xd6
      pi5= xd1*xd2*xd3*xd4*xd6
      pi6= xd1*xd2*xd3*xd4*xd5
      terp=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4+
     + y5*pi5/pii5+y6*pi6/pii6
      return
      end


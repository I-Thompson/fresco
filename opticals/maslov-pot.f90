!***********************************************************************
!     
!    Copyright (c) 2012, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Ian Thompson, I-Thompson@llnl.gov
!     
!    LLNL-CODE-XXXXX All rights reserved.
!
!    Copyright 2012, I.J. Thompson
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
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
!     
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!     
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!        
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!        
!***********************************************************************
!######### OPTICAL MODEL PARAMETERS for MASLOV03 ##########

!          neutron  on 238U 

! Energy    V     rv    av      W     rw    aw      Vd   rvd   avd      Wd   rwd   awd     Vso   rvso  avso    Wso   rwso  awso  rc

!f7.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f5.3f
	character*4 NAME
	character*8 POTL
	character*100 fname
	real elevels(20), jlevels(20)
	
	Z = 92
	A = 238.052163
	NAME = '238U'
	POTL = 'MASLOV03'
	kp = 1
	nex = 4
	elevels(1:10) =(/ 0., .044916, .14838, 0.30718, 0.5181, 0.7759, 1.0767, 1.4155, 1.7884, 2.1911/)
	jlevels(1:10) =(/ 0., 2., 4., 6., 8., 10., 12., 14., 16., 18. /)
	
	EMIN = 2.5
	EMAX = 2.5
	NE = 1
	
	write(6,1) POTL,NAME
1	format('####### OPTICAL PARAMETERS for ',A8,' ########'//'      neutron on ',a4,// &
     & ' Energy    V     rv    av      W     rw    aw     ', & 
     &           ' Vd   rvd   avd      Wd   rwd   awd     ', &
     &           'Vso   rvso  avso    Wso   rwso  awso  rc')
	
	DE = 1.
	IF(NE>1) DE = (EMAX-EMIN)/(NE-1)
	AC = real(nint(A))**(1./3.)

	DO IE=1,NE
	E = EMIN + (IE-1)*DE
	
	VR = 45.93 - 0.28*E + 0.000573*E*E
	RR = 1.26
	AR = 0.63
	VD = 0.0; RVD=0.; AVD = 0.0
	if(E<8.) then
		WD = 3.14 + 0.436*E
	   else
	    	WD = 6.628
	   endif
	RD = 1.26
	AD = 0.52
	W = 0.0
	RW = RR
	AW = AR
	
	
	VSO = 6.2
	RSO =1.120
	ASO = 0.47
	WSO = 0.0
	
	BETA2 = 0.195
	BETA4 = 0.078
	
	RC = 0.0
	RVOLV = AC * RR
	RVOLW = AC * RW
	RSURF = AC * RD
	
	
	
	write(6,10) E,VR,RR,AR, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, WSO,RSO,ASO, RC
10	format(f7.3, 6(f8.3,2f6.3),f6.3)
      fname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nex)//'-E0000000.in'
        write(fname(8:8),'(i1)') mod(nint(z),10)
        write(fname(9:9),'(i1)') mod(nint(a),10)
        write(fname(24:30),'(f7.3)') e
        write(fname(24:26),'(i3.3)') int(e)
        write(0,*) 'Create file <'//trim(fname)//'>'
        open(1,form='formatted',file=trim(fname))
	write(1,'(a)') 'n+'//trim(NAME)//' with '//trim(POTL)//', s='//CHAR(ICHAR('0')+nex)//' at E ='//fname(24:30)
	write(1,'(a)') 'NAMELIST'
	write(1,'(a)') ' &Fresco  hcm= 0.1 rmatch=  20.000'
	write(1,'(a)') '    jtmin=   0.0 jtmax= 20 absend= 0.000001 '
	write(1,14) nex
14	format('    thmin=0.0 thinc=2 thmax=180. iblock=',i3)
	write(1,'(a)') '    chans= 1 smats= 2 xstabl= 1'
	write(1,15) E
15	format('    elab=',f10.5,'  pel=1 exl=1 lab=1 lin=1 lex=1 /')
	write(1,*) 
 	write(1,16) nex
16	format('&Partition namep=''n       '' massp=  1.008665 zp=  0 nex=',i3)
	write(1,17) NAME,A,Z
17	format('            namet=''',a8,''' masst=',f10.6,' zt=',f5.1,' qval=  0.000/')
 	write(1,'(a)') '&States jp= 0.5 ptyp= 1 ep=  0.000000  cpot=  1 jt= 0.0 ptyt= 1 et= 0.000000/'
 	do 20 il = 2,nex
20 	write(1,21) kp,jlevels(il),elevels(il)
21	format('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt= 1 et=',f8.4,'/')
 	write(1,'(a)') '&Partition /'
	write(1,*) 

	write(1,30) kp,0,0,real(nint(A)),0.0,RC
30	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f9.4,'/')
31	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f9.4,'/')
32	format('&POT /'/)
	write(1,30) kp,1,0,VR,RR,AR
	DELTA2 = BETA2 * RVOLV
	DELTA4 = BETA4 * RVOLV
	write(1,31) kp,11,13, 0., DELTA2, 0., DELTA4, 0., 0.
	if(abs(W)>1e-10) then
	write(1,31) kp,1,0,0.,0.,0., W,RW,AW
	DELTA2 = BETA2 * RVOLW
	DELTA4 = BETA4 * RVOLW
	write(1,31) kp,11,13, 0., DELTA2, 0., DELTA4, 0., 0.
	endif
	write(1,31) kp,2,0,VD,RVD,AVD, WD,RD,AD
	DELTA2 = BETA2 * RSURF  ! assume VD=0
	DELTA4 = BETA4 * RSURF
	write(1,31) kp,11,13, 0., DELTA2, 0., DELTA4, 0., 0.
	write(1,31) kp,3,0,VSO,RSO,ASO, WSO,RSO,ASO
	write(1,32)
	
	write(1,*) '&Overlap /'
	write(1,*) '&Coupling /'

	close(1)
	enddo
	end

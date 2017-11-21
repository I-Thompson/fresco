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
!    along with FRESCO. If not, see <http://www.gnu.org/licenses/>.
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
c  Simple program to read FRESCO bound state wave functions
c  output when NFL<0.
c
c  If job.wf is the output file, use this program as
c	grwf < job.wf > job.plot
c	xmgr job.plot


      parameter(MAXN=10000)
      COMPLEX D(MAXN),VERT(MAXN)
      character*100 COMMENT
C
1        READ(5,'(a)',end=999)  COMMENT
         READ(5,*) NPOINTS,DR,RFIRST
         WRITE(0,71) COMMENT,NPOINTS,DR,RFIRST
71      FORMAT(/'  Input:',a80,/'  Reading ',I4,' data points at '
     X   ,F8.4,' intervals, starting from r =',F8.4,' COMPLEX')
        if(NPOINTS.gt.MAXN) then
          write(0,*) ' Only room to read in ',MAXN,' potential',
     x     ' points, not ',NPOINTS
          stop
          endif
        READ(5,*) (D(I),I=1,NPOINTS)
        READ(5,*) (VERT(I),I=1,NPOINTS)
C
      totd = 0.0
      rms = 0.0
      do 15 i=1,NPOINTS
        r = (i-1)*DR
	d(i) = d(i) * r
        de =d(i)**2  * dr
        totd = totd + de
15      rms = rms + de * r*r 
      rms = sqrt(rms/totd)
      write(0,16) totd,rms
16    format(' Wfntn volume integral =',F8.4,', rms radius =',F8.3)

      DO 200 I=1,NPOINTS
      R = (I-1)*DR

      WRITE(6,190) real(R),D(I),VERT(I)
190 	format(f8.3,4f12.5)
200   CONTINUE
      write(6,*) '&'
	go to 1
C
999   STOP
      END

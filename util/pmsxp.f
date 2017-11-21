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
C
       dimension XSEC(5,5)
      a = 1/sqrt(5.)
      b = sqrt(2./7.)
      c = sqrt(18./35.)
      d = 1/sqrt(14.)
      e = sqrt(8./35.)
      f = sqrt(1./70.)
        read(5,*) KQ1PR
C       write(6,*) KQ1PR
1     read(5,*,err=20,end=20) th,xs,
     X     ((XSEC(KQ1,LQ1),LQ1=2-MOD(KQ1,2),KQ1),KQ1=2,KQ1PR)
C     write(6,*) th,xs,
C    X     ((XSEC(KQ1,LQ1),LQ1=2-MOD(KQ1,2),KQ1),KQ1=2,KQ1PR)
       t20 = XSEC(3,1)
       t40 = XSEC(5,1)
       pm0 = a * ( a - b * t20 + c * t40)
       pm1 = a * ( a - d * t20 - e * t40) * 2.0
       pm2 = a * ( a + b * t20 + f * t40) * 2.0
      print *,th,pm0,pm1,pm2
      go to 1
  
20    STOP
      end

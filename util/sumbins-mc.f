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

* compile with
* ifc7 -Vaxlib -o sumbins-mc{,.f}
* then can pass aguments to program
* angle energy
* program to sum up angular breakup cross section for specific
* sets of energy bins for multichannel bins
      program sumbins2
      implicit real*8(a-h,o-z)
      parameter (ibm=500,mset=3000,nth=1001)
      
      real*8 enlab,beproj,enex(ibm),ysum(nth),ycoresum(mset,nth)
      real*8 jexp,enexp0,enext0,sigt,jpi,jold,jex,jext
      real*8 e(mset,mset),sige(ibm),sigr(ibm),enexp(ibm),enext(ibm)
      real*8 sig(mset,mset),theta(nth),sigch(ibm,nth)
      real*8 siga1(nth),sigi(nth),spin(mset)
      real*8 yjpisum(-100:100,nth),sigjpi(-100:100),sigli(ibm)
      real*8 sigjpic(-100:100,ibm)
      real*8 sigia(ibm),sigiai(ibm,nth),jexc(ibm)
      integer pel,exl,nchan,bandp,bandt
      integer ne(mset),in(ibm),setjpi(ibm),set(ibm,3)
      integer iaa(ibm),bandc(ibm),lib(ibm)
      character line*80,tmp,argv*100
      character*1 sign(-1:1)
      logical incl(mset)

      data sign/'-','?','+'/

!      call getarg(1,argv)
!      read(argv,*,end=10,err=10)angmax
!      call getarg(2,argv)
!      read(argv,*,end=11,err=11)enmax
!      goto 12
!10    angmax=0.
!      goto 12
!11    enmax=0.
!12    continue
      angmax=0.
      enmax=0.

      if(angmax>eps)write(0,'(a,f5.1)')'c.m. angular cutoff = ',angmax
      if(enmax>eps) write(0,'(a,f5.1)')'energy cutoff = ',enmax

      set(1:ibm,1:3)=-99
      
      read(13,*)line
* read in reaction info
      read(13,'(3i3,1p,e12.4,0p,f8.4)') pel,exl,nchan,enlab,beproj
      be=beproj
c      print*,'using binding energy=',be     

      nbex=0
      do ic=1,1
* read in number of channels (including gs elastic=1)
         read(13,'(2i4)') i,nex
* read in ground state cross section (ib=1)
         read(13,'(2(f5.1,i3,f8.4),1p,e12.4)') 
     x        jexp,bandp,enexp(1),jext,bandt,enext(1),sigr(1)
         enold=0
         nset=1
         nen=0
         ne(nset)=1
         in(1)=1
         nin=1
         spin(1)=bandp*jexp
         setjpi(1)=1
         set(1,1)=nset
         set(1,2)=nint(2*jexp)
         set(1,3)=in(1)
         if(mod(nint(2*jexp),2)==1)jsn=1
         if(mod(nint(2*jexp),2)==0)jsn=2
* loop over all exciation channels ib = channel number
         do ib=2,nex         
c     excited states
           read(13,'(2(f5.1,i3,f8.4),1p,e12.4)') 
     x           jexp,bandp,enexp(ib),jext,bandt,enext(ib),sigr(ib)
!           write(0,'(2(f5.1,i3,f8.4),1p,e12.4)') 
!     x           jexp,bandp,enexp(ib),jext,bandt,enext(ib),sigr(ib)

           if (enexp(ib).lt.0.0001) then
             enex(ib)=enext(ib)
             jex=jext
             jpi=jex*bandt
           else
             enex(ib)=enexp(ib)
             jex=jexp
             jpi=jex*bandp
           endif
            
           if (ib==2)jold=jpi
* check if exciation channel is bound state
           if(enex(ib)<be)then
*              bound excited state
             nbex=nbex+1
             nset=nset-1
           endif
           
* figure out how many channels in each jpi set
* iset = jpi channel ref   nset number of jpi combinations
* ne(nset) = number of channels in jpi set
* e(ne(nset),nset) = energy of channel
* sig(ne(nset),nset) = total cross section for channel
* if jpi changes then start new nset and store ne(nset) and reset
            if ((enex(ib).lt.enold).or.(jold.ne.jpi)) then
               ne(nset)=ne(nset)-1
               nset=nset+1
               nen=max(ne(nset),nen)
               ne(nset)=1
c               print*,'new set'
            endif
            if(jold.ne.jpi)then
              in(ib)=1
              setjpi(nset)=setjpi(nset-1)+1
            else
              in(ib)=in(ib-1)
            endif
            if(jold.eq.jpi.and.(enex(ib).le.enold))then
              in(ib)=in(ib)+1
            endif
            nin=max(in(ib),nin)
            setjpi(ib)=nset
            spin(nset)=jpi
         set(ib,1)=nset
         set(ib,2)=nint(2*jpi)
         set(ib,3)=in(ib)

         jpimax=max(nint(2*jpi),jpimax)
         jpimin=min(nint(2*jpi),jpimin)

* store energy and total cross section for channel
            e(ne(nset),nset)=enex(ib)-be
            sig(ne(nset),nset)=sigr(ib)
* count on for next channel
!          if(ib==2)write(0,'(a)')'ib  jpi spin in sjp st ne   e     sig'
!            write(0,'(i3,1x,2f4.1,4i3,2f7.2)') 
!     X       ib,jpi,spin(nset),in(ib),setjpi(nset),nset,ne(nset),
!     X       e(ne(nset),nset),sig(ne(nset),nset)
            ne(nset)=ne(nset)+1
            enold=enex(ib)
            jold=jpi
            inold=in(ib)
         enddo
         ne(nset)=ne(nset)-1
         nen=max(ne(nset),nen)
c         print*,nen
      enddo ! ic
      call get_state(nex,iaa,lib)
      do ic=2,2
* read in number of core states
         read(13,'(2i4)') i,nexc
         write(0,'(a,i3,a)')' reading in ',nexc,' core states'
c     excited states
         do ia=1,nexc
!           read(13,'(2(f5.1,i3,f8.4),1p,e12.4)',end=101) 
!     x           jexc(ia),bandc(ia),e1,jext,bandt,enext(ia),sig
!101        write(0,'(2(f5.1,i3,f8.4),1p,e12.4)') 
!     x           jexc(ia),bandc(ia),e1,jext,bandt,enext(ia),sig
           read(13,'(f5.1,i3,f8.4)') 
     x           jexc(ia),bandc(ia)
!           write(0,'(f5.1,i3,f8.4)') 
!     x           jexc(ia),bandc(ia)
         enddo
      enddo

      sigbu=0.
      sigia(:)=0.
      sigli(:)=0.
      sigjpi(:)=0.
      sigjpic(:,:)=0.
      do ib=2+nbex,nex
        ia=iaa(ib)
        sigia(ia)=sigia(ia)+sigr(ib)
        sigli(set(ib,1))=sigli(set(ib,1))+sigr(ib)
        sigjpi(set(ib,2))=sigjpi(set(ib,2))+sigr(ib)
        sigjpic(set(ib,2),ia)=sigjpic(set(ib,2),ia)+sigr(ib)
        sigbu=sigbu+sigr(ib)
      enddo

      write(0,99)sigbu
99    format('Total breakup cross section = ',f10.3,' mb')
      do ia=1,nexc
        write(0,100)jexc(ia),sign(bandc(ia)),sigia(ia)
      enddo
100   format('Total breakup cross section for core state ',
     &       f4.1,a1,' = ',f10.3,' mb')

      do jpi2=jpimin,jpimax,jsn
        if(sigjpi(jpi2)>0d0)write(0,102)jpi2/2d0,sigjpi(jpi2)
      enddo
102   format('Total breakup cross section for jpi ',
     &       f4.1,' = ',f10.3,' mb')
      
      iset=0
      isetold=0
      do ib=2,nex
        iset=set(ib,1)
        if(iset==isetold)cycle
        ia=iaa(ib)
        write(0,104)
     &          lib(iset),jexc(ia),sign(bandc(ia)),sigli(iset)
        isetold=iset
      enddo
104   format('Total breakup cross section for l I ',
     &       f4.1,1x,f4.1,a1,' = ',f10.3,' mb')


      read(16,'(a2,i3)')tmp,na
      
      read(16,'(//////////)')
      do i=1,na
        read(16,*)theta(i),sigch(1,i)
      enddo
      dtheta=theta(3)-theta(2)

c      print'(a,i3)','number of angles = ',na
c      print'(a,f8.4,a)','in steps of ',dtheta,' deg'
c      print'(a,f8.4,a)','max angle = ',theta(na),' deg'
* now read in angular cross section for each excitation channel
      do ich=2,nex
         read(16,'(////)')
         do i=1,na
           read(16,*) theta(i),sigch(ich,i)
         enddo
      enddo
      
!      write(0,*)'Enter max energy of bins to include or 0 for all'
!      read*,emax
      emax=0.
      if(enmax>eps)then
        emax=enmax
      endif
* now work out which bins to include in sum
      incl(1:nbex+1)=.false.
      incl(nbex+2:nex)=.true.
      if(emax>0.)then
        do ib=2+nbex,nex
          if(enex(ib)-be>emax)incl(ib)=.false.
!          write(0,*)'e,incl ',ich,enex(ich),enex(ich)-be,incl(ich)
        enddo
      endif
      
      write(0,'(a,i4,a,i4,a)')
     X      'Read in ',nex,' sets, each of ',na,' values'

      if(angmax>eps)then
        na1=nint(angmax/dtheta)+1
        write(0,'(a,i4,a)')'  using upto ',na1,' values for ang cutoff'
      endif
!      write(0,'(a,500i3)') 'Inwaves:',in(2+nbex:nex)
!      write(0,'(a,500i3)') 'Jpi set:',setjpi(2+nbex:nex)

      sigia(:)=0.
      ysum(:)=0.
      ycoresum(:,:)=0.
      yjpisum(:,:)=0.
      pi=4d0*atan(1d0)
      deg=180d0/pi
      do ib=2+nbex,nex
        do i=1,na
          thrad=(i-1)*dtheta/deg
          ysum(i)=ysum(i)+sigch(ib,i)
          ycoresum(iaa(ib),i)=ycoresum(iaa(ib),i)+sigch(ib,i)
          yjpisum(set(ib,2),i)=yjpisum(set(ib,2),i)+sigch(ib,i)
!          if(i<=na1 .and. incl(ib))
!     &    sigiai(iaa(ib),i)=sigiai(iaa(ib),i)+sin(thrad)*sigch(ib,i)
        enddo
      enddo

!      do ia=1,nexc
!        call sim(sigiai(ia,:),sigia(ia),1,na1,dtheta/deg)
!        sigia(ia)=2d0*pi*sigia(ia)
!      enddo

      do ia=1,nexc
        if(angmax>eps.or.emax>eps)write(0,101)
     &                  angmax,emax,jexc(ia),sign(bandc(ia)),sigia(ia)
      enddo
101   format('Breakup cross section (theta_max =',f5.1,'  Emax = ',f5.1,
     &')  ',f4.1,a1,' = ',f10.3,' mb')

      write(6,'(a,$)')'#Theta      Total breakup xsec'
      do iset=1,nset
       write(6,'(a,f5.1,4x,$)')'      j/pi=',spin(iset)
      enddo
      write(6,'(a)')''
      do i=1,na
        write(6,'(f10.6,500e14.4)')
     X  theta(i),ysum(i),yjpisum(jpimin:jpimax:jsn,i)
      enddo
      write(6,*)'&'
      write(6,'(a)')'#Theta    corestates'
      do i=1,na
        write(6,'(f10.6,10e14.4)')
     X  theta(i),ycoresum(1:nexc,i)
      enddo
      end
      
      subroutine get_state(n,ia,l)
      implicit none
      integer:: ib,j,n,ia(n),ib1,l(n)
      real:: jp
      character:: line*10
      ia(:)=0
      j=0
      l=0
      ib=1
      do while(j<2)
       read(1,'(a)')line
       if(line(1:10)=='          ')j=j+1
      enddo
      do while(ib<=n)
       if(ib<100)then
        read(1,'(18x,i2,6x,i2,f4.1,i2)')l(ib),ia(ib),jp,ib1
        if(ib1/=ib.and.ib1/=0)write(0,*)'ERROR: IB MISMATCH',ib,ib1
       else
        read(1,'(18x,i2,6x,i2,f4.1)')l(ib),ia(ib),jp
       endif
       if(ia(ib)==0)ia(ib)=1
       ib=ib+1
      enddo
      end
      
*     ---------------------------------------------------------------
      subroutine sim(fa,res,m,n,h)
      implicit real*8(a-h,o-z)
      dimension fa(*),dq(1001)
      do 90 i=m,n
      dq(i)=fa(i)
   90 continue
      rq1=dq(m+1)
      rq2=dq(m+2)
      i=m+3
   98 continue
      if(i.ge.n) go to 99
      rq1=rq1+dq(i)
      rq2=rq2+dq(i+1)
      i=i+2
      go to 98
   99 continue
      res=1d0/3d0*h*(dq(m)+4.d0*rq1+2.d0*rq2-dq(n))
      return
      end

      

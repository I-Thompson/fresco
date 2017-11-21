        program angular_correlation

c see Makefile for compilation info   

      implicit real*8 (a-h,o-z)
      character*30 aname
c
c name of fresco amplitude file
c
      aname='f18dn-zr-l.amp'
c nen=number of energies in amplitude file
      nen=1
c nex=number of final-state excitation energies in amplitude file
      nex=2
c loadf: loads amplitude file
      call loadf(aname,nen,nex)
c output differential cross section for iex=1
      call sig_output(1)
c output slice of correlation function for iex=2
      call output(2)

      do iex=1,nex
        print '(1x,a,i5)','iex=',iex
        th=0.0d0
        tha=0.0d0
        phi=0.0d0
        call wcor(th,tha,phi,iex,w,sig)
        print '(1x,a,f15.5)','sig(0)=',sig
        print '(1x,a,f15.5)','w(0,0)=',w
        call search(iex,wmax,wmin)
        print '(1x,a,f15.5)','max w =',wmax
        print '(1x,a,f15.5)','min w =',wmin
      enddo

      stop
      end

      subroutine output(iex)
c
c routine to output a sample slice of the correlation function
c
      implicit real*8 (a-h,o-z)
      open(10,file='ang-cor.out',status='new')
      th=0.0d0
      phi=0.0d0
      do ia=0,180
        tha=float(ia)
        call wcor(th,tha,phi,iex,w,sig)
        write(10,'(1x,f15.5,e15.5)') tha,w
      enddo
      close(10)
      return
      end

      subroutine search(iex,wmax,wmin)
c
c routine to brute-force search parameter space for max and min values for
c correl function.
c iex = index of intermediate state
c
      implicit real*8 (a-h,o-z)
c
      wmax=1.0d0
      wmin=1.0d0
      do ith=0,360,10
        th=0.5d0*float(ith)
        do itha=0,360,10
          tha=0.5d0*float(itha)
          do iphi=0,180,20
            phi=float(iphi)
            call wcor(th,tha,phi,iex,w,sig)
            if (w.gt.wmax) then
              wmax=w
c              print '(1x,a,f12.5,a,3f12.5)','new wmax =',wmax,
c     &          '  at',th,tha,phi
            endif
            if (w.lt.wmin) then
              wmin=w
c              print '(1x,a,f12.5,a,3f12.5)','new wmin =',wmin,
c     &          '  at',th,tha,phi
            endif
          enddo
        enddo
      enddo
      return
      end

      subroutine wcor(th,tha,phi,iex,w,sig)
c
c routine to calculate angular correlation w
c
c th  =  center-of-mass angle of product nucleus c in degrees
c tha =  emission angle of alpha in c.m. system
c phi =  azimuthal angle between nucleus c and alpha
c iex =  index to intermediate state
c w   =  angular correlation correlation
c sig =  cm differential cross section for th
c
c example for 18F[d,n]19Ne(3/2+) -> alpha(0+) + 15O(1/2-)
c
      implicit real*8 (a-h,o-z)
      include 'ang-cor.inc'
      complex*16 tkq(maxk,2*maxk+1)
      complex*16 f(maxamp),qi,ykq,ckq,wsum
      real*8     rk(maxk)
      real*8     pkq(0:maxk,-maxk:maxk)
      logical    debug
      data debug/.false./

      pi=4.0d0*atan(1.0d0)
      qi=dcmplx(0.0d0,1.0d0)
      call interpf(th,iex,f)

      kmax=nint(2.0d0*qjto(iex))-mod(nint(2.0d0*qjto(iex)),2)

      if (kmax.gt.nint(2.0d0*ql(iex))) kmax=2*nint(ql(iex))
      do k=0,kmax,2
        ik=k/2+1
        qk=real(k)
        rk(ik)=sqrt(2.*qjto(iex)+1.)*(2.*ql(iex)+1.0d0)*
     &    (-1)**(k+nint(qic(iex)-qjto(iex)))*
     &    dclebg(ql(iex),ql(iex),qk,0.0d0,0.0d0,0.0d0)*
     &    dracaw(ql(iex),ql(iex),qjto(iex),qjto(iex),qk,qic(iex))
      enddo
c
      if (debug) then
        write(6,'(1x,a,f5.2)') 'I_B =',qjto(iex)
        write(6,'(1x,a,f5.2)') 'S   =',qic(iex)
        write(6,'(1x,a,10f12.5)') 'R_k = ',(rk(ik),ik=1,kmax/2+1)
      endif
c
      do ik=1,kmax/2+1
        do iq=1,2*kmax+1
          tkq(ik,iq)=0.0d0
        enddo
      enddo
      sum=0.0d0
      do mpi=1,njpi(iex)
        qmpi=-qjpi(iex)+float(mpi-1)
        do mti=1,njti(iex)
          qmti=-qjti(iex)+float(mti-1)
          do mpo=1,njpo(iex)
            qmpo=-qjpo(iex)+float(mpo-1)
            do mto=1,njto(iex)
              qmto=-qjto(iex)+float(mto-1)
              ip1=njto(iex)*(njpo(iex)*(njti(iex)*
     &          (mpi-1)+mti-1)+mpo-1)+mto
              sum=sum+real(f(ip1)*dconjg(f(ip1)))
              do k=0,kmax,2
                ik=k/2+1
                qk=float(k)
                do kq=-k,k
                  iq=kq+k+1
                  qq=float(kq)
                  qmq=qmto+qq
                  mtq=nint(qmq+qjto(iex)+1.0d0)
                  ip2=njto(iex)*(njpo(iex)*(njti(iex)*(mpi-1)+
     &              mti-1)+mpo-1)+mtq
                  if ((ip2.ge.1).and.(ip2.le.namp(iex))) then
                    tkq(ik,iq)=tkq(ik,iq)+
     &                dconjg(f(ip1))*
     &                f(ip2)*
     &                sqrt(2.0d0*qk+1.0d0)*
     &                dclebg(qjto(iex),qk,qjto(iex),qmto,qq,qmq)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      do ik=1,kmax/2+1
        do iq=1,2*kmax+1
          tkq(ik,iq)=tkq(ik,iq)/sum
        enddo
      enddo
      sig=sum/float(njpi(iex)*njti(iex))*10.0d0
      if (debug) then
        do k=0,kmax,2
          ik=k/2+1
          do kq=-k,k
            iq=kq+k+1
            print '(2i5,a,2f15.5)',k,kq,': ',tkq(ik,iq)
          enddo
        enddo
      endif

      cta=cos(tha*pi/180.0d0)
      do kq=-kmax,kmax
c get associated Legendres, w/o Condon/Shortley phase
        call daslgf(2,cta,kq,kmax,pkq(0,kq))
      enddo
      phir=phi*pi/180.0d0
      wsum=0.0d0
      do k=0,kmax,2
        ik=k/2+1
        qk=float(k)
        do kq=-k,k
          iq=kq+k+1
          qq=float(kq)
c calculate spherical harmonic, WITH Condon/Shortley phase
          ykq=pkq(k,kq)/sqrt(2.0d0*pi)*exp(qi*qq*phir)
          if (kq.gt.0) ykq=ykq*(-1)**kq
          ckq=sqrt(4.0d0*pi/float(2*k+1))*ykq
          wsum=wsum+tkq(ik,iq)*rk(ik)*ckq
        enddo
      enddo
      if (abs(imag(wsum)).gt.1.0e-10) then
        print '(1x,a)','problematic imag component of w...'
        stop
      endif
      w=real(wsum)
      return
      end

      subroutine sig_output(iex)

      implicit real*8 (a-h,o-z)

      open(10,file='ang-cor.dsdo',status='new')
      do ith=0,360
        theta=0.5d0*float(ith)
        sig=dsdo(theta,iex)
        write(10,'(f8.2,e18.8)') theta,sig
      enddo
      close(10)

      return
      end

      real*8 function dsdo(theta,iex)
c
c Calculate differential cross section from the scattering amplitudes
c
      implicit real*8 (a-h,o-z)
      include 'ang-cor.inc'
      complex*16 f(maxamp)

      call interpf(theta,iex,f)
      sum=0.0d0
      do mpi=1,njpi(iex)
        qmpi=-qjpi(iex)+float(mpi-1)
        do mti=1,njti(iex)
          qmti=-qjti(iex)+float(mti-1)
          do mpo=1,njpo(iex)
            qmpo=-qjpo(iex)+float(mpo-1)
            do mto=1,njto(iex)
              qmto=-qjto(iex)+float(mto-1)
              ip1=njto(iex)*(njpo(iex)*(njti(iex)*
     &          (mpi-1)+mti-1)+mpo-1)+mto
              sum=sum+real(f(ip1)*dconjg(f(ip1)))
            enddo
          enddo
        enddo
      enddo
      dsdo=sum/float(njpi(iex)*njti(iex))*10.0d0

      return
      end

      subroutine loadf(aname,nen,nx)
c
c Routine to read fresco scattering amplitude file.
c Note that this file includes all of the spin information.
c Note that only the last energy (ien=nen) is
c actually saved for later calculation.
c Finally, note that 2 extra angles are added on the front and back to
c allow for easier interpolation.
c
      implicit real*8 (a-h,o-z)
      include 'ang-cor.inc'
      character*30 aname
      logical pinfo
      data pinfo/.true./
c
      nex=nx
      open(10,file=aname,status='old')
      do ien=1,nen
        do iex=1,nex
c
c must set ql (orbital L of breakup) and qic (breakup channel spin)
c by hand for each final state
c
          if (iex.eq.1) then
            ql(iex)=1.0d0
            qic(iex)=0.5d0
          elseif (iex.eq.2) then
            ql(iex)=2.0d0
            qic(iex)=0.5d0
          else
            ql(iex)=-1
            qic(iex)=-1
          endif
          read(10,*,end=400) qjpi(iex),qjti(iex),qjpo(iex),qjto(iex),
     &      nang
          njpi(iex)=nint(2.0d0*qjpi(iex)+1.0d0)
          njti(iex)=nint(2.0d0*qjti(iex)+1.0d0)
          njpo(iex)=nint(2.0d0*qjpo(iex)+1.0d0)
          njto(iex)=nint(2.0d0*qjto(iex)+1.0d0)
          namp(iex)=njpi(iex)*njti(iex)*njpo(iex)*njto(iex)
          if (namp(iex).gt.maxamp) then
            print '(1x,a)','maxamp exceeded'
            stop
          endif
          do iang=3,nang+2
            read(10,*,end=400) ang(iang)
            read(10,'(6e12.4)',end=400) (fam(iex,iang,i),i=1,namp(iex))
          enddo
          ang(1)=-ang(5)
          ang(2)=-ang(4)
          ang(nang+3)=ang(nang+2)+ang(nang+2)-ang(nang+1)
          ang(nang+4)=ang(nang+2)+ang(nang+2)-ang(nang)
          do i=1,namp(iex)
            fam(iex,1,i)=fam(iex,5,i)
            fam(iex,2,i)=fam(iex,4,i)
            fam(iex,nang+3,i)=fam(iex,nang+1,i)
            fam(iex,nang+4,i)=fam(iex,nang,i)
          enddo
          nang=nang+4
          if (pinfo) then
            print '(1x,a,i3)','iex = ',iex
            print '(1x,a,f8.1)','  qjpi:',qjpi(iex)
            print '(1x,a,f8.1)','  qjti:',qjti(iex)
            print '(1x,a,f8.1)','  qjpo:',qjpo(iex)
            print '(1x,a,f8.1)','  qjto:',qjto(iex)
            print '(1x,a,f8.1)','  ql:  ',ql(iex)
            print '(1x,a,f8.1)','  qic: ',qic(iex)
            print '(1x,a,i3)',  '  nang = ',nang
            do i=1,5
              print '(1x,a,i3,a,3f9.3)','  ang(',i,'), fam(iex,i,1):',
     &          ang(i),fam(iex,i,1)
            enddo
            do i=nang-4,nang
              print '(1x,a,i3,a,3f9.3)','  ang(',i,'), fam(iex,i,1):',
     &          ang(i),fam(iex,i,1)
            enddo
          endif
        enddo
      enddo
      goto 500
c
400   continue
      print '(1x,a)','premature end of file encountered'
      stop
c
500   continue
      close(10)
      return
      end

      subroutine interpf(theta,iex,f)
c
c routine to interpolate scattering amplitudes for a specific value
c of theta
c
      implicit real*8 (a-h,o-z)
      include 'ang-cor.inc'
      complex*16 f(maxamp),x,y

c
      IF(theta.GE.ang(nang-1)) THEN
        PRINT '(1X,A)','angle out of range'
        STOP
      ENDIF
      IF(theta.LE.ang(2)) THEN
        print '(1x,a)','angle out of range'
        STOP
      ENDIF
      I=0
10    I=I+1
      IF (theta.GT.ang(i)) GOTO 10
      do j=1,namp(iex)
        x= fam(iex,I-1,j)*(theta   -ang(I)   )*
     &                   ( theta   -ang(I+1)) /
     &                   ((ang(I-1)-ang(I)   )*
     &                   ( ang(I-1)-ang(I+1)))+
     &     fam(iex,I,j)  *(theta   -ang(I-1) )*
     &                   ( theta   -ang(I+1)) /
     &                   (( ang(I) -ang(I-1) )*
     &                   ( ang(I)  -ang(I+1)))+
     &     fam(iex,I+1,j)*(theta   -ang(I-1) )*
     &                   ( theta   -ang(I)  ) /
     &                   ((ang(I+1)-ang(I-1) )*
     &                   ( ang(I+1)-ang(I)  ))
        y= fam(iex,I-2,j)*(theta   -ang(I-1) )*
     &                   ( theta   -ang(I)  ) /
     &                   ((ang(I-2)-ang(I-1) )*
     &                   ( ang(I-2)-ang(I)  ))+
     &     fam(iex,I-1,j)*(theta   -ang(I-2) )*
     &                   ( theta   -ang(I)  ) /
     &                   ((ang(I-1)-ang(I-2) )*
     &                   ( ang(I-1)-ang(I)  ))+
     &     fam(iex,I,j)  *(theta   -ang(I-2) )*
     &                   ( theta   -ang(I-1)) /
     &                   ((ang(I)  -ang(I-2) )*
     &                   ( ang(I)  -ang(I-1)))
        f(j)=0.5d0*(x+y)
      enddo
c
      return
      end


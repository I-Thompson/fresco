
Date: Fri, 16 Nov 2001 12:16:06 GMT
From: Dr Natalia Timofeyuk <N.Timofeyuk@surrey.ac.uk>
To: I.Thompson@surrey.ac.uk


      complex f1,f2
      real ja,jb,jap,jbp
      dimension f1(100,500),f2(100,500)
      open(unit=1,file='amplitudes72')
      open(unit=2,file='amplitudes73')
      read(1,*)ja,jb,jap,jbp,nangl
      read(2,*)ja,jb,jap,jbp,nangl
      ni1=int((2*ja+1)*(2*jb+1)*(2*jap+1)*(2*jbp+1))
c     print *,ni1
      do ith=1,nangl
      read(1,*)angl
      read(2,*)angl
      if(ith.eq.1) angl1=angl
      if(ith.eq.2) angl2=angl
      read(1,90)(f1(i,ith),i=1,ni1)
      read(2,90)(f2(i,ith),i=1,ni1)
90    format( 6e12.5)
      enddo
      dth=angl2-angl1
      do ith=1,nangl
      jth=180./dth-ith+2
      sn1= 1.   
      sn2=-1.   
      cs =0
      cs1=0
      cs2=0
      do i=1,ni1        
      cs=cs+abs(sn1*f1(i,ith)+sn2*f2(i,jth))**2*10/((2*ja+1)*(2*jap+1))
      cs1=cs1+abs(sn1*f1(i,ith))**2*10/((2*ja+1)*(2*jap+1))
      cs2=cs2+abs(sn2*f2(i,jth))**2*10/((2*ja+1)*(2*jap+1))
      enddo
c     print *,ith,jth
      print *,(ith-1)*dth,cs 
c     print 9,(ith-1)*dth,(f1(i,ith),i=1,ni1),cs1
c     print 9,(ith-1)*dth,(f2(i,jth),i=1,ni1),cs2
9     format(f8.4,8d12.4,d15.4)
      enddo
      end

Dr. N. Timofeyuk
Department of Physics
University of Surrey
Guildford
Surrey GU2 7XH 
U.K.

Tel: 44-1483-68.67.95   
Fax: 44-1483-68.67.81
E-mail: N.Timofeyuk@surrey.ac.uk

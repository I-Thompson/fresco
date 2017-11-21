!	write(48,*) ' lanecoup =',lanecoup
	if(lanecoup) then
           IF(LISTCC.GT.1) write(6,*) 'LANE COUPLINGS: '
      DO 1325 IC1=1,MXP
	NEXK = NCHPART(IC1)
	ETATG = (MASS(2,IC1)-2*MASS(4,IC1))/MASS(2,IC1)
	ETAPR = (MASS(1,IC1)-2*MASS(3,IC1))/MASS(1,IC1)
	ETAS = - ETATG * ETAPR
	 T = 2*sqrt(abs(ETAS)/MASS(2,IC1))   ! coefficient for charge-exchange coupling
	I0 = CHBASE(IC1)
	 write(48,*) ' Partition ',IC1,' has chs ',I0+1,I0+NEXK
         DO 1322 C=I0+1,I0+NEXK
         IA = EXCIT(C,1)
         KPA = CPOT(IC1,IA)
	 write(48,*) ' Look for Lane couplings IC,IA =',IC1,IA, 
     x 				' with KP =',KPA
	   do JF=1,NF0
	   if (PTYPE(1,JF)==KPA) then
            JFTT = PTYPE(2,JF)
	   if (JFTT==26.or.JFTT==27) then   ! Vol or Surface Lane form factors
                KK = PTYPE(3,JF)
       		if(KK==7) KK=0
	 	write(48,*) ' Found Lane couplings IC,IA =',IC1,IA, 
     x 				' with KP =',KPA,', T,K =',JFTT,KK
	DO 1320 IC2=1,MXP   
	 IB = IA 		!  assume identical IA values coupled by the Lane coupling
	 KPA2= CPOT(IC2,IA)   
	 if(KPA==KPA2 .and.IC2.ne.IC1) then    ! only off-diagonal couplings here: others already done
	 	write(48,*) ' Match Lane couplings IC2,IB =',IC2,IB, ' cf ',IC1
	  I2 = CHBASE(IC2)
	  do 1312 C2=I2+1,I2+NCHPART(IC2)
	  if(ABS(JPROJ(C2)-JPROJ(C))>.1 .or. ABS(JTARG(C2)-JTARG(C))>.1
     x   .or.LVAL(C2).ne.LVAL(C)) go to 1312 			! no change of projectile or target spins!

          S =  TENSOR(JFTT,LVAL(C),JPROJ(C),JVAL(C),JTARG(C),
     X                JTOTAL,LVAL(C2),JPROJ(C2),JVAL(C2),JTARG(C2),
     X                KK,JEX(3,IC1,IA),JEX(4,IC1,IA),
     X                 JEX(3,IC2,IB),JEX(4,IC2,IB),
     X          MASS(3,IC1),MASS(4,IC1),PTYPE(5,JF),T)
          S = S * CI**NINT( - ABS(JPROJ(C2)-JPROJ(C))
     X                   +     JPROJ(C2)-JPROJ(C)
     X                   - ABS(JTARG(C2)-JTARG(C))
     X                   +     JTARG(C2)-JTARG(C))
C           The above phase factors with JPROJ & JTARG etc.,
C           are there only because of definition of M(Ek) matrix element
        S = S * CI**(LVAL(C2)-LVAL(C))   !* SCALEL

          NC = NCLIST(C,C2)+1
           if(NC>MCLIST) then
                write(KO,*) 'For channels ',C,C2
                write(KO,*) 'Need NC=',NC,'>MCLIST=',MCLIST,' !!!'
                write(KO,*) 'So far use forms ',NFLIST(C,C2,:),' need',IF
                call check(NC,MCLIST,30)
                stop
                endif
           CLIST(C,C2,NC) = S
           NFLIST(C,C2,NC) = JF
           NCLIST(C,C2) =  NC
           IF(C.NE.C2) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1)
     X  WRITE(KO,1330) C,C2,JF,(PTYPE(II,JF),II=1,6),NC,0,0,1.,S

1312 	CONTINUE   ! C2

	 endif     ! KPA==KPA2
1320 	CONTINUE   ! IC2
	endif
	endif
	enddo      ! JF
1322 	CONTINUE   ! C
1325 	CONTINUE   ! IC1
	endif ! lanecoup

      SUBROUTINE E01DAF(MX,MY,X,Y,F,PX,PY,LAMDA,MU,C,WRK,IFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER           IFAIL, MX, MY, PX, PY
      DOUBLE PRECISION  C(MX*MY), F(MX*MY), LAMDA(MX+4), MU(MY+4),
     *                  WRK((MX+6)*(MY+6)), X(MX), Y(MY)

      write(0,*) 'Nag library not available for 2D interpolation. Stop'
	stop
	end

      SUBROUTINE E02DFF(MX,MY,PX,PY,X,Y,LAMDA,MU,C,FF,WRK,LWRK,IWRK,
     *                  LIWRK,IFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER           IFAIL, MX, MY, PX, PY
      DOUBLE PRECISION  C((PX-4)*(PY-4)), FF(MX*MY), LAMDA(PX), MU(PY),
     *                  WRK(LWRK), X(MX), Y(MY)
      INTEGER           IWRK(LIWRK)

      write(0,*) 'Nag library not available for 2D interpolation. Stop'
	stop
	end


SUBROUTINE SMOOTHER(msk,X,lx,ly,Xout)
  !!
  IMPLICIT none
  !!
  INTEGER, INTENT(in) :: lx, ly
  INTEGER, DIMENSION(lx,ly), INTENT(in) :: msk
  REAL(4), DIMENSION(lx,ly), INTENT(in) :: X
  REAL(4), DIMENSION(lx,ly), INTENT(out) :: Xout
  INTEGER i, j
  REAL(4) :: w1, w2, w3, w4, sumw
  !!
  Xout = X
  DO i=2, lx-1
     DO j=2, ly-1
        !!
        w1 = 1. ;   w2 = 1. ;  w3 = 1. ;  w4 = 1.
        !!
        IF ( msk(i,j) == 1 ) THEN
           !!
           IF ( msk(i+1,j) == 0 ) w1 = 0.
           IF ( msk(i,j+1) == 0 ) w2 = 0.
           IF ( msk(i-1,j) == 0 ) w3 = 0.
           IF ( msk(i,j-1) == 0 ) w4 = 0.
           !!
           sumw = w1 + w2 + w3 + w4
           !!
           IF ( sumw /= 0. ) THEN
              Xout(i,j) = 0.5*Xout(i,j)  &
                   & + 0.5/sumw*(w1*Xout(i+1,j) + w2*Xout(i,j+1) + w3*Xout(i-1,j) + w4*Xout(i,j-1))
           END IF
           !!
        END IF
     END DO
  END DO
  !!
END SUBROUTINE SMOOTHER

SUBROUTINE BLENDER(lx, ly, l_min, l_max, vl, X)
  !!
  IMPLICIT none
  !!
  INTEGER, INTENT(in) :: lx, ly
  REAL(4), INTENT(in) :: l_min, l_max
  !!
  REAL(4), DIMENSION(ly),    INTENT(in) :: vl
  REAL(4), DIMENSION(lx,ly), INTENT(inout) :: X
  !!
  INTEGER i, j, j1, j2
  REAL(4) :: slp
  INTEGER, PARAMETER :: lbld = 10
  !!
  !!
  DO j=1, ly - 1
     IF ( (vl(j)    < l_min).and.(vl(j+1) >= l_min) ) j1 = j+1
     IF ( (vl(j+1) >= l_max).and.(vl(j)   <  l_max) ) j2 = j
  END DO
  !!
  !!
  DO i = 1, lx
     slp = (X(i,j1-lbld) - X(i,j1))/(vl(j1-lbld) - vl(j1))
     DO j=j1-1, j1-lbld, -1
        X(i,j) = X(i,j1) + slp*(vl(j) - vl(j1))
     END DO
     !!
     slp = (X(i,j2+lbld) - X(i,j2))/(vl(j2+lbld) - vl(j2))
     DO j=j2+1, j2+lbld
        X(i,j) = X(i,j2) + slp*(vl(j) - vl(j2))
     END DO
  END DO
  !!
END SUBROUTINE BLENDER

SUBROUTINE BLENDER_PERSIST(lx, ly, l_min, l_max, vl, X)
  !!
  IMPLICIT none
  !!
  INTEGER, INTENT(in) :: lx, ly
  REAL(4), INTENT(in) :: l_min, l_max
  !!
  REAL(4), DIMENSION(ly),    INTENT(in) :: vl
  REAL(4), DIMENSION(lx,ly), INTENT(inout) :: X
  !!
  INTEGER i, j, j1, j2
  REAL(4) :: slp
  INTEGER, PARAMETER :: lbld = 10
  !!
  !!
  DO j=1, ly - 1
     IF ( (vl(j)    < l_min).and.(vl(j+1) >= l_min) ) j1 = j+1
     IF ( (vl(j+1) >= l_max).and.(vl(j)   <  l_max) ) j2 = j
  END DO
  !!
  !!
  DO i = 1, lx
     DO j=j1-1, j1-lbld, -1
        X(i,j) = X(i,j1)
     END DO
     !!
     DO j=j2+1, j2+lbld
        X(i,j) = X(i,j2)
     END DO
  END DO
  !!
END SUBROUTINE BLENDER_PERSIST

SUBROUTINE exp_msk(lx, ly, xmsk_in, xmsk_out)
  !!
  !!
  INTEGER, INTENT(in)  :: lx, ly
  !!
  INTEGER, DIMENSION(lx,ly), INTENT(in)   :: xmsk_in
  INTEGER, DIMENSION(lx,ly), INTENT(out)  :: xmsk_out
  !!
  !!
  INTEGER :: ji, jj
  !!
  xmsk_out = xmsk_in
  !!
  DO ji=2, lx-1
     DO jj=2, ly-1
        !!
        !!
        !! Western coasts :
        IF ( (xmsk_in(ji+1,jj) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji-1,jj) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
        !! Eastern coasts :
        IF ( (xmsk_in(ji-1,jj) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji+1,jj) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
        !!
        !! Southern coasts :
        IF ( (xmsk_in(ji,jj+1) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji,jj-1) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
        !! Northern coasts :
        IF ( (xmsk_in(ji,jj-1) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji,jj+1) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
     END DO
  END DO
  !!
  !!
END SUBROUTINE exp_msk


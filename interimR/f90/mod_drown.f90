MODULE MOD_DROWN

  IMPLICIT none

  PRIVATE

  PUBLIC :: drown, fill_extra_bands, extra_2_east, extra_2_west, test_xyz

  LOGICAL, PARAMETER :: ldebug = .FALSE.

CONTAINS





  SUBROUTINE DROWN(k_ew, Xin, mask, Xout,nb_inc, nb_smooth, ni,nj)

    !!#############################################################################
    !!
    !!  PURPOSE : fill continental areas of field X (defined by mask=0)
    !!  -------   using nearest surrounding sea points (defined by mask=1)
    !!            field X is absoluletly unchanged on mask=1 points
    !!
    !!  k_ew :  east-west periodicity on the input file/grid
    !!          k_ew = -1  --> no periodicity
    !!          k_ew >= 0  --> periodicity with overlap of k_ew points
    !!
    !!  X    :  treated array                             (2D array)
    !!  mask :  land-sea mask    INTEGER !!!!             (2D array)
    !!
    !! Optional:
    !!  * nb_inc : how far in terms of number of grid-point we extrapolate sea value into land
    !!                => default: nb_inc = 400
    !!                     (will normally stop before 400 iterations, when all land points have been treated!!!)
    !!
    !!  * nb_smooth : number of times the smoother is applied on masked region (mask=0)
    !!                => default: nb_smooth = 2
    !!
    !!
    !!                       Author : Laurent BRODEAU, 2014
    !!
    !!#############################################################################

    !! Arguments :
    INTEGER,                         INTENT(in)    :: k_ew, ni,nj
    REAL(4),    DIMENSION(ni,nj),    INTENT(in)    :: Xin
    REAL(4),    DIMENSION(ni,nj),    INTENT(out)   :: Xout
    INTEGER(2), DIMENSION(ni,nj),    INTENT(in)    :: mask

    INTEGER,    OPTIONAL,          INTENT(in)    :: nb_inc, nb_smooth

    !! Local :
    INTEGER(2), ALLOCATABLE, DIMENSION(:,:) :: maskv, mask_coast, mtmp
    REAL(4),    ALLOCATABLE, DIMENSION(:,:) :: dold, xtmp, X

    INTEGER :: &
         &      ninc_max,      &
         &      nsmooth_max,   &
         &      jinc,          &
         &      ji, jj, jci,   &
         &      jim, jip, js

    REAL(4), PARAMETER :: rr = 0.707

    INTEGER, DIMENSION(2) :: ivi, vim_per, vip_per

    INTEGER, PARAMETER :: jinc_debg = 2


    ninc_max = 200   ! will stop before when all land points have been treated!!!
    IF ( present(nb_inc) ) ninc_max = nb_inc

    nsmooth_max = 2
    IF ( present(nb_smooth) ) nsmooth_max = nb_smooth


    IF ( (size(Xin,1) /= size(mask,1)).OR.(size(Xin,2) /= size(mask,2)) ) THEN
       PRINT *, 'ERROR, mod_drown.F90 => DROWN : size of data and mask do not match!!!'; STOP
    END IF


    !! Backing up original mask into mask2(:,:)
    ALLOCATE ( maskv(ni,nj), dold(ni,nj), xtmp(ni,nj), mask_coast(ni,nj), mtmp(ni,nj) )
    ALLOCATE ( X(ni,nj) )


    !lolo
    !! First masked region of X should be filled with average value of the field for security...
    !lolo

    X = Xin


    ivi = (/ 1 , ni /)


    IF (k_ew >= 0) THEN
       vim_per = (/ ni-k_ew ,  ni-1  /)
       vip_per = (/    2    , 1+k_ew /)
    END IF

    jinc = 0
    maskv(:,:) = mask(:,:)



    DO jinc = 1, ninc_max

       !! Quiting if no land point left:
       IF ( .NOT. (ANY(maskv == 0))  ) THEN
          IF ( ldebug ) PRINT *, 'DROWN: No land points left! Leaving incursion loop at jinc =', jinc
          EXIT
       END IF

       dold(:,:) = X(:,:)




       !! Building mask of the coast-line (belonging to land points)

       mask_coast(:,:) = 0

       mask_coast(2:ni-1,2:nj-1) = (maskv(3:ni,2:nj-1) + maskv(2:ni-1,3:nj) + maskv(1:ni-2,2:nj-1) + maskv(2:ni-1,1:nj-2)) &
            &                     *(-(maskv(2:ni-1,2:nj-1)-1))


       !! West and East boundaries with periodicity
       !! ------------------------------------------

       IF (k_ew >= 0) THEN

          DO jci = 1, 2
             jim = vim_per(jci)  ! ji-1
             ji  = ivi(jci)      ! first ji = 1, then ji = ni
             jip = vip_per(jci)  ! ji+1
             mask_coast(ji,2:nj-1) = (maskv(jip,2:nj-1) + maskv(ji,3:nj) + maskv(jim,2:nj-1) + maskv(ji,1:nj-2)) &
                  &                     *(-(maskv(ji,2:nj-1)-1))
          END DO

       ELSE

          !! West LBC:
          mask_coast(1,2:nj-1)  = (maskv(2,2:nj-1) + maskv(1,3:nj)     + maskv( 1,1:nj-2))*(-(maskv( 1,2:nj-1) -1))
          !! East LBC:
          mask_coast(ni,2:nj-1) = (maskv(ni,3:nj) + maskv(ni-1,2:nj-1) + maskv(ni,1:nj-2))*(-(maskv(ni,2:nj-1) -1))

       END IF



       ! -------
       ! jj=1
       ! -------
       mask_coast(2:ni-1,1) = (maskv(3:ni,1) + maskv(2:ni-1,2) + maskv(1:ni-2,1)) &
            &                     *(-(maskv(2:ni-1,1)-1))
       !!
       !! ji=1, jj=1
       IF (k_ew >= 0) THEN
          mask_coast(1,1) = (maskv(2,1) + maskv(1,2) + maskv(ni-k_ew,1))*(-(maskv(1,1)-1))
       ELSE
          mask_coast(1,1) = (maskv(2,1) + maskv(1,2)                   )*(-(maskv(1,1)-1))
       END IF
       ! ji=ni, jj=1
       IF (k_ew >= 0) THEN
          mask_coast(ni,1) = (maskv(1+k_ew,1) + maskv(ni,2) + maskv(ni,1))*(-(maskv(ni,1)-1))
       ELSE
          mask_coast(ni,1) = (                  maskv(ni,2) + maskv(ni,1))*(-(maskv(ni,1)-1))
       END IF

       ! jj=nj
       ! -------
       mask_coast(2:ni-1,nj) = (maskv(3:ni,nj) + maskv(1:ni-2,nj) + maskv(2:ni-1,nj-1)) &
            &                     *(-(maskv(2:ni-1,nj)-1))
       !! ji=1, jj=nj
       IF (k_ew >= 0) THEN
          mask_coast(1,nj)  = (maskv(2,nj)            + maskv(ni-k_ew,nj)  + maskv(1,nj-1)   )*(-(maskv(1,nj) -1))
       ELSE
          mask_coast(1,nj)  = (maskv(2,nj)                                 + maskv(1,nj-1)   )*(-(maskv(1,nj) -1))
       END IF
       !! ji=ni, jj=nj
       IF (k_ew >= 0) THEN
          mask_coast(ni,nj) = (maskv(k_ew+1,nj) + maskv(ni-1,nj)     + maskv(ni-1,nj-1))    *(-(maskv(ni,nj) -1))
       ELSE
          mask_coast(ni,nj) = (          maskv(ni-1,nj)     + maskv(ni-1,nj-1))             *(-(maskv(ni,nj) -1))
       END IF


       !! mask_coast is fine now
       mtmp(:,:) = mask_coast(:,:)
       mask_coast(:,:) = 0

       WHERE ( mtmp(:,:) > 0 )
          mask_coast = 1
       END WHERE

       !! mask_coast done, time to fill the coastline points with values from the nearest se points
       !! -----------------------------------------------------------------------------------------

       !! Center of the domain:
       DO ji = 2, ni-1
          DO jj = 2, nj-1
             IF ( mask_coast(ji,jj) == 1 ) THEN
                X(ji,jj) = 1./(maskv(ji+1,jj)+maskv(ji,jj+1)+maskv(ji-1,jj)+maskv(ji,jj-1) + &
                     & rr*(maskv(ji+1,jj+1)+maskv(ji-1,jj+1)+maskv(ji-1,jj-1)+maskv(ji+1,jj-1)))*( &
                     & maskv(ji+1,jj)*dold(ji+1,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                     & maskv(ji-1,jj)*dold(ji-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                     & rr*maskv(ji+1,jj+1)*dold(ji+1,jj+1) + rr*maskv(ji-1,jj+1)*dold(ji-1,jj+1) + &
                     & rr*maskv(ji-1,jj-1)*dold(ji-1,jj-1) + rr*maskv(ji+1,jj-1)*dold(ji+1,jj-1)  )
             END IF
          END DO
       END DO



       DO jci = 1, 2


          ji  = ivi(jci)      ! first ji = 1, then ji = ni


          IF (k_ew >= 0) THEN

             !! West and East boundaries with periodicity


             jim = vim_per(jci)  ! ji-1
             jip = vip_per(jci)  ! ji+1

             DO jj = 2, nj-1
                IF ( mask_coast(ji,jj) == 1 ) THEN
                   X(ji,jj) = 1./(maskv(jip,jj)+maskv(ji,jj+1)+maskv(jim,jj)+maskv(ji,jj-1) + &
                        & rr*(maskv(jip,jj+1)+maskv(jim,jj+1)+maskv(jim,jj-1)+maskv(jip,jj-1)))*( &
                        & maskv(jip,jj)*dold(jip,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & maskv(jim,jj)*dold(jim,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & rr*maskv(jip,jj+1)*dold(jip,jj+1) + rr*maskv(jim,jj+1)*dold(jim,jj+1) + &
                        & rr*maskv(jim,jj-1)*dold(jim,jj-1) + rr*maskv(jip,jj-1)*dold(jip,jj-1)  )
                END IF
             END DO

          ELSE

             !! West & East LBCs when not east-west periodicity, extrapolating lineraly
             IF ( ji == 1 ) THEN
                DO jj = 2, nj-1
                   IF ( mask_coast(ji,jj) == 1 ) THEN
                      X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj+1)+maskv(ji,jj-1) + &
                           & rr*maskv(2,jj+1)+rr*maskv(2,jj-1))*( &
                           & maskv(2,jj)*dold(2,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                           & maskv(ji,jj-1)*dold(ji,jj-1) + &
                           & rr*maskv(2,jj+1)*dold(2,jj+1) + &
                           & rr*maskv(2,jj-1)*dold(2,jj-1)  )
                   END IF
                END DO
             END IF
             IF ( ji == ni ) THEN
                DO jj = 2, nj-1
                   IF ( mask_coast(ji,jj) == 1 ) THEN
                      X(ji,jj) = 1./( maskv(ji,jj+1)+maskv(ni-1,jj)+maskv(ji,jj-1) + &
                           & rr*maskv(ni-1,jj+1)+rr*maskv(ni-1,jj-1))*( &
                           & maskv(ji,jj+1)*dold(ji,jj+1) + &
                           & maskv(ni-1,jj)*dold(ni-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                           & rr*maskv(ni-1,jj+1)*dold(ni-1,jj+1) + &
                           & rr*maskv(ni-1,jj-1)*dold(ni-1,jj-1) )
                   END IF
                END DO
             END IF
          END IF
       END DO





       !! Center of Top row:
       jj = nj
       DO ji = 2, ni-1
          IF ( mask_coast(ji,jj) == 1 ) THEN
             X(ji,jj) = 1./( maskv(ji+1,jj)+maskv(ji-1,jj)+maskv(ji,jj-1) + &
                  & rr*maskv(ji-1,jj-1)+rr*maskv(ji+1,jj-1) )*( &
                  & maskv(ji+1,jj)*dold(ji+1,jj) + &
                  & maskv(ji-1,jj)*dold(ji-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                  & rr*maskv(ji-1,jj-1)*dold(ji-1,jj-1) + rr*maskv(ji+1,jj-1)*dold(ji+1,jj-1)  )
          END IF
       END DO

       !! West and East corner of top row:
       DO jci = 1, 2

          ji  = ivi(jci)      ! first ji = 1, then ji = ni

          IF (k_ew >= 0) THEN
             jim = vim_per(jci)  ! ji-1
             jip = vip_per(jci)  ! ji+1
             IF ( mask_coast(ji,jj) == 1 ) THEN
                X(ji,jj) = 1./(maskv(jip,jj)+maskv(jim,jj)+maskv(ji,jj-1) + &
                     & rr*maskv(jim,jj-1)+rr*maskv(jip,jj-1))*( &
                     & maskv(jip,jj)*dold(jip,jj) + &
                     & maskv(jim,jj)*dold(jim,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                     & rr*maskv(jim,jj-1)*dold(jim,jj-1) + rr*maskv(jip,jj-1)*dold(jip,jj-1)  )
             END IF

             ! No E-W periodicity:
          ELSE
             IF ( ji == 1 ) THEN
                IF ( mask_coast(ji,jj) == 1 ) THEN
                   X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj-1) + &
                        & rr*maskv(2,jj-1))*( &
                        & maskv(2,jj)*dold(2,jj) + &
                        & maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & rr*maskv(2,jj-1)*dold(2,jj-1)  )
                END IF
             END IF
             IF ( ji == ni ) THEN
                IF ( mask_coast(ji,jj) == 1 ) THEN
                   X(ji,jj) = 1./(maskv(ni-1,jj)+maskv(ji,jj-1) + &
                        & rr*maskv(ni-1,jj-1))*( &
                        & maskv(ni-1,jj)*dold(ni-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & rr*maskv(ni-1,jj-1)*dold(ni-1,jj-1)  )
                END IF
             END IF

          END IF
       END DO




       !! Center of Bottom row:
       jj = 1
       DO ji = 2, ni-1
          IF ( mask_coast(ji,jj) == 1 ) THEN
             X(ji,jj) = 1./(maskv(ji+1,jj)+maskv(ji,jj+1)+maskv(ji-1,jj) + &
                  & rr*maskv(ji+1,jj+1)+rr*maskv(ji-1,jj+1) )*( &
                  & maskv(ji+1,jj)*dold(ji+1,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                  & maskv(ji-1,jj)*dold(ji-1,jj) + &
                  & rr*maskv(ji+1,jj+1)*dold(ji+1,jj+1) + rr*maskv(ji-1,jj+1)*dold(ji-1,jj+1) )
          END IF
       END DO

       !! West and East corner of bottom row:
       DO jci = 1, 2

          ji  = ivi(jci)      ! first ji = 1, then ji = ni

          IF (k_ew >= 0) THEN
             jim = vim_per(jci)  ! ji-1
             jip = vip_per(jci)  ! ji+1
             IF ( mask_coast(ji,jj) == 1 ) THEN
                X(ji,jj) = 1./(maskv(jip,jj)+maskv(ji,jj+1)+maskv(jim,jj) + &
                     & rr*maskv(jip,jj+1)+rr*maskv(jim,jj+1) )*( &
                     & maskv(jip,jj)*dold(jip,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                     & maskv(jim,jj)*dold(jim,jj) + &
                     & rr*maskv(jip,jj+1)*dold(jip,jj+1) + rr*maskv(jim,jj+1)*dold(jim,jj+1) )
             END IF

             !! No E-W periodicity:
          ELSE
             IF ( ji == 1 ) THEN
                IF ( mask_coast(ji,jj) == 1 ) THEN
                   X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj+1) + &
                        & rr*maskv(2,jj+1) )*( &
                        & maskv(2,jj)*dold(2,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & rr*maskv(2,jj+1)*dold(2,jj+1) )
                END IF
             END IF
             IF ( ji == ni ) THEN
                IF ( mask_coast(ji,jj) == 1 ) THEN
                   X(ji,jj) = 1./(maskv(ji,jj+1)+maskv(ni-1,jj) + &
                        & rr*maskv(ni-1,jj+1) )*( &
                        & maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & maskv(ni-1,jj)*dold(ni-1,jj) + &
                        & rr*maskv(ni-1,jj+1)*dold(ni-1,jj+1) )
                END IF
             END IF
          END IF
       END DO





       !! Loosing land for the next iteration:
       !! -----------------------------------
       maskv = maskv + mask_coast
       !! -----------------------------------


    END DO




    !! Time to smooth what's been drowned:

    dold(:,:) = X(:,:)

    DO js = 1, nsmooth_max
                     
       xtmp(:,:) = X(:,:)

       !! Center of the domain:
       X(2:ni-1,2:nj-1) = 0.35*xtmp(2:ni-1,2:nj-1) &
            &           + 0.65*0.25*(xtmp(3:ni,2:nj-1) + xtmp(2:ni-1,3:nj) + xtmp(1:ni-2,2:nj-1) + xtmp(2:ni-1,1:nj-2) )

       !! we can use east-west periodicity:
       IF (k_ew >= 0) THEN
          DO jci = 1, 2
             jim = vim_per(jci)  ! ji-1
             ji  = ivi(jci)      ! first ji = 1, then ji = ni
             jip = vip_per(jci)  ! ji+1

             X(ji,2:nj-1) = 0.35*xtmp(ji,2:nj-1) &
                  &       + 0.65*0.25*(xtmp(jip,2:nj-1) + xtmp(ji,3:nj) + xtmp(jim,2:nj-1) + xtmp(ji,1:nj-2) )

          END DO
       END IF

       ! Important to put original sea-values back on sea-domain at each
       ! iteration so they constrain correct values on coastal values on the
       ! continent during iteration.
       X(:,2:nj-1) = mask(:,2:nj-1)*dold(:,2:nj-1) - (mask(:,2:nj-1) - 1)*X(:,2:nj-1)
       
    END DO

        
    DEALLOCATE ( maskv, mtmp, xtmp, dold, mask_coast )

    
    IF ( ldebug ) PRINT *, 'DROWN: jinc =', jinc

    Xout = X

    DEALLOCATE( X )
    
  END SUBROUTINE DROWN









  SUBROUTINE FILL_EXTRA_BANDS(k_ew, X, Y, DAT,   XP4, YP4, DATP4)
    !!
    !!============================================================================
    !! Extending input arrays with an extraband of two points at north,south,east
    !! and west boundaries.
    !!
    !! The extension is done thanks to Akima's exptrapolation method.
    !!
    !! East-west periodicity of global map is taken into account through 'k_ew' :
    !!
    !!
    !!  k_ew : east-west periodicity on the input file/grid
    !!         k_ew = -1  --> no periodicity
    !!         k_ew >= 0  --> periodicity with overlap of k_ew points
    !!
    !!
    !!                       Author : Laurent BRODEAU, 2007
    !!============================================================================
    !!
    !!
    INTEGER ,                INTENT(in)  :: k_ew
    !!
    REAL(8), DIMENSION(:,:), INTENT(in)  :: X, Y, DAT
    !!
    REAL(8), DIMENSION(:,:), INTENT(out) :: XP4, YP4, DATP4
    !!
    !! Local
    INTEGER :: lx, ly, lxp4, lyp4
    INTEGER :: ji, jj
    !!
    !!
    IF ( (size(X,1) /= size(Y,1)).OR.(size(X,2) /= size(Y,2)).OR. &
         & (size(X,1) /= size(DAT,1)).OR.(size(X,2) /= size(DAT,2))) THEN
       PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : size of input coor. and data do not match!!!'; STOP
    END IF
    !!
    IF ( (size(XP4,1) /= size(YP4,1)).OR.(size(XP4,2) /= size(YP4,2)).OR. &
         & (size(XP4,1) /= size(DATP4,1)).OR.(size(XP4,2) /= size(DATP4,2))) THEN
       PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : size of output coor. and data do not match!!!'; STOP
    END IF
    !!
    !!
    lx = size(X,1)
    ly = size(X,2)
    !!
    lxp4 = size(XP4,1)
    lyp4 = size(XP4,2)
    !!
    IF ( lxp4 /= lx + 4 ) THEN
       PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : target x dim is not ni+4!!!'; STOP
    END IF
    IF ( lyp4 /= ly + 4 ) THEN
       PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : target y dim is not nj+4!!!'; STOP
    END IF
    !!
    !!
    !!   C r e a t i n g   e x t e n d e d   a r r a y s  :
    !!   --------------------------------------------------
    !!
    !! Initialising :
    !! --------------
    XP4   = 0.
    YP4   = 0.
    DATP4 = 0.
    !!
    !! Filling centers :
    !! -----------------
    XP4(3:lxp4-2, 3:lyp4-2)     = X(:,:)
    YP4(3:lxp4-2, 3:lyp4-2)     = Y(:,:)
    DATP4(3:lxp4-2, 3:lyp4-2)   = DAT(:,:)
    !!
    !!
    !! X array :
    !! ---------
    !!
    IF (k_ew /= -1) THEN   ! we can use east-west periodicity of input file to
       !!                   ! fill extra bands :
       XP4( 1     , 3:lyp4-2) = X(lx - 1 - k_ew , :) - 360.
       XP4( 2     , 3:lyp4-2) = X(lx - k_ew     , :) - 360.
       XP4(lxp4   , 3:lyp4-2) = X( 2 + k_ew     , :) + 360.
       XP4(lxp4-1 , 3:lyp4-2) = X( 1 + k_ew     , :) + 360.
       !!
    ELSE
       !!
       !! Left side :
       XP4(2, 3:lyp4-2) = X(2,:) - (X(3,:) - X(1,:))
       XP4(1, 3:lyp4-2) = X(1,:) - (X(3,:) - X(1,:))
       !!
       !! Right side :
       XP4(lxp4-1, 3:lyp4-2) = X(lx-1,:) + X(lx,:) - X(lx-2,:)
       XP4(lxp4  , 3:lyp4-2) = X(lx,:)   + X(lx,:) - X(lx-2,:)
       !!
    END IF
    !!
    !!
    !! Bottom side :
    XP4(:, 2) = XP4(:,4) - (XP4(:,5) - XP4(:,3))
    XP4(:, 1) = XP4(:,3) - (XP4(:,5) - XP4(:,3))
    !!
    !! Top side :
    XP4(:,lyp4-1) = XP4(:,lyp4-3) + XP4(:,lyp4-2) - XP4(:,lyp4-4)
    XP4(:,lyp4)   = XP4(:,lyp4-2) + XP4(:,lyp4-2) - XP4(:,lyp4-4)
    !!
    !!
    !!
    !! Y array :
    !! ---------
    !!
    !! Top side :
    YP4(3:lxp4-2, lyp4-1) = Y(:, ly-1) + Y(:,ly) - Y(:,ly-2)
    YP4(3:lxp4-2, lyp4)   = Y(:, ly)   + Y(:,ly) - Y(:,ly-2)
    !! Bottom side :
    YP4(3:lxp4-2, 2) = Y(:,2) - (Y(:,3) - Y(:,1))
    YP4(3:lxp4-2, 1) = Y(:,1) - (Y(:,3) - Y(:,1))
    !!
    !!
    IF (k_ew /= -1) THEN   ! we can use east-west periodicity
       !!                   ! fill extra bands :
       YP4( 1     , :) = YP4(lx - 1 - k_ew + 2, :)
       YP4( 2     , :) = YP4(lx - k_ew     + 2, :)
       YP4(lxp4   , :) = YP4( 2 + k_ew     + 2, :)
       YP4(lxp4-1 , :) = YP4( 1 + k_ew     + 2, :)
       !!
    ELSE
       !!
       !! Left side :
       YP4(2, :) = YP4(4,:) - (YP4(5,:) - YP4(3,:))
       YP4(1, :) = YP4(3,:) - (YP4(5,:) - YP4(3,:))
       !! Right side :
       YP4(lxp4-1,:) = YP4(lxp4-3,:) + YP4(lxp4-2, :) - YP4(lxp4-4, :)
       YP4(lxp4,:)   = YP4(lxp4-2,:) + YP4(lxp4-2,:)  - YP4(lxp4-4, :)
       !!
    END IF
    !!
    !!
    !! Data array :
    !! ------------
    !!
    IF (k_ew /= -1) THEN   ! we can use east-west periodicity of input file to
       !!                   ! fill extra bands :
       DATP4( 1     , 3:lyp4-2) = DAT(lx - 1 - k_ew , :)
       DATP4( 2     , 3:lyp4-2) = DAT(lx - k_ew     , :)
       DATP4(lxp4   , 3:lyp4-2) = DAT( 2 + k_ew     , :)
       DATP4(lxp4-1 , 3:lyp4-2) = DAT( 1 + k_ew     , :)
       !!
       !!
    ELSE
       !!
       !! Left side :
       DO jj = 3, lyp4-2
          CALL extra_2_east(XP4(lxp4-4,jj),XP4(lxp4-3,jj),XP4(lxp4-2,jj),        &
               &          XP4(lxp4-1,jj),XP4(lxp4,jj),                         &
               &          DATP4(lxp4-4,jj),DATP4(lxp4-3,jj),DATP4(lxp4-2,jj),  &
               &          DATP4(lxp4-1,jj),DATP4(lxp4,jj) )
       END DO
       !!
       !! Right side :
       DO jj = 3, lyp4-2
          CALL extra_2_west(XP4(5,jj),XP4(4,jj),XP4(3,jj),                    &
               &          XP4(2,jj),XP4(1,jj),                               &
               &          DATP4(5,jj),DATP4(4,jj),DATP4(3,jj),               &
               &          DATP4(2,jj),DATP4(1,jj) )
       END DO
       !!
       !!
    END IF
    !!
    !!
    !! Top side :
    DO ji = 1, lxp4
       CALL extra_2_east(YP4(ji,lyp4-4),YP4(ji,lyp4-3),YP4(ji,lyp4-2),        &
            &          YP4(ji,lyp4-1),YP4(ji,lyp4),                         &
            &          DATP4(ji,lyp4-4),DATP4(ji,lyp4-3),DATP4(ji,lyp4-2),  &
            &          DATP4(ji,lyp4-1),DATP4(ji,lyp4) )
    END DO
    !!
    !! Bottom side :
    DO ji = 1, lxp4
       CALL extra_2_west(YP4(ji,5),YP4(ji,4),YP4(ji,3),        &
            &          YP4(ji,2),YP4(ji,1),                    &
            &          DATP4(ji,5),DATP4(ji,4),DATP4(ji,3),    &
            &          DATP4(ji,2),DATP4(ji,1) )
    END DO
    !!
    !!
  END SUBROUTINE FILL_EXTRA_BANDS
  !!
  !!
  !!
  SUBROUTINE extra_2_east(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)
    !!
    !!============================================================================
    !!
    !! Extrapolates 2 extra east (or north) points of a curve with Akima's 1D method
    !!
    !! Input  : x1, x2, x3, x4, x5, y1, y2, y3
    !! Output : y4, y5
    !!
    !!                       Author : Laurent BRODEAU, 2007
    !!============================================================================
    !!
    !!
    REAL(8), INTENT(in)  :: x1, x2, x3, x4, x5, y1, y2, y3
    REAL(8), INTENT(out) :: y4, y5
    !!
    !! Local :
    REAL(8) :: A, B, C, D, ALF, BET
    !!
    !!
    A    = x2 - x1
    B    = x3 - x2
    C    = x4 - x3
    D    = x5 - x4
    !!
    ALF  = y2 - y1
    BET  = y3 - y2
    !!
    IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
       y4 = y3 ; y5 = y3
    ELSE
       y4   = C*(2*BET/B - ALF/A) + y3
       y5   = y4 + y4*D/C + BET*D/B - ALF*D/A - y3*D/C
    END IF
    !!
    !!
  END SUBROUTINE extra_2_east
  !!
  !!
  !!
  !!
  SUBROUTINE extra_2_west(x5, x4, x3, x2, x1, y5, y4, y3, y2, y1)
    !!
    !!============================================================================
    !!
    !! Extrapolates 2 extra west (or south) points of a curve with Akima's 1D method
    !!
    !! Input  : x1, x2, x3, x4, x5, y1, y2, y3
    !! Output : y4, y5
    !!
    !!                       Author : Laurent BRODEAU, 2007
    !!============================================================================
    !!
    !!
    REAL(8), INTENT(in)  :: x1, x2, x3, x4, x5, y5, y4, y3
    REAL(8), INTENT(out) :: y1, y2
    REAL(8) :: A, B, C, D, ALF, BET
    !!
    !! x1 -> x5
    !! x2 -> x4
    !! x3 -> x3
    !! x4 -> x2
    !! x5 -> x1
    !!
    A    = x4 - x5
    B    = x3 - x4
    C    = x2 - x3
    D    = x1 - x2
    !!
    ALF  = y4 - y5
    BET  = y3 - y4
    !!
    IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
       y2 = y3; y1 = y3
    ELSE
       y2   = C*(2*BET/B - ALF/A) + y3
       y1   = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C
    END IF
    !!
    !!
  END SUBROUTINE extra_2_west
  !!
  !!
  FUNCTION TEST_XYZ(rx, ry, rz)
    !!
    !! Testing if 2D coordinates or 1D, and if match shape of data...
    !!
    CHARACTER(len=2) :: TEST_XYZ
    !!
    REAL(8), DIMENSION(:,:), INTENT(in) :: rx, ry
    REAL(4), DIMENSION(:,:), INTENT(in) :: rz
    !!
    INTEGER :: ix1, ix2, iy1, iy2, iz1, iz2
    !!
    ix1 = size(rx,1) ; ix2 = size(rx,2)
    iy1 = size(ry,1) ; iy2 = size(ry,2)
    iz1 = size(rz,1) ; iz2 = size(rz,2)
    !!
    IF ( (ix2 == 1).AND.(iy2 == 1) ) THEN
       !!
       IF ( (ix1 == iz1).AND.(iy1 == iz2) ) THEN
          TEST_XYZ = '1d'
       ELSE
          PRINT *, 'ERROR, mod_drown.F90 = >TEST_XYZ 1 : longitude and latitude array do not match data!'
          PRINT *, ''; STOP
       END IF
       !!
    ELSE
       IF ( (ix1 == iz1).AND.(iy1 == iz1).AND.(ix2 == iz2).AND.(iy2 == iz2) ) THEN
          TEST_XYZ = '2d'
       ELSE
          PRINT *, 'ERROR, mod_drown.F90 = >TEST_XYZ 2 : longitude and latitude array do not match data!'
          PRINT *, ''; STOP
       END IF
    END IF
    !!
    !!
  END FUNCTION TEST_XYZ
  !!
  !!
  !!
END MODULE MOD_DROWN

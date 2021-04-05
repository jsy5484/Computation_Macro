SUBROUTINE linspace(d1,d2,n,grid)
    IMPLICIT NONE

    INTEGER, INTENT(IN)                         :: n
    DOUBLE PRECISION, INTENT(IN)                :: d1, d2
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: grid
    INTEGER :: idx

    grid(1) = d1
    DO idx = 0,n-2
       grid(idx+1) = d1+(DBLE(idx)*(d2-d1))/DBLE(n-1)
    END DO
    grid(n) = d2

END SUBROUTINE



PROGRAM VFI_determ
  IMPLICIT NONE
  external linspace
  !GRID SIZES
  INTEGER, parameter :: nkgrid  = 101
  INTEGER  :: iter = 0
  INTEGER  :: i, j

  !PARAMETERS
  REAL(8), parameter   :: alpha = 0.33333333333
  REAL(8), parameter   :: beta = 0.975
  REAL(8), parameter   :: A = 0.4
  REAL(8), parameter   :: delta  = 0.1
  REAL(8), parameter   :: tol  = 1.0e-7
  REAL(8) :: maxDiff  = 10.0
  REAL(8) :: kmin, kmax, rtrn, kss, yss, css


  REAL(8), DIMENSION(nkgrid) 	:: k_grid, kp_grid, V_old, V_new, k_pol, argvalue, argmax
  REAL(8), DIMENSION(nkgrid, nkgrid) :: temp


  ! Find Steady-State
  kss   = ((1 + beta * (delta - 1))/(alpha * beta * A))**(1 / (alpha - 1));
  yss   = A * kss**alpha + (1 - delta)*kss;
  css   = yss - kss;

  print *, 'Steady States are'
  print *, 'Output     : ', yss
  print *, 'Capital    : ', kss
  print *, 'Consumption: ', css
  print *

  ! Construct Capital GRID
  kmin = 0.5*kss; kmax = 1.5*kss;
  call linspace(kmin, kmax, nkgrid, k_grid);
  call linspace(kmin, kmax, nkgrid, kp_grid);

  DO WHILE (maxDiff > tol)
    iter = iter + 1;

    DO i = 1, nkgrid
        DO j = 1, nkgrid
          rtrn = A*(k_grid(i)**alpha) + (1-delta)*k_grid(i) - kp_grid(j);

          IF (rtrn < 0.0) THEN
            temp(i, j) = -100.0;
          ELSE
            temp(i, j) = log(rtrn) + beta*V_old(j);
          END IF

        END DO
    END DO

    ! Find Max and argmax
    ![argvalue, argmax]= max(temp,[],2);
    argvalue = maxval(temp, dim=2)
    V_new = argvalue;
    argmax = maxloc(temp, dim=2)

     DO i = 1, nkgrid
          k_pol(i) = k_grid(argmax(i));
     END DO

     maxDiff = maxval((abs(V_new - V_old)))
     V_old = V_new;
  END DO

  print *, 'Opitmal Value Function is: '
  print *, V_new
  print *
  print *, 'Opitmal Policy Function is: '
  print *, k_pol



END PROGRAM

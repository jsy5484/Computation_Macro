SUBROUTINE linspace(d1,d2,n,grid)
    ! a Fortran code same as linspace in Matlab
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

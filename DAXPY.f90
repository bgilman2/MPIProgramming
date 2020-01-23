! DAXPY function - adds a scalar multiple of a double precision
! vector to another double precision vector
! performs the following operation:
! y <-- alpha*x + y
! incx and incy specify the increment between two consecutive elements
! of respectively vector x and vector y

! ARGUMENTS 
! n	 	INTEGER (input)
!		Number of elements in the vectors. If n <= 0, these routines
!		return without any computation
!
! alpha 	DOUBLE PRECISION (input)
!			
! x 		DOUBLE PRECISION (input)
!		Array of dimension (n-1)* |incx| + 1. Contains the vector
!		to be scaled before summation
!
! incx		INTEGER (input)
!		Increment between elements of x.
!		If incx = 0, the results will be unpredictable
!
! y		DOUBLE PRECISION, (input and output)
!		array of dimension (n-1)*|incy| + 1
!		Before calling the routine, y contains the vector to be summed
!		After the routine ends, y contains the result of the summation
!
! incy		INTEGER (input)
!		Increment between elements of y.
!		If incy = 0, the results will be unpredictable

! source: http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga8f99d6a644d3396aa32db472e0cfc91c.html#ga8f99d6a644d3396aa32db472e0cfc91c


SUBROUTINE DAXPY(n, alpha, x, incx, y, incy)

     ! Arguments
     INTEGER, intent(in) :: n, incx, incy
     REAL, intent(in) :: alpha
     REAL, intent(in) :: x(*)
     REAL, intent(inout) :: y(*)

     INTEGER m, mp1, ix, iy ! local variables

     IF ((n <= 0) .OR. (alpha == 0.0)) RETURN ! if the size for the vectors passed in are too small or if alpha is 0, then terminate the function
     IF ((incx == 1) .AND. (incy == 1)) THEN ! code for if the increments are the same 
          m = mod(n,4) ! take the modulus of n for loop unrolling
          IF (m /= 0) THEN ! if the loop cannot be unrolled go element by element
               DO i = 1,m ! perform the DAXPY operation with no loop unrolling
                    y(i) = y(i) + alpha * x(i)
               END DO
          END IF

          IF (n < 4) RETURN
          mp1 = m + 1
          DO i = mp1, n, 4 ! perform the operation with loop unrolling (unrolled 4 times)
               y(i) = y(i) + alpha * x(i)
               y(i+1) = y(i+1) + alpha * x(i+1)
               y(i+2) = y(i+2) + alpha * x(i+2)
               y(i+3) = y(i+3) + alpha * x(i+3)
          END DO

     ELSE ! code for if the increments are unequal
          ix = 1 ! indexing variables for if the vectors are not equal in length
          iy = 1
          IF (incx < 0) ix = (-n + 1)*incx + 1 ! if the stride variable for vector x is negative, compute the new indexing variable for accessing the end of the vector first
          IF (incy < 0) iy = (-n + 1)*incy + 1 ! if the stride variable for vector y is negative, compute the new indexing variable for accessing the end of the vector first
          DO i = 1, n ! compute the DAXPY function and compute the new indices for each vector
               y(iy) = y(iy) + alpha * x(ix)
               ix = ix + incx
               iy = iy + incy
          END DO
     END IF
     RETURN ! terminate the function
END SUBROUTINE DAXPY

! Testing program for the above subroutine
program DAXPY_Test
real, dimension(5) :: x
real, dimension(5) :: y
integer :: n, incx, incy
real :: alpha

n = 5
incx = 1
incy = 1
alpha = 2.00
x = (/1.0, 1.0, 1.0, 1.0, 1.0 /)
y = (/1.0, 2.0, 3.0, 4.0, 5.0 /)

! Print the values of Y before calling DAXPY
Print *, "Y before DAXPY:"
do i = 1, 5
     Print *, y(i)
end do

call DAXPY(n, alpha, x, incx, y, incy)

Print *, "Y after DAXPY:"
! Print the values of Y after calling DAXPY
do i = 1, 5
     Print *, y(i)
end do

stop
end
